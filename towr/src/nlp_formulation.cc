/******************************************************************************
Copyright (c) 2018, Alexander W. Winkler. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************/

#include <towr/nlp_formulation.h>

#include <towr/variables/phase_durations.h>
#include <towr/variables/variable_names.h>

#include <towr/constraints/base_motion_constraint.h>
#include <towr/constraints/dynamic_constraint.h>
#include <towr/constraints/force_constraint.h>
#include <towr/constraints/range_of_motion_constraint.h>
#include <towr/constraints/spline_acc_constraint.h>
#include <towr/constraints/swing_constraint.h>
#include <towr/constraints/terrain_constraint.h>
#include <towr/constraints/total_duration_constraint.h>

#include <towr/costs/node_cost.h>
#include <towr/variables/nodes_variables_all.h>

#include <iostream>

namespace towr {

NlpFormulation::NlpFormulation() {}

NlpFormulation::VariablePtrVec NlpFormulation::GetVariableSets(
    SplineHolder &spline_holder,
    const std::vector<PointsOnFrames::Ptr> &ee_motion,
    const ValuesOnFrames::Ptr &torques, const ValuesOnFrames::Ptr &com_init,
    const std::vector<std::vector<Sphere>> &spheres_vec) {
  VariablePtrVec vars;

  auto base_motion = MakeBaseVariables();
  vars.insert(vars.end(), base_motion.begin(), base_motion.end());

  auto ee_force = MakeForceVariables();
  vars.insert(vars.end(), ee_force.begin(), ee_force.end());

  auto contact_schedule = MakeContactScheduleVariables();
  // can also just be fixed timings that aren't optimized over, but still
  // added to spline_holder.

  // stores these readily constructed spline
  spline_holder =
      SplineHolder(base_motion.at(0), ee_motion, // linear
                   params_.GetBasePolyDurations(), ee_force, torques, com_init,
                   spheres_vec, contact_schedule);
  return vars;
}

std::vector<NodesVariables::Ptr> NlpFormulation::MakeBaseVariables() const {
  std::vector<NodesVariables::Ptr> vars;

  int n_nodes = params_.GetBasePolyDurations().size() + 1;

  auto spline_lin =
      std::make_shared<NodesVariablesAll>(n_nodes, k3D, id::base_lin_nodes);

  spline_lin->SetByLinearInterpolation(
      initial_base_.lin.p(), final_base_.lin.p(), params_.GetTotalTime());
  spline_lin->AddStartBound(kPos, {X, Y, Z}, initial_base_.lin.p());
  Eigen::Vector3d zero = Eigen::Vector3d::Zero();
  spline_lin->AddStartBound(kVel, {X, Y, Z}, zero);
  spline_lin->AddFinalBound(kPos, {X, Y, Z}, final_base_.lin.p());
  spline_lin->AddFinalBound(kVel, {X, Y, Z}, zero);
  vars.push_back(spline_lin);

  return vars;
}

std::vector<NodesVariablesPhaseBased::Ptr>
NlpFormulation::MakeForceVariables() const {
  std::vector<NodesVariablesPhaseBased::Ptr> vars;

  double T = params_.GetTotalTime();
  for (int ee = 0; ee < params_.GetEECount(); ee++) {
    auto nodes = std::make_shared<NodesVariablesEEForce>(
        params_.ee_phase_durations_.at(ee),
        params_.ee_in_contact_at_start_.at(ee), id::EEForceNodes(ee),
        params_.duration_force_polynomial_);

    // initialize with mass of robot distributed equally on all legs
    double m = dynamic_model_->m();
    double g = dynamic_model_->g();

    Vector3d f_stance(0.0, m * g / params_.GetEECount(), 0.0);
    nodes->SetByLinearInterpolation(f_stance, f_stance,
                                    T); // stay constant
    vars.push_back(nodes);
  }

  return vars;
}

NlpFormulation::ContraintPtrVec
NlpFormulation::MakeBaseRangeOfMotionConstraint(const SplineHolder &s) const {
  return {std::make_shared<BaseMotionConstraint>(
      params_.GetTotalTime(), params_.dt_constraint_base_motion_, s)};
}

std::vector<PhaseDurations::Ptr>
NlpFormulation::MakeContactScheduleVariables() const {
  std::vector<PhaseDurations::Ptr> vars;

  for (int ee = 0; ee < params_.GetEECount(); ee++) {
    auto var =
        std::make_shared<PhaseDurations>(ee, params_.ee_phase_durations_.at(ee),
                                         params_.ee_in_contact_at_start_.at(ee),
                                         params_.bound_phase_duration_.first,
                                         params_.bound_phase_duration_.second);
    vars.push_back(var);
  }

  return vars;
}

NlpFormulation::ContraintPtrVec
NlpFormulation::GetConstraints(const SplineHolder &spline_holder) const {
  ContraintPtrVec constraints;
  for (auto name : params_.constraints_)
    for (auto c : GetConstraint(name, spline_holder))
      constraints.push_back(c);

  return constraints;
}

NlpFormulation::ContraintPtrVec
NlpFormulation::GetConstraint(Parameters::ConstraintName name,
                              const SplineHolder &s) const {
  switch (name) {
  case Parameters::Dynamic:
    return MakeDynamicConstraint(s);
    //        case Parameters::EndeffectorRom:
    //            return MakeRangeOfMotionBoxConstraint(s);
  case Parameters::BaseRom:
    return MakeBaseRangeOfMotionConstraint(s);
    //        case Parameters::TotalTime:
    //            return MakeTotalTimeConstraint();
    //        case Parameters::Terrain:
    //            return MakeTerrainConstraint();
  case Parameters::Force:
    return MakeForceConstraint();
    //        case Parameters::Swing:
    //            return MakeSwingConstraint();
    //        case Parameters::BaseAcc:
    //            return MakeBaseAccConstraint(s);
  default:
    throw std::runtime_error("constraint not defined!");
  }
}

NlpFormulation::ContraintPtrVec
NlpFormulation::MakeDynamicConstraint(const SplineHolder &s) const {
  auto constraint = std::make_shared<DynamicConstraint>(
      dynamic_model_, params_.GetTotalTime(), params_.dt_constraint_dynamic_,
      params_.torque_threshold_, s);
  return {constraint};
}

NlpFormulation::ContraintPtrVec NlpFormulation::MakeForceConstraint() const {
  ContraintPtrVec constraints;

  for (int ee = 0; ee < params_.GetEECount(); ee++) {
    auto c = std::make_shared<ForceConstraint>(
        terrain_, params_.force_limit_in_normal_direction_, ee);
    constraints.push_back(c);
  }

  return constraints;
}

NlpFormulation::ContraintPtrVec NlpFormulation::GetCosts() const {
  ContraintPtrVec costs;
  for (const auto &pair : params_.costs_)
    for (auto c : GetCost(pair.first, pair.second))
      costs.push_back(c);

  return costs;
}

NlpFormulation::CostPtrVec
NlpFormulation::GetCost(const Parameters::CostName &name, double weight) const {
  switch (name) {
  case Parameters::ForcesCostID:
    return MakeForcesCost(weight);
  default:
    throw std::runtime_error("cost not defined!");
  }
}

NlpFormulation::CostPtrVec NlpFormulation::MakeForcesCost(double weight) const {
  CostPtrVec cost;

  for (int ee = 0; ee < params_.GetEECount(); ee++)
    cost.push_back(
        std::make_shared<NodeCost>(id::EEForceNodes(ee), kPos, Y, weight));

  return cost;
}

} /* namespace towr */
