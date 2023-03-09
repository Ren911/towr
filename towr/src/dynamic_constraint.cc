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

#include <towr/constraints/dynamic_constraint.h>

#include <towr/variables/variable_names.h>
#include <towr/variables/cartesian_dimensions.h>

namespace towr {

DynamicConstraint::DynamicConstraint(const DynamicModel::Ptr& m, double T,
                                     double dt, double torque_threshold,
                                     const SplineHolder& spline_holder)
    : TimeDiscretizationConstraint(T, dt, "dynamic") {
    model_ = m;

    // link with up-to-date spline variables
    base_linear_ = spline_holder.base_linear_;
    ee_forces_ = spline_holder.ee_force_;
    ee_motion_ = spline_holder.ee_motion_;
    torques_ = spline_holder.torques_;
    phase_durations_ = spline_holder.phase_durations_;
    torque_bounds_ = ifopt::Bounds(-torque_threshold, torque_threshold);

  SetRows(GetNumberOfNodes()*k6D);
}

int
DynamicConstraint::GetRow (int k, Dim6D dimension) const
{
  return k6D*k + dimension;
}

void
DynamicConstraint::UpdateConstraintAtInstance(double t, int k, VectorXd& g) const
{
  UpdateModel(t);
  g.segment(GetRow(k,AX), k6D) = model_->GetDynamicViolation();
}

void DynamicConstraint::UpdateBoundsAtInstance(double t, int k,
                                               VecBound& bounds) const {
    Bounds bounds_value;
    if (!IsBallistic(t)) {
        bounds_value = torque_bounds_;
    } else {
        bounds_value = ifopt::NoBound;
    }
    for (auto dim : {AX, AY, AZ}) bounds.at(GetRow(k, dim)) = bounds_value;
    for (auto dim : {LX, LY, LZ}) bounds.at(GetRow(k, dim)) = ifopt::BoundZero;
}

bool DynamicConstraint::IsBallistic(double t) const {
    for (const auto& phase_duration : phase_durations_) {
        if (phase_duration->IsContactPhase(t)) {
            return false;
        }
    }
    return true;
}

void DynamicConstraint::UpdateJacobianAtInstance(double t, int k,
                                                 std::string var_set,
                                                 Jacobian& jac) const {
    UpdateModel(t);

    int n = jac.cols();
    Jacobian jac_model(k6D, n);

    // sensitivity of dynamic constraint w.r.t base variables.
    if (var_set == id::base_lin_nodes) {
        Jacobian jac_base_lin_pos = base_linear_->GetJacobianWrtNodes(t, kPos);
        Jacobian jac_base_lin_acc = base_linear_->GetJacobianWrtNodes(t, kAcc);

        jac_model =
            model_->GetJacobianWrtBaseLin(jac_base_lin_pos, jac_base_lin_acc);
    }

    // sensitivity of dynamic constraint w.r.t. endeffector variables
    for (int ee = 0; ee < model_->GetEECount(); ++ee) {
        if (var_set == id::EEForceNodes(ee)) {
            Jacobian jac_ee_force =
                ee_forces_.at(ee)->GetJacobianWrtNodes(t, kPos);
            jac_model = model_->GetJacobianWrtForce(jac_ee_force, ee);
        }
    }

    jac.middleRows(GetRow(k, AX), k6D) = jac_model;
}

void DynamicConstraint::UpdateModel(double t) const {
    auto com = base_linear_->GetPoint(t);

    int n_ee = model_->GetEECount();
    std::vector<Eigen::Vector3d> ee_pos;
    std::vector<Eigen::Vector3d> ee_force;
    for (int ee = 0; ee < n_ee; ++ee) {
        ee_force.push_back(ee_forces_.at(ee)->GetPoint(t).p());
        ee_pos.push_back(ee_motion_.at(ee)->GetPoint(t));
    }
    Eigen::Vector3d tau_anim = torques_->GetPoint(t);

    model_->SetCurrent(com.p(), com.a(), tau_anim, ee_force, ee_pos);
}

} /* namespace towr */
