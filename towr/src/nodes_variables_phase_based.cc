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

#include <towr/variables/cartesian_dimensions.h>
#include <towr/variables/nodes_variables_phase_based.h>

#include <iostream>

namespace towr {

// split changing phase in n_polys
std::vector<NodesVariablesPhaseBased::PolyInfo>
BuildPolyInfos(const std::vector<double> &phase_durations,
               bool first_phase_constant, int n_polys_in_changing_phase) {
  using PolyInfo = NodesVariablesPhaseBased::PolyInfo;
  std::vector<PolyInfo> polynomial_info;

  bool phase_constant = first_phase_constant;
  int phase_count = phase_durations.size();
  for (int i = 0; i < phase_count; ++i) {
    if (phase_constant)

      polynomial_info.push_back(PolyInfo(i, 0, 1, true, phase_durations.at(i)));
    else
      for (int j = 0; j < n_polys_in_changing_phase; ++j) {
        double duration = phase_durations.at(i) / n_polys_in_changing_phase;
        polynomial_info.push_back(
            PolyInfo(i, j, n_polys_in_changing_phase, false, duration));
      }

    phase_constant =
        !phase_constant; // constant and non-constant phase alternate
  }

  return polynomial_info;
}

// build poly infos with duration <= dt_max
std::vector<NodesVariablesPhaseBased::PolyInfo>
BuildPolyInfos(const std::vector<double> &phase_durations,
               bool first_phase_constant, double dt_max) {
  using PolyInfo = NodesVariablesPhaseBased::PolyInfo;
  std::vector<PolyInfo> polynomial_info;

  bool phase_constant = first_phase_constant;
  int phase_count = phase_durations.size();
  for (int i = 0; i < phase_count; ++i) {
    double duration = phase_durations.at(i);
    if (phase_constant)
      polynomial_info.push_back(PolyInfo(i, 0, 1, true, duration));
    else {
      double t_left = phase_durations.at(i);
      int n_polys_in_changing_phase = std::ceil(t_left / dt_max);
      for (int j = 0; j < n_polys_in_changing_phase; ++j) {
        double duration = t_left > dt_max ? dt_max : t_left;
        polynomial_info.push_back(
            PolyInfo(i, j, n_polys_in_changing_phase, false, duration));
        t_left -= dt_max;
      }
    }

    phase_constant =
        !phase_constant; // constant and non-constant phase alternate
  }

  return polynomial_info;
}

NodesVariablesPhaseBased::NodesVariablesPhaseBased(
    const std::vector<double> &phase_durations, bool first_phase_constant,
    const std::string &name, int n_polys_in_changing_phase)
    : NodesVariables(name) {
  polynomial_info_ = BuildPolyInfos(phase_durations, first_phase_constant,
                                    n_polys_in_changing_phase);

  n_dim_ = k3D;
  int n_nodes = polynomial_info_.size() + 1;
  nodes_ = std::vector<Node>(n_nodes, Node(n_dim_));
}

NodesVariablesPhaseBased::NodesVariablesPhaseBased(
    const std::vector<double> &phase_durations, bool first_phase_constant,
    const std::string &name, double dt_max)
    : NodesVariables(name) {
  polynomial_info_ =
      BuildPolyInfos(phase_durations, first_phase_constant, dt_max);

  n_dim_ = k3D;
  int n_nodes = polynomial_info_.size() + 1;
  nodes_ = std::vector<Node>(n_nodes, Node(n_dim_));
}

NodesVariablesPhaseBased::VecDurations
NodesVariablesPhaseBased::ConvertPhaseToPolyDurations(
    const VecDurations &phase_durations) const {
  VecDurations poly_durations;

  for (int i = 0; i < GetPolynomialCount(); ++i) {
    auto info = polynomial_info_.at(i);
    poly_durations.push_back(phase_durations.at(info.phase_) /
                             info.n_polys_in_phase_);
  }

  return poly_durations;
}

NodesVariablesPhaseBased::VecDurations
NodesVariablesPhaseBased::GetPolyDurations() const {
  VecDurations poly_durations(polynomial_info_.size());

  for (size_t i = 0; i < polynomial_info_.size(); ++i) {
    poly_durations[i] = polynomial_info_.at(i).duration_;
  }
  return poly_durations;
}

bool NodesVariablesPhaseBased::IsConstantNode(int node_id) const {
  bool is_constant = false;

  // node is considered constant if either left or right polynomial
  // belongs to a constant phase
  for (int poly_id : GetAdjacentPolyIds(node_id))
    if (IsInConstantPhase(poly_id))
      is_constant = true;

  return is_constant;
}

bool NodesVariablesPhaseBased::IsInConstantPhase(int poly_id) const {
  return polynomial_info_.at(poly_id).is_constant_;
}

NodesVariablesPhaseBased::NodeIds
NodesVariablesPhaseBased::GetIndicesOfNonConstantNodes() const {
  NodeIds node_ids;

  for (size_t id = 0; id < GetNodes().size(); ++id)
    if (!IsConstantNode(id))
      node_ids.push_back(id);

  return node_ids;
}

int NodesVariablesPhaseBased::GetPhase(int node_id) const {
  assert(!IsConstantNode(node_id)); // because otherwise it has two phases

  int poly_id = GetAdjacentPolyIds(node_id).front();
  return polynomial_info_.at(poly_id).phase_;
}

int NodesVariablesPhaseBased::GetPolyIDAtStartOfPhase(int phase) const {
  for (size_t i = 0; i < polynomial_info_.size(); ++i)
    if (polynomial_info_.at(i).phase_ == phase)
      return i;
  return 0;
}

Eigen::Vector3d
NodesVariablesPhaseBased::GetValueAtStartOfPhase(int phase) const {
  int node_id = GetNodeIDAtStartOfPhase(phase);
  return GetNodes().at(node_id).p();
}

int NodesVariablesPhaseBased::GetNodeIDAtStartOfPhase(int phase) const {
  int poly_id = GetPolyIDAtStartOfPhase(phase);
  return GetNodeId(poly_id, Side::Start);
}

std::vector<int>
NodesVariablesPhaseBased::GetAdjacentPolyIds(int node_id) const {
  std::vector<int> poly_ids;
  int last_node_id = GetNodes().size() - 1;

  if (node_id == 0)
    poly_ids.push_back(0);
  else if (node_id == last_node_id)
    poly_ids.push_back(last_node_id - 1);
  else {
    poly_ids.push_back(node_id - 1);
    poly_ids.push_back(node_id);
  }

  return poly_ids;
}

NodesVariablesPhaseBased::PolyInfo::PolyInfo(int phase, int poly_id_in_phase,
                                             int num_polys_in_phase,
                                             bool is_constant, double duration)
    : phase_(phase), poly_in_phase_(poly_id_in_phase),
      n_polys_in_phase_(num_polys_in_phase), is_constant_(is_constant),
      duration_(duration) {}

void NodesVariablesPhaseBased::SetNumberOfVariables(int n_variables) {
  bounds_ = VecBound(n_variables, ifopt::NoBound);
  SetRows(n_variables);
}

NodesVariablesEEForce::NodesVariablesEEForce(
    const std::vector<double> &phase_durations, bool is_in_contact_at_start,
    const std::string &name, int n_polys_in_changing_phase)
    : NodesVariablesPhaseBased(
          phase_durations,
          !is_in_contact_at_start, // contact phase for force is non-constant
          name, n_polys_in_changing_phase) {
  index_to_node_value_info_ = GetPhaseBasedEEParameterization();
  SetNumberOfVariables(index_to_node_value_info_.size());
}

NodesVariablesEEForce::NodesVariablesEEForce(
    const std::vector<double> &phase_durations, bool is_in_contact_at_start,
    const std::string &name, double dt_max)
    : NodesVariablesPhaseBased(
          phase_durations,
          !is_in_contact_at_start, // contact phase for force is non-constant
          name, dt_max) {
  index_to_node_value_info_ = GetPhaseBasedEEParameterization();
  SetNumberOfVariables(index_to_node_value_info_.size());
}

NodesVariablesEEForce::OptIndexMap
NodesVariablesEEForce::GetPhaseBasedEEParameterization() {
  OptIndexMap index_map;

  int idx = 0; // index in variables set
  for (size_t id = 0; id < nodes_.size(); ++id) {
    // stance node:
    // forces can be created during stance, so these nodes are optimized
    // over.
    if (!IsConstantNode(id)) {
      for (int dim = 0; dim < GetDim(); ++dim) {
        index_map[idx++].push_back(NodeValueInfo(id, kPos, dim));
        index_map[idx++].push_back(NodeValueInfo(id, kVel, dim));
      }
    }
    // swing node (next one will also be swing, so handle that one too)
    else {
      // forces can't exist during swing phase, so no need to be optimized
      // -> all node values simply set to zero.
      nodes_.at(id).at(kPos).setZero();
      nodes_.at(id + 1).at(kPos).setZero();

      nodes_.at(id).at(kVel).setZero();
      nodes_.at(id + 1).at(kVel).setZero();

      id += 1; // already added next constant node, so skip
    }
  }

  return index_map;
}

} /* namespace towr */
