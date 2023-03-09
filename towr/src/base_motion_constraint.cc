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

#include <towr/constraints/base_motion_constraint.h>
#include <towr/variables/cartesian_dimensions.h>
#include <towr/variables/spline_holder.h>
#include <towr/variables/variable_names.h>

namespace towr {

BaseMotionConstraint::BaseMotionConstraint(double T, double dt,
                                           const SplineHolder &spline_holder)
    : TimeDiscretizationConstraint(T, dt, "baseMotion") {
  base_linear_ = spline_holder.base_linear_;
  phase_durations_ = spline_holder.phase_durations_;
  com_init_ = spline_holder.com_init_;
  spheres_vec_ = spline_holder.spheres_vec_;
  node_ballistic_bounds_.resize(k3D, ifopt::NoBound);
  node_ballistic_bounds_.at(AY) = ifopt::BoundGreaterZero;
  SetRows(BuildConstraintInfos());
}

const std::vector<Sphere> &BaseMotionConstraint::GetSpheresAt(double t) const {
  double dt_between_frames = 1.0;
  size_t index = round(t / dt_between_frames);
  return spheres_vec_.at(index);
}

int BaseMotionConstraint::BuildConstraintInfos() {
  int row_id = 0;
  int n_ticks = dts_.size();
  constraint_infos_.resize(n_ticks);
  for (int i = 0; i < n_ticks; i++) {
    int n_spheres = GetSpheresAt(dts_[i]).size();
    int n_constraints = n_spheres; // ballistic constraints
    bool is_ballistic = false;
    if (n_spheres == 0) {
      n_constraints = k3D;
      is_ballistic = true;
    }
    constraint_infos_[i] = {is_ballistic, row_id, n_constraints};
    row_id += n_constraints;
  }

  return row_id;
}

void BaseMotionConstraint::UpdateConstraintAtInstance(double t, int k,
                                                      VectorXd &g) const {
  auto info = constraint_infos_.at(k);
  if (info.is_ballistic) {
    g.middleRows(GetRow(k, AX), k3D) = base_linear_->GetPoint(t).p();
  } else {
    const auto &spheres = GetSpheresAt(t);
    int row_id = info.row_id;
    for (size_t i = 0; i < spheres.size(); i++, row_id++) {
      Eigen::Vector3d dr =
          (base_linear_->GetPoint(t).p() - spheres.at(i).center);
      g(row_id) = dr.squaredNorm() - pow(spheres.at(i).radius, 2);
    }
  }
}

bool BaseMotionConstraint::IsBallistic(double t) const {
  for (const auto &phase_duration : phase_durations_) {
    if (phase_duration->IsContactPhase(t)) {
      return false;
    }
  }
  return true;
}

void BaseMotionConstraint::UpdateBoundsAtInstance(double t, int k,
                                                  VecBound &bounds) const {
  auto info = constraint_infos_.at(k);
  if (info.is_ballistic) {
    for (size_t dim = 0; dim < node_ballistic_bounds_.size(); ++dim)
      bounds.at(GetRow(k, dim)) = node_ballistic_bounds_.at(dim);
  } else {
    for (int dim = 0; dim < info.n_constraints; ++dim)
      bounds.at(GetRow(k, dim)) = ifopt::BoundSmallerZero;
  }
}

void BaseMotionConstraint::UpdateJacobianAtInstance(double t, int k,
                                                    std::string var_set,
                                                    Jacobian &jac) const {
  if (var_set == id::base_lin_nodes) {
    auto info = constraint_infos_.at(k);
    NodeSpline::Jacobian jac_pos = base_linear_->GetJacobianWrtNodes(t, kPos);
    if (info.is_ballistic) {
      jac.middleRows(GetRow(k, AX), k3D) = jac_pos;
    } else {
      const auto &spheres = GetSpheresAt(t);
      int row_id = info.row_id;
      for (size_t i = 0; i < spheres.size(); i++, row_id++) {
        Eigen::Vector3d twice_dr = 2 * (base_linear_->GetPoint(t).p() -
                                        spheres.at(i).center.cast<double>());
        jac.row(GetRow(k, i)) = twice_dr(0) * jac_pos.row(0) +
                                twice_dr(1) * jac_pos.row(1) +
                                twice_dr(2) * jac_pos.row(2);
      }
    }
  }
}

int BaseMotionConstraint::GetRow(int node, int dim) const {
  return constraint_infos_.at(node).row_id + dim;
}

} /* namespace towr */
