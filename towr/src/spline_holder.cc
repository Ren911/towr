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

#include <towr/variables/phase_spline.h>
#include <towr/variables/spline_holder.h>

namespace towr {
SplineHolder::SplineHolder(
    NodesVariables::Ptr base_lin_nodes,
    std::vector<PointsOnFrames::Ptr> ee_motion,
    const std::vector<double> &base_poly_durations,
    std::vector<NodesVariablesPhaseBased::Ptr> ee_force_nodes,
    const ValuesOnFrames::Ptr &torques, const ValuesOnFrames::Ptr &com_init,
    const std::vector<std::vector<Sphere>> &spheres_vec,
    std::vector<PhaseDurations::Ptr> phase_durations) {
  base_linear_ =
      std::make_shared<NodeSpline>(base_lin_nodes.get(), base_poly_durations);
  phase_durations_ = phase_durations;
  ee_motion_ = ee_motion;
  torques_ = torques;
  com_init_ = com_init;
  spheres_vec_ = spheres_vec;
  for (uint ee = 0; ee < ee_force_nodes.size(); ++ee) {
    // spline without changing the polynomial durations
    auto ee_force_poly_durations = ee_force_nodes.at(ee)->GetPolyDurations();
    ee_force_.push_back(std::make_shared<NodeSpline>(
        ee_force_nodes.at(ee).get(), ee_force_poly_durations));
  }
} // namespace towr

PointsOnFrames::PointsOnFrames(
    const std::map<size_t, Eigen::Vector3d> &positions, double dt_frames)
    : positions_(positions), dt_(dt_frames){};
Eigen::Vector3d PointsOnFrames::GetPoint(double t) const {
  size_t index = round(t / dt_);
  Eigen::Vector3d result = Eigen::Vector3d::Zero();
  if (positions_.find(index) != positions_.end()) {
    result = positions_.at(index);
  }
  return result;
}

ValuesOnFrames::ValuesOnFrames(const PhysArray &values, double dt_frames)
    : values_(values), dt_(dt_frames){};
Eigen::Vector3d ValuesOnFrames::GetPoint(double t) const {
  size_t index = round(t / dt_);
  return values_.row(index).cast<double>();
}
} /* namespace towr */
