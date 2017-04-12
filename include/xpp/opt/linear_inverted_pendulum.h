/**
 @file    linear_inverted_pendulum.h
 @author  Alexander W. Winkler (winklera@ethz.ch)
 @date    Dec 5, 2016
 @brief   Brief description
 */

#ifndef XPP_XPP_OPT_INCLUDE_XPP_OPT_LINEAR_INVERTED_PENDULUM_H_
#define XPP_XPP_OPT_INCLUDE_XPP_OPT_LINEAR_INVERTED_PENDULUM_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <xpp/cartesian_declarations.h>
#include <xpp/endeffectors.h>

namespace xpp {
namespace opt {

class BaseMotion;

class LinearInvertedPendulum {
public:
  using JacobianRow = Eigen::SparseVector<double, Eigen::RowMajor>;

  using ComPos = Eigen::Vector2d;
//  using ComVel = Eigen::Vector2d;
  using ComAcc = Eigen::Vector2d;
  using Cop    = Eigen::Vector2d;

  using EELoad = Endeffectors<double>;
  using EEPos  = EndeffectorsPos;

  LinearInvertedPendulum ();
  virtual ~LinearInvertedPendulum ();

  void SetCurrent(const ComPos&, double height, const EELoad&, const EEPos&);

  ComAcc GetAcceleration() const;

  /** Approximates the acceleration with small angle assumption and calculates
    * Jacobian w.r.t. spline coefficients.
    */
  JacobianRow GetJacobianApproxWrtSplineCoeff(const BaseMotion&, double t_global,
                                              Coords3D dim) const;

  // zmp_ remove this
  double GetDerivativeOfAccWrtCop(d2::Coords dim) const;


  double GetDerivativeOfAccWrtLoad(EndeffectorID, d2::Coords dim) const;
  double GetDerivativeOfAccWrtEEPos(EndeffectorID) const; // same for x and y direction

private:
  ComPos pos_;
//  ComVel vel_;
  double h_;

  Cop CalculateCop(const EELoad&, const EEPos&) const;

  EELoad ee_load_;
  EEPos ee_pos_;
};

} /* namespace opt */
} /* namespace xpp */

#endif /* XPP_XPP_OPT_INCLUDE_XPP_OPT_LINEAR_INVERTED_PENDULUM_H_ */
