/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef RBDL_SPATIALALGEBRAOPERATORS_H
#define RBDL_SPATIALALGEBRAOPERATORS_H

#include <iostream>
#include <cmath>

namespace RigidBodyDynamics
{
namespace Math
{
template <class T>
inline Matrix3<T> VectorCrossMatrix(const Vector3<T> &vector)
{
  return Matrix3<T>(0., -vector[2], vector[1], vector[2], 0., -vector[0], -vector[1],
                    vector[0], 0.);
}

/** \brief Compact representation for Spatial Inertia. */
template <class T>
struct RBDL_DLLAPI SpatialRigidBodyInertia
{
  SpatialRigidBodyInertia<T>()
    : m(0.), h(Vector3<T>::Zero(3, 1)), Ixx(0.), Iyx(0.), Iyy(0.), Izx(0.), Izy(0.), Izz(0.)
  {
  }
  SpatialRigidBodyInertia<T>(T mass, const Vector3<T> &com_mass, const Matrix3<T> &inertia)
    : m(mass)
    , h(com_mass)
    , Ixx(inertia(0, 0))
    , Iyx(inertia(1, 0))
    , Iyy(inertia(1, 1))
    , Izx(inertia(2, 0))
    , Izy(inertia(2, 1))
    , Izz(inertia(2, 2))
  {
  }
  SpatialRigidBodyInertia<T>(T _m, const Vector3<T> &_h, const T &_Ixx, const T &_Iyx,
                             const T &_Iyy, const T &_Izx, const T &_Izy, const T &_Izz)
    : m(_m), h(_h), Ixx(_Ixx), Iyx(_Iyx), Iyy(_Iyy), Izx(_Izx), Izy(_Izy), Izz(_Izz)
  {
  }

  SpatialVector<T> operator*(const SpatialVector<T> &mv)
  {
    Vector3<T> mv_lower(mv[3], mv[4], mv[5]);

    Vector3<T> res_upper = Vector3<T>(Ixx * mv[0] + Iyx * mv[1] + Izx * mv[2],
                                      Iyx * mv[0] + Iyy * mv[1] + Izy * mv[2],
                                      Izx * mv[0] + Izy * mv[1] + Izz * mv[2]) +
                           h.cross(mv_lower);
    Vector3<T> res_lower = m * mv_lower - h.cross(Vector3<T>(mv[0], mv[1], mv[2]));

    return SpatialVector<T>(res_upper[0], res_upper[1], res_upper[2], res_lower[0],
                            res_lower[1], res_lower[2]);
  }

  template <class C>
  SpatialVector<C> operator*(const SpatialVector<C> &mv)
  {
    Vector3<C> mv_lower(mv[3], mv[4], mv[5]);

    Vector3<C> res_upper = Vector3<C>(Ixx * mv[0] + Iyx * mv[1] + Izx * mv[2],
                                      Iyx * mv[0] + Iyy * mv[1] + Izy * mv[2],
                                      Izx * mv[0] + Izy * mv[1] + Izz * mv[2]) +
                           h.template cast<C>().cross(mv_lower);
    Vector3<C> res_lower =
        m * mv_lower - h.template cast<C>().cross(Vector3<C>(mv[0], mv[1], mv[2]));

    return SpatialVector<C>(res_upper[0], res_upper[1], res_upper[2], res_lower[0],
                            res_lower[1], res_lower[2]);
  }

  SpatialRigidBodyInertia<T> operator+(const SpatialRigidBodyInertia<T> &rbi)
  {
    return SpatialRigidBodyInertia<T>(m + rbi.m, h + rbi.h, Ixx + rbi.Ixx, Iyx + rbi.Iyx,
                                      Iyy + rbi.Iyy, Izx + rbi.Izx, Izy + rbi.Izy,
                                      Izz + rbi.Izz);
  }

  void createFromMatrix(const SpatialMatrix<T> &Ic)
  {
    m = Ic(3, 3);
    h.set(-Ic(1, 5), Ic(0, 5), -Ic(0, 4));
    Ixx = Ic(0, 0);
    Iyx = Ic(1, 0);
    Iyy = Ic(1, 1);
    Izx = Ic(2, 0);
    Izy = Ic(2, 1);
    Izz = Ic(2, 2);
  }

  SpatialMatrix<T> toMatrix() const
  {
    SpatialMatrix<T> result;
    result(0, 0) = Ixx;
    result(0, 1) = Iyx;
    result(0, 2) = Izx;
    result(1, 0) = Iyx;
    result(1, 1) = Iyy;
    result(1, 2) = Izy;
    result(2, 0) = Izx;
    result(2, 1) = Izy;
    result(2, 2) = Izz;

    result.template block<3, 3>(0, 3) = VectorCrossMatrix(h);
    result.template block<3, 3>(3, 0) = -VectorCrossMatrix<T>(h);
    result.template block<3, 3>(3, 3) = Matrix3<T>::Identity(3, 3) * m;

    return result;
  }

  void setSpatialMatrix(SpatialMatrix<T> &mat) const
  {
    mat(0, 0) = Ixx;
    mat(0, 1) = Iyx;
    mat(0, 2) = Izx;
    mat(1, 0) = Iyx;
    mat(1, 1) = Iyy;
    mat(1, 2) = Izy;
    mat(2, 0) = Izx;
    mat(2, 1) = Izy;
    mat(2, 2) = Izz;

    mat(3, 0) = 0.;
    mat(3, 1) = h[2];
    mat(3, 2) = -h[1];
    mat(4, 0) = -h[2];
    mat(4, 1) = 0.;
    mat(4, 2) = h[0];
    mat(5, 0) = h[1];
    mat(5, 1) = -h[0];
    mat(5, 2) = 0.;

    mat(0, 3) = 0.;
    mat(0, 4) = -h[2];
    mat(0, 5) = h[1];
    mat(1, 3) = h[2];
    mat(1, 4) = 0.;
    mat(1, 5) = -h[0];
    mat(2, 3) = -h[1];
    mat(2, 4) = h[0];
    mat(2, 5) = 0.;

    mat(3, 3) = m;
    mat(3, 4) = 0.;
    mat(3, 5) = 0.;
    mat(4, 3) = 0.;
    mat(4, 4) = m;
    mat(4, 5) = 0.;
    mat(5, 3) = 0.;
    mat(5, 4) = 0.;
    mat(5, 5) = m;
  }

  static SpatialRigidBodyInertia<T> createFromMassComInertiaC(T mass, const Vector3<T> &com,
                                                              const Matrix3<T> &inertia_C)
  {
    SpatialRigidBodyInertia<T> result;
    result.m = mass;
    result.h = com * mass;
    Matrix3<T> I = inertia_C + VectorCrossMatrix(com) * VectorCrossMatrix(com).transpose() * mass;
    result.Ixx = I(0, 0);
    result.Iyx = I(1, 0);
    result.Iyy = I(1, 1);
    result.Izx = I(2, 0);
    result.Izy = I(2, 1);
    result.Izz = I(2, 2);
    return result;
  }

  template <class C>
  SpatialRigidBodyInertia<C> cast() const
  {
    SpatialRigidBodyInertia<C> casted;
    casted.m = C(m);
    casted.h = h.template cast<C>();
    casted.Ixx = C(Ixx);
    casted.Iyx = C(Iyx);
    casted.Iyy = C(Iyy);
    casted.Izx = C(Izx);
    casted.Izy = C(Izy);
    casted.Izz = C(Izz);

    return casted;
  }

  /// Mass
  T m;
  /// Coordinates of the center of mass
  Vector3<T> h;
  /// Inertia expressed at the origin
  T Ixx, Iyx, Iyy, Izx, Izy, Izz;
};

/** \brief Compact representation of spatial transformations.
 *
 * Instead of using a verbose 6x6 matrix, this structure only stores a 3x3
 * matrix and a 3-d vector to store spatial transformations. It also
 * encapsulates efficient operations such as concatenations and
 * transformation of spatial vectors.
 */
template <class T>
struct RBDL_DLLAPI SpatialTransform
{
  SpatialTransform() : E(Matrix3<T>::Identity(3, 3)), r(Vector3<T>::Zero(3, 1))
  {
  }
  SpatialTransform(const Matrix3<T> &rotation, const Vector3<T> &translation)
    : E(rotation), r(translation)
  {
  }

  /** Same as X * v.
   *
   * \returns (E * w, - E * rxw + E * v)
   */
  SpatialVector<T> apply(const SpatialVector<T> &v_sp)
  {
    Vector3<T> v_rxw(v_sp[3] - r[1] * v_sp[2] + r[2] * v_sp[1],
                     v_sp[4] - r[2] * v_sp[0] + r[0] * v_sp[2],
                     v_sp[5] - r[0] * v_sp[1] + r[1] * v_sp[0]);
    return SpatialVector<T>(E(0, 0) * v_sp[0] + E(0, 1) * v_sp[1] + E(0, 2) * v_sp[2],
                            E(1, 0) * v_sp[0] + E(1, 1) * v_sp[1] + E(1, 2) * v_sp[2],
                            E(2, 0) * v_sp[0] + E(2, 1) * v_sp[1] + E(2, 2) * v_sp[2],
                            E(0, 0) * v_rxw[0] + E(0, 1) * v_rxw[1] + E(0, 2) * v_rxw[2],
                            E(1, 0) * v_rxw[0] + E(1, 1) * v_rxw[1] + E(1, 2) * v_rxw[2],
                            E(2, 0) * v_rxw[0] + E(2, 1) * v_rxw[1] + E(2, 2) * v_rxw[2]);
  }

  /** Same as X^T * f.
   *
   * \returns (E^T * n + rx * E^T * f, E^T * f)
   */
  SpatialVector<T> applyTranspose(const SpatialVector<T> &f_sp)
  {
    Vector3<T> E_T_f(E(0, 0) * f_sp[3] + E(1, 0) * f_sp[4] + E(2, 0) * f_sp[5],
                     E(0, 1) * f_sp[3] + E(1, 1) * f_sp[4] + E(2, 1) * f_sp[5],
                     E(0, 2) * f_sp[3] + E(1, 2) * f_sp[4] + E(2, 2) * f_sp[5]);

    return SpatialVector<T>(E(0, 0) * f_sp[0] + E(1, 0) * f_sp[1] + E(2, 0) * f_sp[2] -
                                r[2] * E_T_f[1] + r[1] * E_T_f[2],
                            E(0, 1) * f_sp[0] + E(1, 1) * f_sp[1] + E(2, 1) * f_sp[2] +
                                r[2] * E_T_f[0] - r[0] * E_T_f[2],
                            E(0, 2) * f_sp[0] + E(1, 2) * f_sp[1] + E(2, 2) * f_sp[2] -
                                r[1] * E_T_f[0] + r[0] * E_T_f[1],
                            E_T_f[0], E_T_f[1], E_T_f[2]);
  }

  /** Same as X^* I X^{-1}
  */
  SpatialRigidBodyInertia<T> apply(const SpatialRigidBodyInertia<T> &rbi)
  {
    return SpatialRigidBodyInertia<T>(
        rbi.m, E * (rbi.h - rbi.m * r),
        E * (Matrix3<T>(rbi.Ixx, rbi.Iyx, rbi.Izx, rbi.Iyx, rbi.Iyy, rbi.Izy, rbi.Izx,
                        rbi.Izy, rbi.Izz) +
             VectorCrossMatrix<T>(r) * VectorCrossMatrix<T>(rbi.h) +
             (VectorCrossMatrix<T>(rbi.h - rbi.m * r) * VectorCrossMatrix<T>(r))) *
            E.transpose());
  }

  /** Same as X^T I X
  */
  SpatialRigidBodyInertia<T> applyTranspose(const SpatialRigidBodyInertia<T> &rbi)
  {
    Vector3<T> E_T_mr = E.transpose() * rbi.h + rbi.m * r;
    return SpatialRigidBodyInertia<T>(
        rbi.m, E_T_mr,
        E.transpose() * Matrix3<T>(rbi.Ixx, rbi.Iyx, rbi.Izx, rbi.Iyx, rbi.Iyy, rbi.Izy,
                                   rbi.Izx, rbi.Izy, rbi.Izz) *
                E -
            VectorCrossMatrix<T>(r) * VectorCrossMatrix<T>(E.transpose() * rbi.h) -
            VectorCrossMatrix<T>(E_T_mr) * VectorCrossMatrix<T>(r));
  }

  SpatialVector<T> applyAdjoint(const SpatialVector<T> &f_sp)
  {
    Vector3<T> En_rxf = E * (Vector3<T>(f_sp[0], f_sp[1], f_sp[2]) -
                             r.cross(Vector3<T>(f_sp[3], f_sp[4], f_sp[5])));
    //		Vector3<T>  En_rxf = E * (Vector3<T>  (f_sp[0], f_sp[1], f_sp[2]) -
    //r.cross(Eigen::Map<Vector3> (&(f_sp[3]))));

    return SpatialVector<T>(En_rxf[0], En_rxf[1], En_rxf[2],
                            E(0, 0) * f_sp[3] + E(0, 1) * f_sp[4] + E(0, 2) * f_sp[5],
                            E(1, 0) * f_sp[3] + E(1, 1) * f_sp[4] + E(1, 2) * f_sp[5],
                            E(2, 0) * f_sp[3] + E(2, 1) * f_sp[4] + E(2, 2) * f_sp[5]);
  }

  SpatialMatrix<T> toMatrix() const
  {
    Matrix3<T> _Erx = E * Matrix3<T>(0., -r[2], r[1], r[2], 0., -r[0], -r[1], r[0], 0.);
    SpatialMatrix<T> result;
    result.template block<3, 3>(0, 0) = E;
    result.template block<3, 3>(0, 3) = Matrix3<T>::Zero(3, 3);
    result.template block<3, 3>(3, 0) = -_Erx;
    result.template block<3, 3>(3, 3) = E;

    return result;
  }

  SpatialMatrix<T> toMatrixAdjoint() const
  {
    Matrix3<T> _Erx = E * Matrix3<T>(0., -r[2], r[1], r[2], 0., -r[0], -r[1], r[0], 0.);
    SpatialMatrix<T> result;
    result.template block<3, 3>(0, 0) = E;
    result.template block<3, 3>(0, 3) = -_Erx;
    result.template block<3, 3>(3, 0) = Matrix3<T>::Zero(3, 3);
    result.template block<3, 3>(3, 3) = E;

    return result;
  }

  SpatialMatrix<T> toMatrixTranspose() const
  {
    Matrix3<T> _Erx = E * Matrix3<T>(0., -r[2], r[1], r[2], 0., -r[0], -r[1], r[0], 0.);
    SpatialMatrix<T> result;
    result.template block<3, 3>(0, 0) = E.transpose();
    result.template block<3, 3>(0, 3) = -_Erx.transpose();
    result.template block<3, 3>(3, 0) = Matrix3<T>::Zero(3, 3);
    result.template block<3, 3>(3, 3) = E.transpose();

    return result;
  }

  SpatialTransform<T> inverse() const
  {
    return SpatialTransform(E.transpose(), -E * r);
  }

  SpatialTransform<T> operator*(const SpatialTransform<T> &XT) const
  {
    return SpatialTransform(E * XT.E, XT.r + XT.E.transpose() * r);
  }

  template <class C>
  SpatialTransform<T> operator*(const SpatialTransform<C> &XT) const
  {
    return SpatialTransform(E * XT.E.template cast<T>(),
                            XT.r.template cast<T>() + XT.E.transpose().template cast<T>() * r);
  }

  void operator*=(const SpatialTransform &XT)
  {
    r = XT.r + XT.E.transpose() * r;
    E *= XT.E;
  }

  template <class C>
  SpatialTransform<C> cast() const
  {
    SpatialTransform<C> res;
    res.E = E.template cast<C>();
    res.r = r.template cast<C>();

    return res;
  }

  Matrix3<T> E;
  Vector3<T> r;
};

template <class T>
inline std::ostream &operator<<(std::ostream &output, const SpatialRigidBodyInertia<T> &rbi)
{
  output << "rbi.m = " << rbi.m << std::endl;
  output << "rbi.h = " << rbi.h.transpose();
  output << "rbi.Ixx = " << rbi.Ixx << std::endl;
  output << "rbi.Iyx = " << rbi.Iyx << " rbi.Iyy = " << rbi.Iyy << std::endl;
  output << "rbi.Izx = " << rbi.Izx << " rbi.Izy = " << rbi.Izy
         << " rbi.Izz = " << rbi.Izz << std::endl;
  return output;
}

template <class T>
inline std::ostream &operator<<(std::ostream &output, const SpatialTransform<T> &X)
{
  output << "X.E = " << std::endl << X.E << std::endl;
  output << "X.r = " << X.r.transpose();
  return output;
}

template <class T>
inline SpatialTransform<T> Xrot(T angle_rad, const Vector3<T> &axis)
{
  T s, c;
  s = sin(angle_rad);
  c = cos(angle_rad);

  return SpatialTransform<T>(Matrix3<T>(axis[0] * axis[0] * (T(1.0) - c) + c,
                                        axis[1] * axis[0] * (T(1.0) - c) + axis[2] * s,
                                        axis[0] * axis[2] * (T(1.0) - c) - axis[1] * s,

                                        axis[0] * axis[1] * (T(1.0) - c) - axis[2] * s,
                                        axis[1] * axis[1] * (T(1.0) - c) + c,
                                        axis[1] * axis[2] * (T(1.0) - c) + axis[0] * s,

                                        axis[0] * axis[2] * (T(1.0) - c) + axis[1] * s,
                                        axis[1] * axis[2] * (T(1.0) - c) - axis[0] * s,
                                        axis[2] * axis[2] * (T(1.0) - c) + c

                                        ),
                             Vector3<T>(T(0.), T(0.), T(0.)));
}

template <class T>
inline SpatialTransform<T> Xrotx(const T &xrot)
{
  T s, c;
  s = sin(xrot);
  c = cos(xrot);
  return SpatialTransform<T>(Matrix3<T>(T(1.), T(0.), T(0.), T(0.), c, s, T(0.), -s, c),
                             Vector3<T>(T(0.), T(0.), T(0.)));
}

template <class T>
inline SpatialTransform<T> Xroty(const T &yrot)
{
  T s, c;
  s = sin(yrot);
  c = cos(yrot);
  return SpatialTransform<T>(Matrix3<T>(c, T(0.), -s, T(0.), T(1.), T(0.), s, T(0.), c),
                             Vector3<T>(T(0.), T(0.), T(0.)));
}

template <class T>
inline SpatialTransform<T> Xrotz(const T &zrot)
{
  T s, c;
  s = sin(zrot);
  c = cos(zrot);
  return SpatialTransform<T>(Matrix3<T>(c, s, T(0.), -s, c, T(0.), T(0.), T(0.), T(1.)),
                             Vector3<T>(T(0.), T(0.), T(0.)));
}

template <class T>
inline SpatialTransform<T> Xtrans(const Vector3<T> &r)
{
  return SpatialTransform<T>(Matrix3<T>::Identity(3, 3), r);
}

template <class T>
inline SpatialMatrix<T> crossm(const SpatialVector<T> &v)
{
  return SpatialMatrix<T>(0, -v[2], v[1], 0, 0, 0, v[2], 0, -v[0], 0, 0, 0, -v[1], v[0],
                          0, 0, 0, 0, 0, -v[5], v[4], 0, -v[2], v[1], v[5], 0, -v[3],
                          v[2], 0, -v[0], -v[4], v[3], 0, -v[1], v[0], 0);
}

template <class T>
inline SpatialVector<T> crossm(const SpatialVector<T> &v1, const SpatialVector<T> &v2)
{
  return SpatialVector<T>(-v1[2] * v2[1] + v1[1] * v2[2], v1[2] * v2[0] - v1[0] * v2[2],
                          -v1[1] * v2[0] + v1[0] * v2[1],
                          -v1[5] * v2[1] + v1[4] * v2[2] - v1[2] * v2[4] + v1[1] * v2[5],
                          v1[5] * v2[0] - v1[3] * v2[2] + v1[2] * v2[3] - v1[0] * v2[5],
                          -v1[4] * v2[0] + v1[3] * v2[1] - v1[1] * v2[3] + v1[0] * v2[4]);
}

template <class T>
inline SpatialMatrix<T> crossf(const SpatialVector<T> &v)
{
  return SpatialMatrix<T>(0, -v[2], v[1], 0, -v[5], v[4], v[2], 0, -v[0], v[5], 0, -v[3],
                          -v[1], v[0], 0, -v[4], v[3], 0, 0, 0, 0, 0, -v[2], v[1], 0, 0,
                          0, v[2], 0, -v[0], 0, 0, 0, -v[1], v[0], 0);
}

template <class T>
inline SpatialVector<T> crossf(const SpatialVector<T> &v1, const SpatialVector<T> &v2)
{
  return SpatialVector<T>(-v1[2] * v2[1] + v1[1] * v2[2] - v1[5] * v2[4] + v1[4] * v2[5],
                          v1[2] * v2[0] - v1[0] * v2[2] + v1[5] * v2[3] - v1[3] * v2[5],
                          -v1[1] * v2[0] + v1[0] * v2[1] - v1[4] * v2[3] + v1[3] * v2[4],
                          -v1[2] * v2[4] + v1[1] * v2[5], v1[2] * v2[3] - v1[0] * v2[5],
                          -v1[1] * v2[3] + v1[0] * v2[4]);
}

} /* Math */

} /* RigidBodyDynamics */

/* RBDL_SPATIALALGEBRAOPERATORS_H*/
#endif
