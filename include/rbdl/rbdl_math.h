/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef RBDL_MATH_H
#define RBDL_MATH_H

#include <rbdl/rbdl_config.h>

#ifdef RBDL_USE_SIMPLE_MATH
#include "rbdl/SimpleMath/SimpleMathFixed.h"
#include "rbdl/SimpleMath/SimpleMathDynamic.h"
#include "rbdl/SimpleMath/SimpleMathMixed.h"
#include "rbdl/SimpleMath/SimpleMathQR.h"
#include "rbdl/SimpleMath/SimpleMathCholesky.h"
#include "rbdl/SimpleMath/SimpleMathCommaInitializer.h"
#include <vector>

typedef SimpleMath::Fixed::Matrix<double, 3, 1> Vector3_t;
typedef SimpleMath::Fixed::Matrix<double, 3, 3> Matrix3_t;
typedef SimpleMath::Fixed::Matrix<double, 4, 1> Vector4_t;

typedef SimpleMath::Fixed::Matrix<double, 6, 1> SpatialVector_t;
typedef SimpleMath::Fixed::Matrix<double, 6, 6> SpatialMatrix_t;

typedef SimpleMath::Fixed::Matrix<double, 6, 3> Matrix63_t;

typedef SimpleMath::Dynamic::Matrix<double> MatrixN_t;
typedef SimpleMath::Dynamic::Matrix<double> VectorN_t;

#else
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/QR>

#include "rbdl/rbdl_eigenmath.h"

template <class T>
using Matrix63_t = Eigen::Matrix<T, 6, 3>;

template <typename T>
using VectorN_t = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template <typename T>
using MatrixN_t = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

#endif

namespace RigidBodyDynamics
{
/** \brief Math types such as vectors and matrices and utility functions. */
namespace Math
{
template <class T>
using Vector3 = Vector3_t<T>;
using Vector3d = Vector3<double>;

template <class T>
using Vector4 = Vector4_t<T>;
using Vector4d = Vector4<double>;

template <class T>
using Matrix3 = Matrix3_t<T>;
using Matrix3d = Matrix3<double>;

template <class T>
using SpatialVector = SpatialVector_t<T>;
using SpatialVectord = SpatialVector<double>;

template <class T>
using SpatialMatrix = SpatialMatrix_t<T>;
using SpatialMatrixd = SpatialMatrix<double>;

template <class T>
using Matrix63 = Matrix63_t<T>;
using Matrix63d = Matrix63<double>;

template <class T>
using VectorN = VectorN_t<T>;
using VectorNd = VectorN<double>;

template <class T>
using MatrixN = MatrixN_t<T>;
using MatrixNd = MatrixN<double>;

template <class T>
using Isometry3 = Isometry3_t<T>;
using Isometry3d = Isometry3<double>;

} /* Math */

} /* RigidBodyDynamics */

#include "rbdl/Quaternion.h"
#include "rbdl/SpatialAlgebraOperators.h"

// If we use Eigen3 we have to create specializations of the STL
// std::vector such that the alignment is done properly.
#ifndef RBDL_USE_SIMPLE_MATH
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RigidBodyDynamics::Math::SpatialVectord)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RigidBodyDynamics::Math::SpatialMatrixd)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RigidBodyDynamics::Math::Matrix63d)

namespace RigidBodyDynamics
{
namespace Math
{
using SpatialTransformd = SpatialTransform<double>;
using SpatialRigidBodyInertiad = SpatialRigidBodyInertia<double>;
}
}

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RigidBodyDynamics::Math::SpatialTransformd)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RigidBodyDynamics::Math::SpatialRigidBodyInertiad)
#endif

/* RBDL_MATH_H_H */
#endif
