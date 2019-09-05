#include <gtest/gtest.h>
#include <eigen_checks/gtest.h>

#include <iostream>
#include <iomanip>

#include "rbdl/Body.h"
#include "rbdl/rbdl_math.h"
#include "rbdl/rbdl_mathutils.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

const double TEST_PREC = 1.0e-14;

SpatialMatrixd spatial_adjoint(const SpatialMatrixd &m) {
  SpatialMatrixd res (m);
  res.block<3,3>(3,0) = m.block<3,3>(0,3);
  res.block<3,3>(0,3) = m.block<3,3>(3,0);
  return res;
}

SpatialMatrixd spatial_inverse(const SpatialMatrixd &m) {
  SpatialMatrixd res(m);
  res.block<3,3>(0,0) = m.block<3,3>(0,0).transpose();
  res.block<3,3>(3,0) = m.block<3,3>(3,0).transpose();
  res.block<3,3>(0,3) = m.block<3,3>(0,3).transpose();
  res.block<3,3>(3,3) = m.block<3,3>(3,3).transpose();
  return res;
}

Matrix3d get_rotation (const SpatialMatrixd &m) {
  return m.block<3,3>(0,0);
}

Vector3d get_translation (const SpatialMatrixd &m) {
  return Vector3d (-m(4,2), m(3,2), -m(3,1));
}

/// \brief Checks the multiplication of a SpatialMatrixd with a SpatialVector
TEST(SpatialAlgebraTests, TestSpatialMatrixTimesSpatialVector) {
  SpatialMatrixd s_matrix (
      1., 0., 0., 0., 0., 7.,
      0., 2., 0., 0., 8., 0.,
      0., 0., 3., 9., 0., 0.,
      0., 0., 6., 4., 0., 0.,
      0., 5., 0., 0., 5., 0.,
      4., 0., 0., 0., 0., 6.
      );
  SpatialVectord s_vector (
      1., 2., 3., 4., 5., 6.
      );

  SpatialVectord result;
  result = s_matrix * s_vector;

  SpatialVectord test_result (
      43., 44., 45., 34., 35., 40.
      );
  EXPECT_EQ (test_result, result);
}

/// \brief Checks the multiplication of a scalar with a SpatialVector
TEST(SpatialAlgebraTests, TestScalarTimesSpatialVector) {
  SpatialVectord s_vector (
      1., 2., 3., 4., 5., 6.
      );

  SpatialVectord result;
  result = 3. * s_vector;

  SpatialVectord test_result(3., 6., 9., 12., 15., 18.);

  EXPECT_EQ (test_result, result);
}

/// \brief Checks the multiplication of a scalar with a SpatialMatrix
TEST(SpatialAlgebraTests, TestScalarTimesSpatialMatrix) {
  SpatialMatrixd s_matrix (
      1., 0., 0., 0., 0., 7.,
      0., 2., 0., 0., 8., 0.,
      0., 0., 3., 9., 0., 0.,
      0., 0., 6., 4., 0., 0.,
      0., 5., 0., 0., 5., 0.,
      4., 0., 0., 0., 0., 6.
      );

  SpatialMatrixd result;
  result = 3. * s_matrix;

  SpatialMatrixd test_result(
      3., 0., 0., 0., 0., 21.,
      0., 6., 0., 0., 24., 0.,
      0., 0., 9., 27., 0., 0.,
      0., 0., 18., 12., 0., 0.,
      0., 15., 0., 0., 15., 0.,
      12., 0., 0., 0., 0., 18.
      );

  EXPECT_EQ (test_result, result);
}

/// \brief Checks the multiplication of a scalar with a SpatialMatrix
TEST(SpatialAlgebraTests, TestSpatialMatrixTimesSpatialMatrix) {
  SpatialMatrixd s_matrix (
      1., 0., 0., 0., 0., 7.,
      0., 2., 0., 0., 8., 0.,
      0., 0., 3., 9., 0., 0.,
      0., 0., 6., 4., 0., 0.,
      0., 5., 0., 0., 5., 0.,
      4., 0., 0., 0., 0., 6.
      );

  SpatialMatrixd result;
  result = s_matrix * s_matrix;

  SpatialMatrixd test_result(
      29., 0., 0., 0., 0., 49.,
      0., 44., 0., 0., 56., 0.,
      0., 0., 63., 63., 0., 0.,
      0., 0., 42., 70., 0., 0.,
      0., 35., 0., 0., 65., 0.,
      28., 0., 0., 0., 0., 64.
      );

  EXPECT_EQ (test_result, result);
}

/// \brief Checks the adjoint method
//
// This method computes a spatial force transformation from a spatial
// motion transformation and vice versa
TEST(SpatialAlgebraTests, TestSpatialMatrixTransformAdjoint) {
  SpatialMatrixd s_matrix (
      1.,  2.,  3.,  4.,  5.,  6.,
      7.,  8.,  9., 10., 11., 12.,
      13., 14., 15., 16., 17., 18.,
      19., 20., 21., 22., 23., 24.,
      25., 26., 27., 28., 29., 30.,
      31., 32., 33., 34., 35., 36.
      );

  SpatialMatrixd result = spatial_adjoint(s_matrix);

  SpatialMatrixd test_result_matrix (
      1.,  2.,  3., 19., 20., 21.,
      7.,  8.,  9., 25., 26., 27.,
      13., 14., 15., 31., 32., 33.,
      4.,  5.,  6., 22., 23., 24.,
      10., 11., 12., 28., 29., 30.,
      16., 17., 18., 34., 35., 36.);

  EXPECT_EQ (test_result_matrix, result);
}

TEST(SpatialAlgebraTests, TestSpatialMatrixInverse) {
  SpatialMatrixd s_matrix (
      0, 1, 2, 0, 1, 2,
      3, 4, 5, 3, 4, 5,
      6, 7, 8, 6, 7, 8,
      0, 1, 2, 0, 1, 2,
      3, 4, 5, 3, 4, 5,
      6, 7, 8, 6, 7, 8
      );

  SpatialMatrixd test_inv (
      0, 3, 6, 0, 3, 6,
      1, 4, 7, 1, 4, 7,
      2, 5, 8, 2, 5, 8,
      0, 3, 6, 0, 3, 6,
      1, 4, 7, 1, 4, 7,
      2, 5, 8, 2, 5, 8
      );

  EXPECT_EQ (test_inv, spatial_inverse(s_matrix));
}

TEST(SpatialAlgebraTests, TestSpatialMatrixGetRotation) {
  SpatialMatrixd spatial_transform (
      1.,  2.,  3.,  0.,  0.,  0.,
      4.,  5.,  6.,  0.,  0.,  0.,
      7.,  8.,  9.,  0.,  0.,  0.,
      0.,  0.,  0.,  0.,  0.,  0.,
      0.,  0.,  0.,  0.,  0.,  0.,
      0.,  0.,  0.,  0.,  0.,  0.
      );

  //	Matrix3d rotation = spatial_transform.block<3,3>(0,0);
  Matrix3d rotation = get_rotation (spatial_transform);
  Matrix3d test_result (
      1., 2., 3.,
      4., 5., 6.,
      7., 8., 9.
      );

  EXPECT_EQ(test_result, rotation);
}

TEST(SpatialAlgebraTests, TestSpatialMatrixGetTranslation) {
  SpatialMatrixd spatial_transform (
      0.,  0.,  0.,  0.,  0.,  0.,
      0.,  0.,  0.,  0.,  0.,  0.,
      0.,  0.,  0.,  0.,  0.,  0.,
      0., -3.,  2.,  0.,  0.,  0.,
      0.,  0., -1.,  0.,  0.,  0.,
      0.,  0.,  0.,  0.,  0.,  0.
      );

  Vector3d translation = get_translation(spatial_transform);
  Vector3d test_result (
      1., 2., 3.
      );

  EXPECT_EQ( test_result, translation);
}

TEST(SpatialAlgebraTests,TestSpatialVectorCross) {
  SpatialVectord s_vec (1., 2., 3., 4., 5., 6.);

  SpatialMatrixd test_cross (
      0., -3.,  2.,  0.,  0.,  0.,
      3.,  0., -1.,  0.,  0.,  0.,
      -2.,  1.,  0.,  0.,  0.,  0.,
      0., -6.,  5.,  0., -3.,  2.,
      6.,  0., -4.,  3.,  0., -1.,
      -5.,  4.,  0., -2.,  1.,  0.
      );

  SpatialMatrixd s_vec_cross (crossm(s_vec));
  EXPECT_EQ (test_cross, s_vec_cross);

  SpatialMatrixd s_vec_crossf (crossf(s_vec));
  SpatialMatrixd test_crossf = -1. * crossm(s_vec).transpose();

  EXPECT_EQ (test_crossf, s_vec_crossf);
}

TEST(SpatialAlgebraTests, TestSpatialVectorCrossmCrossf) {
  SpatialVectord s_vec (1., 2., 3., 4., 5., 6.);
  SpatialVectord t_vec (9., 8., 7., 6., 5., 4.);

  // by explicitly building the matrices (crossm/f with only one vector)
  SpatialVectord crossm_s_x_t = crossm(s_vec) * t_vec;
  SpatialVectord crossf_s_x_t = crossf(s_vec) * t_vec;

  // by using direct computation that avoids building of the matrix
  SpatialVectord crossm_s_t = crossm(s_vec, t_vec);
  SpatialVectord crossf_s_t = crossf(s_vec, t_vec);

  /*
     cout << crossm_s_x_t << endl;
     cout << "---" << endl;
     cout << crossf_s_x_t << endl;
     cout << "---" << endl;
     cout << crossf_s_t << endl;
     */

  EXPECT_EQ (crossm_s_x_t, crossm_s_t);
  EXPECT_EQ (crossf_s_x_t, crossf_s_t);
}

TEST(SpatialAlgebraTests, TestSpatialTransformApply) {
  Vector3d rot (1.1, 1.2, 1.3);
  Vector3d trans (1.1, 1.2, 1.3);

  SpatialTransformd X_st;
  X_st.r = trans;

  SpatialMatrixd X_66_matrix (SpatialMatrixd::Zero());
  X_66_matrix = Xrotz_mat (rot[2]) * Xroty_mat (rot[1]) * Xrotx_mat (rot[0]) * Xtrans_mat(trans);
  X_st.E = X_66_matrix.block<3,3>(0,0);

  // cout << X_66_matrix << endl;
  // cout << X_st.E << endl;
  // cout << X_st.r.transpose() << endl;

  SpatialVectord v (1.1, 2.1, 3.1, 4.1, 5.1, 6.1);
  SpatialVectord v_66_res = X_66_matrix * v;
  SpatialVectord v_st_res = X_st.apply(v);

  // cout << (v_66_res - v_st_res).transpose() << endl;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (v_66_res, v_st_res,TEST_PREC));
}

TEST(SpatialAlgebraTests, TestSpatialTransformApplyTranspose) {
  Vector3d rot (1.1, 1.2, 1.3);
  Vector3d trans (1.1, 1.2, 1.3);

  SpatialTransformd X_st;
  X_st.r = trans;

  SpatialMatrixd X_66_matrix (SpatialMatrixd::Zero());
  X_66_matrix = Xrotz_mat (rot[2]) * Xroty_mat (rot[1]) * Xrotx_mat (rot[0]) * Xtrans_mat(trans);
  X_st.E = X_66_matrix.block<3,3>(0,0);

  // cout << X_66_matrix << endl;
  // cout << X_st.E << endl;
  // cout << X_st.r.transpose() << endl;

  SpatialVectord v (1.1, 2.1, 3.1, 4.1, 5.1, 6.1);
  SpatialVectord v_66_res = X_66_matrix.transpose() * v;
  SpatialVectord v_st_res = X_st.applyTranspose(v);

  // cout << (v_66_res - v_st_res).transpose() << endl;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (v_66_res, v_st_res, TEST_PREC));
}

TEST(SpatialAlgebraTests, TestSpatialTransformApplyAdjoint) {
  SpatialTransformd X (
      Xrotz (0.5) *
      Xroty (0.9) *
      Xrotx (0.2) *
      Xtrans (Vector3d (1.1, 1.2, 1.3))
      );

  SpatialMatrixd X_adjoint = X.toMatrixAdjoint();

  SpatialVectord f (1.1, 2.1, 4.1, 9.2, 3.3, 0.8);
  SpatialVectord f_apply = X.applyAdjoint(f);
  SpatialVectord f_matrix = X_adjoint * f;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (f_matrix, f_apply, TEST_PREC));
}

TEST(SpatialAlgebraTests, TestSpatialTransformToMatrix) {
  Vector3d rot (1.1, 1.2, 1.3);
  Vector3d trans (1.1, 1.2, 1.3);

  SpatialMatrixd X_matrix (SpatialMatrixd::Zero());
  X_matrix = Xrotz_mat (rot[2]) * Xroty_mat (rot[1]) * Xrotx_mat (rot[0]) * Xtrans_mat(trans);

  SpatialTransformd X_st;
  X_st.E = X_matrix.block<3,3>(0,0);
  X_st.r = trans;

  //	SpatialMatrixd X_diff = X_st.toMatrix() - X_matrix;
  //	cout << "Error: " << endl << X_diff << endl;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (X_matrix, X_st.toMatrix(), TEST_PREC));
}

TEST(SpatialAlgebraTests, TestSpatialTransformToMatrixAdjoint) {
  Vector3d rot (1.1, 1.2, 1.3);
  Vector3d trans (1.1, 1.2, 1.3);

  SpatialMatrixd X_matrix (SpatialMatrixd::Zero());
  X_matrix = Xrotz_mat (rot[2]) * Xroty_mat (rot[1]) * Xrotx_mat (rot[0]) * Xtrans_mat(trans);

  SpatialTransformd X_st;
  X_st.E = X_matrix.block<3,3>(0,0);
  X_st.r = trans;

  //	SpatialMatrixd X_diff = X_st.toMatrixAdjoint() - spatial_adjoint(X_matrix);
  //	cout << "Error: " << endl << X_diff << endl;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (spatial_adjoint(X_matrix), X_st.toMatrixAdjoint(), TEST_PREC));
}

TEST(SpatialAlgebraTests, TestSpatialTransformToMatrixTranspose) {
  Vector3d rot (1.1, 1.2, 1.3);
  Vector3d trans (1.1, 1.2, 1.3);

  SpatialMatrixd X_matrix (SpatialMatrixd::Zero());
  X_matrix = Xrotz_mat (rot[2]) * Xroty_mat (rot[1]) * Xrotx_mat (rot[0]) * Xtrans_mat(trans);

  SpatialTransformd X_st;
  X_st.E = X_matrix.block<3,3>(0,0);
  X_st.r = trans;

  // we have to copy the matrix as it is only transposed via a flag and
  // thus data() does not return the proper data.
  SpatialMatrixd X_matrix_transposed = X_matrix.transpose();
  //	SpatialMatrixd X_diff = X_st.toMatrixTranspose() - X_matrix_transposed;
  //	cout << "Error: " << endl << X_diff << endl;
  //	cout << "X_st: " << endl << X_st.toMatrixTranspose() << endl;
  //	cout << "X: " << endl << X_matrix_transposed() << endl;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (X_matrix_transposed, X_st.toMatrixTranspose(), TEST_PREC));
}

TEST(SpatialAlgebraTests, TestSpatialTransformMultiply) {
  Vector3d rot (1.1, 1.2, 1.3);
  Vector3d trans (1.1, 1.2, 1.3);

  SpatialMatrixd X_matrix_1 (SpatialMatrixd::Zero());
  SpatialMatrixd X_matrix_2 (SpatialMatrixd::Zero());

  X_matrix_1 = Xrotz_mat (rot[2]) * Xroty_mat (rot[1]) * Xrotx_mat (rot[0]) * Xtrans_mat(trans);
  X_matrix_2 = Xrotz_mat (rot[2]) * Xroty_mat (rot[1]) * Xrotx_mat (rot[0]) * Xtrans_mat(trans);

  SpatialTransformd X_st_1;
  SpatialTransformd X_st_2;

  X_st_1.E = X_matrix_1.block<3,3>(0,0);
  X_st_1.r = trans;
  X_st_2.E = X_matrix_2.block<3,3>(0,0);
  X_st_2.r = trans;

  SpatialTransformd X_st_res = X_st_1 * X_st_2;
  SpatialMatrixd X_matrix_res = X_matrix_1 * X_matrix_2;

  //	SpatialMatrixd X_diff = X_st_res.toMatrix() - X_matrix_res;
  //	cout << "Error: " << endl << X_diff << endl;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (X_matrix_res, X_st_res.toMatrix(), TEST_PREC));
}

TEST(SpatialAlgebraTests, TestSpatialTransformMultiplyEqual) {
  Vector3d rot (1.1, 1.2, 1.3);
  Vector3d trans (1.1, 1.2, 1.3);

  SpatialMatrixd X_matrix_1 (SpatialMatrixd::Zero());
  SpatialMatrixd X_matrix_2 (SpatialMatrixd::Zero());

  X_matrix_1 = Xrotz_mat (rot[2]) * Xroty_mat (rot[1]) * Xrotx_mat (rot[0]) * Xtrans_mat(trans);
  X_matrix_2 = Xrotz_mat (rot[2]) * Xroty_mat (rot[1]) * Xrotx_mat (rot[0]) * Xtrans_mat(trans);

  SpatialTransformd X_st_1;
  SpatialTransformd X_st_2;

  X_st_1.E = X_matrix_1.block<3,3>(0,0);
  X_st_1.r = trans;
  X_st_2.E = X_matrix_2.block<3,3>(0,0);
  X_st_2.r = trans;

  SpatialTransformd X_st_res = X_st_1;
  X_st_res *= X_st_2;
  SpatialMatrixd X_matrix_res = X_matrix_1 * X_matrix_2;

  //	SpatialMatrixd X_diff = X_st_res.toMatrix() - X_matrix_res;
  //	cout << "Error: " << endl << X_diff << endl;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (X_matrix_res, X_st_res.toMatrix(), TEST_PREC));
}

TEST(SpatialAlgebraTests, TestXrotAxis) {
  SpatialTransformd X_rotX = Xrotx (M_PI * 0.15);
  SpatialTransformd X_rotX_axis = Xrot (M_PI * 0.15, Vector3d (1., 0., 0.));

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (X_rotX.toMatrix(), X_rotX_axis.toMatrix(), TEST_PREC));

  // all the other axes
  SpatialTransformd X_rotX_90 = Xrotx (M_PI * 0.5);
  SpatialTransformd X_rotX_90_axis = Xrot (M_PI * 0.5, Vector3d (1., 0., 0.));

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (X_rotX_90.toMatrix(), X_rotX_90_axis.toMatrix(), TEST_PREC));

  SpatialTransformd X_rotY_90 = Xroty (M_PI * 0.5);
  SpatialTransformd X_rotY_90_axis = Xrot (M_PI * 0.5, Vector3d (0., 1., 0.));

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (X_rotY_90.toMatrix(), X_rotY_90_axis.toMatrix(), TEST_PREC));

  SpatialTransformd X_rotZ_90 = Xrotz (M_PI * 0.5);
  SpatialTransformd X_rotZ_90_axis = Xrot (M_PI * 0.5, Vector3d (0., 0., 1.));

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (X_rotZ_90.toMatrix(), X_rotZ_90_axis.toMatrix(), TEST_PREC));
}

TEST(TSpatialAlgebraTests, estSpatialTransformApplySpatialRigidBodyInertiaAdd) {
  SpatialRigidBodyInertiad rbi (
      1.1,
      Vector3d (1.2, 1.3, 1.4),
      Matrix3d (
        1.1, 0.5, 0.3,
        0.5, 1.2, 0.4,
        0.3, 0.4, 1.3
        ));

  SpatialMatrixd rbi_matrix_added = rbi.toMatrix() + rbi.toMatrix();
  SpatialRigidBodyInertiad rbi_added = rbi + rbi;

  // cout << "rbi = " << endl << rbi.toMatrix() << endl;
  // cout << "rbi_added = " << endl << rbi_added.toMatrix() << endl;
  // cout << "rbi_matrix_added = " << endl << rbi_matrix_added << endl;
  // cout << "diff = " << endl <<
  //  	rbi_added.toMatrix() - rbi_matrix_added << endl;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (
      rbi_matrix_added,
      rbi_added.toMatrix(),
           TEST_PREC)
      );
}

TEST(SpatialAlgebraTests,TestSpatialTransformApplySpatialRigidBodyInertiaFull) {
  SpatialRigidBodyInertiad rbi (
      1.1,
      Vector3d (1.2, 1.3, 1.4),
      Matrix3d (
        1.1, 0.5, 0.3,
        0.5, 1.2, 0.4,
        0.3, 0.4, 1.3
        ));

  SpatialTransformd X (
      Xrotz (0.5) *
      Xroty (0.9) *
      Xrotx (0.2) *
      Xtrans (Vector3d (1.1, 1.2, 1.3))
      );

  SpatialRigidBodyInertiad rbi_transformed = X.apply (rbi);
  SpatialMatrixd rbi_matrix_transformed = X.toMatrixAdjoint () * rbi.toMatrix() * X.inverse().toMatrix();

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (
      rbi_matrix_transformed,
      rbi_transformed.toMatrix(),
      TEST_PREC)
      );
}

TEST(SpatialAlgebraTests, TestSpatialTransformApplyTransposeSpatialRigidBodyInertiaFull) {
  SpatialRigidBodyInertiad rbi (
      1.1,
      Vector3d (1.2, 1.3, 1.4),
      Matrix3d (
        1.1, 0.5, 0.3,
        0.5, 1.2, 0.4,
        0.3, 0.4, 1.3
        ));

  SpatialTransformd X (
      Xrotz (0.5) *
      Xroty (0.9) *
      Xrotx (0.2) *
      Xtrans (Vector3d (1.1, 1.2, 1.3))
      );

  SpatialRigidBodyInertiad rbi_transformed = X.applyTranspose (rbi);
  SpatialMatrixd rbi_matrix_transformed = X.toMatrixTranspose() * rbi.toMatrix() * X.toMatrix();

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (
      rbi_matrix_transformed,
      rbi_transformed.toMatrix(),
      TEST_PREC)
      );
}

TEST(SpatialAlgebraTests, TestSpatialRigidBodyInertiaCreateFromMatrix) {
  double mass = 1.1;
  Vector3d com (0., 0., 0.);
  Matrix3d inertia (
      1.1, 0.5, 0.3,
      0.5, 1.2, 0.4,
      0.3, 0.4, 1.3
      );
  SpatialRigidBodyInertiad body_rbi(mass, com , inertia);

  SpatialMatrixd spatial_inertia = body_rbi.toMatrix();

  SpatialRigidBodyInertiad rbi;
  rbi.createFromMatrix (spatial_inertia);

  EXPECT_DOUBLE_EQ (mass, rbi.m);
  EXPECT_TRUE(EIGEN_MATRIX_EQUAL_DOUBLE (Vector3d(mass * com), rbi.h));
  Matrix3d rbi_I_matrix (
      rbi.Ixx, rbi.Iyx, rbi.Izx,
      rbi.Iyx, rbi.Iyy, rbi.Izy,
      rbi.Izx, rbi.Izy, rbi.Izz
      );
  EXPECT_TRUE(EIGEN_MATRIX_EQUAL_DOUBLE (inertia, rbi_I_matrix));
}

#ifdef USE_SLOW_SPATIAL_ALGEBRA
TEST(TestSpatialLinSolve) {
  SpatialVectord b (1, 2, 0, 1, 1, 1);
  SpatialMatrixd A (
      1., 2., 3., 0., 0., 0.,
      3., 4., 5., 0., 0., 0.,
      6., 7., 7., 0., 0., 0.,
      0., 0., 0., 1., 0., 0.,
      0., 0., 0., 0., 1., 0.,
      0., 0., 0., 0., 0., 1.
      );

  SpatialVectord x = SpatialLinSolve (A, b);
  SpatialVectord x_test (3.5, -6.5, 3.5, 1, 1, 1);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (x_test, x, 6, TEST_PREC));
}
#endif


int main( int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
