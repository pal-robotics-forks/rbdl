#include <gtest/gtest.h>
#include <eigen_checks/gtest.h>
#include <iostream>

#include "rbdl/rbdl_mathutils.h"
#include "rbdl/Body.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

const double TEST_PREC = 1.0e-14;

/* Tests whether the spatial inertia matches the one specified by its
 * parameters
 */
TEST(BodyTests, TestComputeSpatialInertiaFromAbsoluteRadiiGyration)
{
  Body body(1.1, Vector3d (1.5, 1.2, 1.3), Vector3d (1.4, 2., 3.));

  Matrix3d inertia_C (
        1.4, 0., 0.,
        0., 2., 0.,
        0., 0., 3.);

  SpatialMatrixd reference_inertia (
        4.843, -1.98, -2.145, 0, -1.43, 1.32,
        -1.98, 6.334, -1.716, 1.43, 0, -1.65,
        -2.145, -1.716, 7.059, -1.32, 1.65, 0,
        0, 1.43, -1.32, 1.1, 0, 0,
        -1.43, 0, 1.65, 0, 1.1, 0,
        1.32, -1.65, 0, 0, 0, 1.1
        );

  //	cout << LogOutput.str() << endl;

  SpatialRigidBodyInertiad body_rbi = SpatialRigidBodyInertiad::createFromMassComInertiaC (body.mMass, body.mCenterOfMass, body.mInertia);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR(reference_inertia, body_rbi.toMatrix(), TEST_PREC));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(inertia_C, body.mInertia, TEST_PREC));
}

TEST (BodyTests, TestBodyConstructorMassComInertia ) {
  double mass = 1.1;
  Vector3d com (1.5, 1.2, 1.3);
  Matrix3d inertia_C (
        8.286, -3.96, -4.29,
        -3.96, 10.668, -3.432,
        -4.29, -3.432, 11.118
        );

  Body body (mass, com, inertia_C);

  SpatialMatrixd reference_inertia (
        11.729, -5.94, -6.435, 0, -1.43, 1.32,
        -5.94, 15.002, -5.148, 1.43, 0, -1.65,
        -6.435, -5.148, 15.177, -1.32, 1.65, 0,
        0, 1.43, -1.32, 1.1, 0, 0,
        -1.43, 0, 1.65, 0, 1.1, 0,
        1.32, -1.65, 0, 0, 0, 1.1
        );

  SpatialRigidBodyInertiad body_rbi = SpatialRigidBodyInertiad::createFromMassComInertiaC (body.mMass, body.mCenterOfMass, body.mInertia);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(reference_inertia, body_rbi.toMatrix(), TEST_PREC));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(inertia_C, body.mInertia, TEST_PREC));
}

TEST (BodyTests, TestBodyJoinNullbody ) {
  ClearLogOutput();
  Body body(1.1, Vector3d (1.5, 1.2, 1.3), Vector3d (1.4, 2., 3.));
  Body nullbody (0., Vector3d (0., 0., 0.), Vector3d (0., 0., 0.));

  Body joined_body = body;
  joined_body.Join (Xtrans(Vector3d (0., 0., 0.)), nullbody);

  SpatialRigidBodyInertiad body_rbi (body.mMass, body.mCenterOfMass, body.mInertia);
  SpatialRigidBodyInertiad joined_body_rbi (joined_body.mMass, joined_body.mCenterOfMass, joined_body.mInertia);

  EXPECT_EQ(1.1, body.mMass);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(body.mCenterOfMass, joined_body.mCenterOfMass, TEST_PREC));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(body_rbi.toMatrix(), joined_body_rbi.toMatrix(), TEST_PREC));
}

TEST (BodyTests, TestBodyJoinTwoBodies ) {
  ClearLogOutput();
  Body body_a(1.1, Vector3d (-1.1, 1.3, 0.), Vector3d (3.1, 3.2, 3.3));
  Body body_b(1.1, Vector3d (1.1, 1.3, 0.), Vector3d (3.1, 3.2, 3.3));

  Body body_joined (body_a);
  body_joined.Join (Xtrans(Vector3d (0., 0., 0.)), body_b);

  SpatialRigidBodyInertiad body_joined_rbi = SpatialRigidBodyInertiad::createFromMassComInertiaC (body_joined.mMass, body_joined.mCenterOfMass, body_joined.mInertia);

  SpatialMatrixd reference_inertia (
        9.918, 0, 0, 0, -0, 2.86,
        0, 9.062, 0, 0, 0, -0,
        0, 0, 12.98, -2.86, 0, 0,
        0, 0, -2.86, 2.2, 0, 0,
        -0, 0, 0, 0, 2.2, 0,
        2.86, -0, 0, 0, 0, 2.2
        );

  EXPECT_EQ (2.2, body_joined.mMass);
  EXPECT_TRUE(EIGEN_MATRIX_EQUAL_DOUBLE(Vector3d (0., 1.3, 0.), body_joined.mCenterOfMass));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(reference_inertia, body_joined_rbi.toMatrix(), TEST_PREC));
}

TEST (BodyTests, TestBodyJoinTwoBodiesDisplaced ) {
  ClearLogOutput();
  Body body_a(1.1, Vector3d (-1.1, 1.3, 0.), Vector3d (3.1, 3.2, 3.3));
  Body body_b(1.1, Vector3d (0., 0., 0.), Vector3d (3.1, 3.2, 3.3));

  Body body_joined (body_a);
  body_joined.Join (Xtrans(Vector3d (1.1, 1.3, 0.)), body_b);

  SpatialRigidBodyInertiad body_joined_rbi = SpatialRigidBodyInertiad::createFromMassComInertiaC (body_joined.mMass, body_joined.mCenterOfMass, body_joined.mInertia);

  SpatialMatrixd reference_inertia (
        9.918, 0, 0, 0, -0, 2.86,
        0, 9.062, 0, 0, 0, -0,
        0, 0, 12.98, -2.86, 0, 0,
        0, 0, -2.86, 2.2, 0, 0,
        -0, 0, 0, 0, 2.2, 0,
        2.86, -0, 0, 0, 0, 2.2
        );

  EXPECT_EQ(2.2, body_joined.mMass);
  EXPECT_TRUE(EIGEN_MATRIX_EQUAL_DOUBLE(Vector3d (0., 1.3, 0.), body_joined.mCenterOfMass));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(reference_inertia, body_joined_rbi.toMatrix(), TEST_PREC));


}

TEST (BodyTests, TestBodyJoinTwoBodiesRotated ) {
  ClearLogOutput();
  Body body_a(1.1, Vector3d (0., 0., 0.), Vector3d (3.1, 3.2, 3.3));
  Body body_b(1.1, Vector3d (0., 0., 0.), Vector3d (3.1, 3.3, 3.2));

  Body body_joined (body_a);
  body_joined.Join (Xrotx(-M_PI*0.5), body_b);

  SpatialRigidBodyInertiad body_joined_rbi (body_joined.mMass, body_joined.mCenterOfMass, body_joined.mInertia);

  SpatialMatrixd reference_inertia (
        6.2, 0., 0., 0., 0., 0.,
        0., 6.4, 0., 0., 0., 0.,
        0., 0., 6.6, 0., 0., 0.,
        0., 0., 0., 2.2, 0., 0.,
        0., 0., 0., 0., 2.2, 0.,
        0., 0., 0., 0., 0., 2.2
        );

  EXPECT_EQ(2.2, body_joined.mMass);
  EXPECT_TRUE(EIGEN_MATRIX_EQUAL_DOUBLE(Vector3d (0., 0., 0.), body_joined.mCenterOfMass));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(reference_inertia, body_joined_rbi.toMatrix(), TEST_PREC));
}

TEST (BodyTests, TestBodyJoinTwoBodiesRotatedAndTranslated ) {
  ClearLogOutput();
  Body body_a(1.1, Vector3d (0., 0., 0.), Vector3d (3.1, 3.2, 3.3));
  Body body_b(1.1, Vector3d (-1., 1., 0.), Vector3d (3.2, 3.1, 3.3));

  Body body_joined (body_a);
  body_joined.Join (Xrotz(M_PI*0.5) * Xtrans(Vector3d (1., 1., 0.)), body_b);

  SpatialRigidBodyInertiad body_joined_rbi (body_joined.mMass, body_joined.mCenterOfMass, body_joined.mInertia);

  SpatialMatrixd reference_inertia (
        6.2, 0., 0., 0., 0., 0.,
        0., 6.4, 0., 0., 0., 0.,
        0., 0., 6.6, 0., 0., 0.,
        0., 0., 0., 2.2, 0., 0.,
        0., 0., 0., 0., 2.2, 0.,
        0., 0., 0., 0., 0., 2.2
        );

  EXPECT_EQ(2.2, body_joined.mMass);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(Vector3d (0., 0., 0.), body_joined.mCenterOfMass, TEST_PREC));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(reference_inertia, body_joined_rbi.toMatrix(), TEST_PREC));
}

TEST (BodyTests, TestBodyConstructorSpatialRigidBodyInertiadMultiplyMotion ) {
  Body body(1.1, Vector3d (1.5, 1.2, 1.3), Vector3d (1.4, 2., 3.));

  SpatialRigidBodyInertiad rbi = SpatialRigidBodyInertiad(
                                  body.mMass,
                                  body.mCenterOfMass * body.mMass,
                                  body.mInertia
                                  );

  SpatialVectord  mv (1.1, 1.2, 1.3, 1.4, 1.5, 1.6);
  SpatialVectord  fv_matrix = rbi.toMatrix() * mv;
  SpatialVectord  fv_rbi = rbi * mv;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR(fv_matrix, fv_rbi, TEST_PREC));
}

TEST (BodyTests, TestBodyConstructorSpatialRigidBodyInertiad ) {
  Body body(1.1, Vector3d (1.5, 1.2, 1.3), Vector3d (1.4, 2., 3.));

  SpatialRigidBodyInertiad rbi = SpatialRigidBodyInertiad(
                                  body.mMass,
                                  body.mCenterOfMass * body.mMass,
                                  body.mInertia
                                  );
  SpatialMatrixd spatial_inertia = rbi.toMatrix();

  EXPECT_TRUE(EIGEN_MATRIX_NEAR(spatial_inertia, rbi.toMatrix(), TEST_PREC));
}

TEST (BodyTests, TestBodyConstructorCopySpatialRigidBodyInertiad ){
  Body body(1.1, Vector3d (1.5, 1.2, 1.3), Vector3d (1.4, 2., 3.));

  SpatialRigidBodyInertiad rbi = SpatialRigidBodyInertiad(
                                  body.mMass,
                                  body.mCenterOfMass * body.mMass,
                                  body.mInertia
                                  );

  SpatialRigidBodyInertiad rbi_from_matrix;
  rbi_from_matrix.createFromMatrix (rbi.toMatrix());

  //	cout << "Spatial Inertia = " << endl << spatial_inertia << endl;
  //	cout << "rbi = " << endl << rbi.toMatrix() << endl;
  //	cout << "rbi.m = " << rbi.m << endl;
  //	cout << "rbi.h = " << rbi.h.transpose() << endl;
  //	cout << "rbi.I = " << endl << rbi.I << endl;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR(rbi.toMatrix(), rbi_from_matrix.toMatrix(), TEST_PREC));
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
