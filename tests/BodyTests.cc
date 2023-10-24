#include <UnitTest++.h>

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
TEST ( TestComputeSpatialInertiaFromAbsoluteRadiiGyration ) {
  Body body(1.1, Vector3d (1.5, 1.2, 1.3), Vector3d (1.4, 2., 3.));

  Matrix3d inertia_C (
      1.4, 0., 0.,
      0., 2., 0., 
      0., 0., 3.);

  SpatialMatrix reference_inertia (
      4.843, -1.98, -2.145, 0, -1.43, 1.32,
      -1.98, 6.334, -1.716, 1.43, 0, -1.65,
      -2.145, -1.716, 7.059, -1.32, 1.65, 0,
      0, 1.43, -1.32, 1.1, 0, 0,
      -1.43, 0, 1.65, 0, 1.1, 0,
      1.32, -1.65, 0, 0, 0, 1.1
      );

  //	cout << LogOutput.str() << endl;

  SpatialRigidBodyInertia body_rbi = SpatialRigidBodyInertia::createFromMassComInertiaC (body.mMass, body.mCenterOfMass, body.mInertia);

  CHECK_ARRAY_CLOSE (reference_inertia.data(), body_rbi.toMatrix().data(), 36, TEST_PREC);
  CHECK_ARRAY_CLOSE (inertia_C.data(), body.mInertia.data(), 9, TEST_PREC);
}

TEST ( TestBodyConstructorMassComInertia ) {
  double mass = 1.1;
  Vector3d com (1.5, 1.2, 1.3);
  Matrix3d inertia_C (
      8.286, -3.96, -4.29,
      -3.96, 10.668, -3.432,
      -4.29, -3.432, 11.118
      );

  Body body (mass, com, inertia_C);

  SpatialMatrix reference_inertia (
      11.729, -5.94, -6.435, 0, -1.43, 1.32,
      -5.94, 15.002, -5.148, 1.43, 0, -1.65,
      -6.435, -5.148, 15.177, -1.32, 1.65, 0,
      0, 1.43, -1.32, 1.1, 0, 0,
      -1.43, 0, 1.65, 0, 1.1, 0,
      1.32, -1.65, 0, 0, 0, 1.1
      );

  SpatialRigidBodyInertia body_rbi = SpatialRigidBodyInertia::createFromMassComInertiaC (body.mMass, body.mCenterOfMass, body.mInertia);
  CHECK_ARRAY_CLOSE (reference_inertia.data(), body_rbi.toMatrix().data(), 36, TEST_PREC);
  CHECK_ARRAY_CLOSE (inertia_C.data(), body.mInertia.data(), 9, TEST_PREC);
}

TEST ( TestBodyJoinNullbody ) {
  ClearLogOutput();
  Body body(1.1, Vector3d (1.5, 1.2, 1.3), Vector3d (1.4, 2., 3.));
  Body nullbody (0., Vector3d (0., 0., 0.), Vector3d (0., 0., 0.));

  Body joined_body = body;
  joined_body.Join (Xtrans(Vector3d (0., 0., 0.)), nullbody);

  SpatialRigidBodyInertia body_rbi (body.mMass, body.mCenterOfMass, body.mInertia);
  SpatialRigidBodyInertia joined_body_rbi (joined_body.mMass, joined_body.mCenterOfMass, joined_body.mInertia);

  CHECK_EQUAL (1.1, body.mMass);
  CHECK_ARRAY_CLOSE (body.mCenterOfMass.data(), joined_body.mCenterOfMass.data(), 3, TEST_PREC);
  CHECK_ARRAY_CLOSE (body_rbi.toMatrix().data(), joined_body_rbi.toMatrix().data(), 36, TEST_PREC);
}

TEST ( TestBodyJoinTwoBodies ) {
  ClearLogOutput();
  Body body_a(1.1, Vector3d (-1.1, 1.3, 0.), Vector3d (3.1, 3.2, 3.3));
  Body body_b(1.1, Vector3d (1.1, 1.3, 0.), Vector3d (3.1, 3.2, 3.3));

  Body body_joined (body_a);
  body_joined.Join (Xtrans(Vector3d (0., 0., 0.)), body_b);

  SpatialRigidBodyInertia body_joined_rbi = SpatialRigidBodyInertia::createFromMassComInertiaC (body_joined.mMass, body_joined.mCenterOfMass, body_joined.mInertia);

  SpatialMatrix reference_inertia (
      9.918, 0, 0, 0, -0, 2.86,
      0, 9.062, 0, 0, 0, -0,
      0, 0, 12.98, -2.86, 0, 0,
      0, 0, -2.86, 2.2, 0, 0,
      -0, 0, 0, 0, 2.2, 0,
      2.86, -0, 0, 0, 0, 2.2
      );

  CHECK_EQUAL (2.2, body_joined.mMass);
  CHECK_ARRAY_EQUAL (Vector3d (0., 1.3, 0.).data(), body_joined.mCenterOfMass.data(), 3);
  CHECK_ARRAY_CLOSE (reference_inertia.data(), body_joined_rbi.toMatrix().data(), 36, TEST_PREC);
}

TEST ( TestBodyJoinTwoBodiesDisplaced ) {
  ClearLogOutput();
  Body body_a(1.1, Vector3d (-1.1, 1.3, 0.), Vector3d (3.1, 3.2, 3.3));
  Body body_b(1.1, Vector3d (0., 0., 0.), Vector3d (3.1, 3.2, 3.3));

  Body body_joined (body_a);
  body_joined.Join (Xtrans(Vector3d (1.1, 1.3, 0.)), body_b);

  SpatialRigidBodyInertia body_joined_rbi = SpatialRigidBodyInertia::createFromMassComInertiaC (body_joined.mMass, body_joined.mCenterOfMass, body_joined.mInertia);

  SpatialMatrix reference_inertia (
      9.918, 0, 0, 0, -0, 2.86,
      0, 9.062, 0, 0, 0, -0,
      0, 0, 12.98, -2.86, 0, 0,
      0, 0, -2.86, 2.2, 0, 0,
      -0, 0, 0, 0, 2.2, 0,
      2.86, -0, 0, 0, 0, 2.2
      );

  CHECK_EQUAL (2.2, body_joined.mMass);
  CHECK_ARRAY_EQUAL (Vector3d (0., 1.3, 0.).data(), body_joined.mCenterOfMass.data(), 3);
  CHECK_ARRAY_CLOSE (reference_inertia.data(), body_joined_rbi.toMatrix().data(), 36, TEST_PREC);


}

TEST ( TestBodyJoinTwoBodiesRotated ) {
  ClearLogOutput();
  Body body_a(1.1, Vector3d (0., 0., 0.), Vector3d (3.1, 3.2, 3.3));
  Body body_b(1.1, Vector3d (0., 0., 0.), Vector3d (3.1, 3.3, 3.2));

  Body body_joined (body_a);
  body_joined.Join (Xrotx(-M_PI*0.5), body_b);

  SpatialRigidBodyInertia body_joined_rbi (body_joined.mMass, body_joined.mCenterOfMass, body_joined.mInertia);

  SpatialMatrix reference_inertia (
      6.2, 0., 0., 0., 0., 0.,
      0., 6.4, 0., 0., 0., 0.,
      0., 0., 6.6, 0., 0., 0.,
      0., 0., 0., 2.2, 0., 0.,
      0., 0., 0., 0., 2.2, 0.,
      0., 0., 0., 0., 0., 2.2
      );

  CHECK_EQUAL (2.2, body_joined.mMass);
  CHECK_ARRAY_EQUAL (Vector3d (0., 0., 0.).data(), body_joined.mCenterOfMass.data(), 3);
  CHECK_ARRAY_CLOSE (reference_inertia.data(), body_joined_rbi.toMatrix().data(), 36, TEST_PREC);
}

TEST ( TestBodyJoinTwoBodiesRotatedAndTranslated ) {
  ClearLogOutput();
  Body body_a(1.1, Vector3d (0., 0., 0.), Vector3d (3.1, 3.2, 3.3));
  Body body_b(1.1, Vector3d (-1., 1., 0.), Vector3d (3.2, 3.1, 3.3));

  Body body_joined (body_a);
  body_joined.Join (Xrotz(M_PI*0.5) * Xtrans(Vector3d (1., 1., 0.)), body_b);

  SpatialRigidBodyInertia body_joined_rbi (body_joined.mMass, body_joined.mCenterOfMass, body_joined.mInertia);

  SpatialMatrix reference_inertia (
      6.2, 0., 0., 0., 0., 0.,
      0., 6.4, 0., 0., 0., 0.,
      0., 0., 6.6, 0., 0., 0.,
      0., 0., 0., 2.2, 0., 0.,
      0., 0., 0., 0., 2.2, 0.,
      0., 0., 0., 0., 0., 2.2
      );

  CHECK_EQUAL (2.2, body_joined.mMass);
  CHECK_ARRAY_CLOSE (Vector3d (0., 0., 0.).data(), body_joined.mCenterOfMass.data(), 3, TEST_PREC);
  CHECK_ARRAY_CLOSE (reference_inertia.data(), body_joined_rbi.toMatrix().data(), 36, TEST_PREC);
}

TEST ( TestBodyConstructorSpatialRigidBodyInertiaMultiplyMotion ) {
  Body body(1.1, Vector3d (1.5, 1.2, 1.3), Vector3d (1.4, 2., 3.));

  SpatialRigidBodyInertia rbi = SpatialRigidBodyInertia(
      body.mMass,
      body.mCenterOfMass * body.mMass,
      body.mInertia
      );

  SpatialVector mv (1.1, 1.2, 1.3, 1.4, 1.5, 1.6);
  SpatialVector fv_matrix = rbi.toMatrix() * mv;
  SpatialVector fv_rbi = rbi * mv;

  CHECK_ARRAY_CLOSE (
      fv_matrix.data(),
      fv_rbi.data(),
      6,
      TEST_PREC
      );
}

TEST ( TestBodyConstructorSpatialRigidBodyInertia ) {
  Body body(1.1, Vector3d (1.5, 1.2, 1.3), Vector3d (1.4, 2., 3.));

  SpatialRigidBodyInertia rbi = SpatialRigidBodyInertia(
      body.mMass,
      body.mCenterOfMass * body.mMass,
      body.mInertia
      );
  SpatialMatrix spatial_inertia = rbi.toMatrix();

  CHECK_ARRAY_CLOSE (
      spatial_inertia.data(),
      rbi.toMatrix().data(),
      36,
      TEST_PREC
      );
}

TEST ( TestBodyConstructorCopySpatialRigidBodyInertia ) {
  Body body(1.1, Vector3d (1.5, 1.2, 1.3), Vector3d (1.4, 2., 3.));

  SpatialRigidBodyInertia rbi = SpatialRigidBodyInertia(
      body.mMass,
      body.mCenterOfMass * body.mMass,
      body.mInertia
      );

  SpatialRigidBodyInertia rbi_from_matrix;
  rbi_from_matrix.createFromMatrix (rbi.toMatrix());

  //	cout << "Spatial Inertia = " << endl << spatial_inertia << endl;
  //	cout << "rbi = " << endl << rbi.toMatrix() << endl;
  //	cout << "rbi.m = " << rbi.m << endl;
  //	cout << "rbi.h = " << rbi.h.transpose() << endl;
  //	cout << "rbi.I = " << endl << rbi.I << endl;

  CHECK_ARRAY_CLOSE (
      rbi.toMatrix().data(),
      rbi_from_matrix.toMatrix().data(),
      36,
      TEST_PREC
      );
}

TEST_CASE (__FILE__"_TestBodySeparateBody", "" ) {
  Body body_init(1.1, Vector3d (1.5, 1.2, 1.3), Vector3d (1.4, 2., 3.));
  const double mass_init = body_init.mMass;
  const Math::Matrix3d inertia_init = body_init.mInertia;
  const Math::Vector3d com_init = body_init.mCenterOfMass;

  Body fixedbody_1(0.7, Vector3d (1.1, 0.8, 1.6), Vector3d (1.5, 1.2, 1.3));
  Body fixedbody_2(0.4, Vector3d (0.7, 1.2, 0.9), Vector3d (0.6,0.7,0.2));
  Math::SpatialTransform body1_transform(Matrix3d::Identity(), Vector3d::Zero());
  Math::SpatialTransform body2_transform(Matrix3d::Identity(), Vector3d::Zero());

  body_init.Join(body1_transform, fixedbody_1);

  CHECK(fabs(body_init.mMass - (mass_init + fixedbody_1.mMass) < TEST_PREC));

  const double mass_joined_1 = body_init.mMass;
  const Math::Matrix3d inertia_joined_1 = body_init.mInertia;
  const Math::Vector3d com_joined_1 = body_init.mCenterOfMass;
  body_init.Join(body2_transform, fixedbody_2);

  CHECK(fabs(body_init.mMass - (mass_init + fixedbody_1.mMass + fixedbody_2.mMass) < TEST_PREC));

  body_init.Separate(body2_transform, fixedbody_2);

  CHECK(fabs(body_init.mMass - (mass_init + fixedbody_1.mMass) < TEST_PREC));
  CHECK_THAT(body_init.mInertia, AllCloseMatrix(inertia_joined_1, TEST_PREC, TEST_PREC));
  CHECK_THAT(body_init.mCenterOfMass, AllCloseVector(com_joined_1, TEST_PREC, TEST_PREC));


  body_init.Separate(body1_transform, fixedbody_1);

  CHECK(fabs(body_init.mMass - (mass_init) < TEST_PREC));
  CHECK_THAT(body_init.mInertia, AllCloseMatrix(inertia_init, TEST_PREC, TEST_PREC));
  CHECK_THAT(body_init.mCenterOfMass, AllCloseVector(com_init, TEST_PREC, TEST_PREC));

}


TEST_CASE (__FILE__"_TestBodySeparateBodyWithTransform", "" ) {
  Body body_init(1.1, Vector3d (1.5, 1.2, 1.3), Vector3d (1.4, 2., 3.));
  const double mass_init = body_init.mMass;
  const Math::Matrix3d inertia_init = body_init.mInertia;
  const Math::Vector3d com_init = body_init.mCenterOfMass;

  Body fixedbody_1(0.7, Vector3d (1.1, 0.8, 1.6), Vector3d (1.5, 1.2, 1.3));
  Body fixedbody_2(0.4, Vector3d (0.7, 1.2, 0.9), Vector3d (0.6,0.7,0.2));
  Math::SpatialTransform body1_transform(
      (Xrotx(M_PI * 0.5) * Xroty(-M_PI * 0.7) * Xrotz(M_PI * 0.1)).E, Vector3d(1.1, 0.6, -0.8));
  Math::SpatialTransform body2_transform(
        (Xrotx(-M_PI * 0.6) * Xroty(-M_PI * 0.2) * Xrotz(M_PI * 0.9)).E, Vector3d(-0.3, 0.1, 1.8));

  body_init.Join(body1_transform, fixedbody_1);

  CHECK(fabs(body_init.mMass - (mass_init + fixedbody_1.mMass) < TEST_PREC));

  const double mass_joined_1 = body_init.mMass;
  const Math::Matrix3d inertia_joined_1 = body_init.mInertia;
  const Math::Vector3d com_joined_1 = body_init.mCenterOfMass;
  body_init.Join(body2_transform, fixedbody_2);

  CHECK(fabs(body_init.mMass - (mass_init + fixedbody_1.mMass + fixedbody_2.mMass) < TEST_PREC));

  body_init.Separate(body2_transform, fixedbody_2);

  CHECK(fabs(body_init.mMass - (mass_init + fixedbody_1.mMass) < TEST_PREC));
  CHECK_THAT(body_init.mInertia, AllCloseMatrix(inertia_joined_1, TEST_PREC, TEST_PREC));
  CHECK_THAT(body_init.mCenterOfMass, AllCloseVector(com_joined_1, TEST_PREC, TEST_PREC));


  body_init.Separate(body1_transform, fixedbody_1);

  CHECK(fabs(body_init.mMass - (mass_init) < TEST_PREC));
  CHECK_THAT(body_init.mInertia, AllCloseMatrix(inertia_init, TEST_PREC, TEST_PREC));
  CHECK_THAT(body_init.mCenterOfMass, AllCloseVector(com_init, TEST_PREC, TEST_PREC));

}

TEST_CASE (__FILE__"_TestBodySeparateBodyWithTransformChangeOrder", "" ) {
  Body body_init(1.1, Vector3d (1.5, 1.2, 1.3), Vector3d (1.4, 2., 3.));
  const double mass_init = body_init.mMass;
  const Math::Matrix3d inertia_init = body_init.mInertia;
  const Math::Vector3d com_init = body_init.mCenterOfMass;

  Body fixedbody_1(0.7, Vector3d (1.1, 0.8, 1.6), Vector3d (1.5, 1.2, 1.3));
  Body fixedbody_2(0.4, Vector3d (0.7, 1.2, 0.9), Vector3d (0.6,0.7,0.2));
  Math::SpatialTransform body1_transform(
      (Xrotx(M_PI * 0.5) * Xroty(-M_PI * 0.7) * Xrotz(M_PI * 0.1)).E, Vector3d(1.1, 0.6, -0.8));
  Math::SpatialTransform body2_transform(
        (Xrotx(-M_PI * 0.6) * Xroty(-M_PI * 0.2) * Xrotz(M_PI * 0.9)).E, Vector3d(-0.3, 0.1, 1.8));

  // create a body with only body_init and fixedbody_2 to check the values after joining and separating fixedbody_1
  Body body_init_2(body_init);
  body_init_2.Join(body2_transform, fixedbody_2);
  CHECK(fabs(body_init_2.mMass - (mass_init + fixedbody_2.mMass) < TEST_PREC));


  // Join fixedbody 1 and 2 in body init
  body_init.Join(body1_transform, fixedbody_1);

  CHECK(fabs(body_init.mMass - (mass_init + fixedbody_1.mMass) < TEST_PREC));

  body_init.Join(body2_transform, fixedbody_2);

  CHECK(fabs(body_init.mMass - (mass_init + fixedbody_1.mMass + fixedbody_2.mMass) < TEST_PREC));

  body_init.Separate(body1_transform, fixedbody_1);

  CHECK(fabs(body_init.mMass - (mass_init + fixedbody_2.mMass) < TEST_PREC));
  CHECK_THAT(body_init.mInertia, AllCloseMatrix(body_init_2.mInertia, TEST_PREC, TEST_PREC));
  CHECK_THAT(body_init.mCenterOfMass, AllCloseVector(body_init_2.mCenterOfMass, TEST_PREC, TEST_PREC));


  body_init.Separate(body2_transform, fixedbody_2);

  CHECK(fabs(body_init.mMass - (mass_init) < TEST_PREC));
  CHECK_THAT(body_init.mInertia, AllCloseMatrix(inertia_init, TEST_PREC, TEST_PREC));
  CHECK_THAT(body_init.mCenterOfMass, AllCloseVector(com_init, TEST_PREC, TEST_PREC));
}
