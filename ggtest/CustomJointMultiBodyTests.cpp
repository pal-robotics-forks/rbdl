
/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 * Copyright (c) 2016 Matthew Millard <matthew.millard@iwr.uni-heidelberg.de>
 */


#include <gtest/gtest.h>
#include <eigen_checks/gtest.h>

#include <iostream>

#include "FixturesTests.h"
#include "Human36FixtureGTest.h"
#include "rbdl/rbdl_mathutils.h"
#include "rbdl/Logging.h"

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"
#include "rbdl/Dynamics.h"
#include "rbdl/Contacts.h"
#include <vector>

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

const double TEST_PREC = 1.0e-11;
const int NUMBER_OF_MODELS = 3;
const int NUMBER_OF_BODIES = 3;

//==============================================================================
/*

  The purpose of this test is to test that all of the code in RBDL
  related to a multibody mechanism that includes a custom joint functions.
  Specifically this test is for the multi-pass algorithms in rbdl ... namely
  the CompositeRigidBodyAlgorithm. However, because these tests have already
  been written for CustomJointSingleBodyTests.cc, we'll run them all on the
  multibody models that we will be testing.

  We will be testing 3 models to get good coverage of the
  CompositeRigidBodyAlgorithm:

  1. Rx   - multidof - custom
  2. Rx   - custom   - multidof
  3. custom - multidof - Rx

  As before, to test that the model works, we will create a model using
  standard RBDL versions (the reference model), and then we will create
  a model using a custom joint in the correct location. The following
  algorithms will be tested in this manner:

  UpdateKinematicsCustom
  Jacobians
  InverseDynamics
  CompositeRigidBodyAlgorithm
  ForwardDynamics
  CalcMInvTimestau
  ForwardDynamicsContactsKokkevis

*/
//==============================================================================

struct CustomJointTypeRevoluteX : public CustomJoint {
  CustomJointTypeRevoluteX(){
    mDoFCount = 1;
    S = MatrixNd::Zero(6,1);
    S(0,0)=1.0;
    d_u = MatrixNd::Zero(mDoFCount,1);
  }

  virtual void jcalc (const Model &model,
                      ModelData &model_data,
                      unsigned int joint_id,
                      const Math::VectorNd &q,
                      const Math::VectorNd &qdot)
  {
    model_data.X_J[joint_id] = Xrotx(q[model.mJoints[joint_id].q_index]);
    model_data.v_J[joint_id][0] = qdot[model.mJoints[joint_id].q_index];
  }

  virtual void jcalc_X_lambda_S ( const Model &model,
                                  ModelData &model_data,
                                  unsigned int joint_id,
                                  const Math::VectorNd &q)
  {
    model_data.X_lambda[joint_id] =
        Xrotx (q[model.mJoints[joint_id].q_index])
        * model.X_T[joint_id];


    const Joint &joint = model.mJoints[joint_id];
    model.mCustomJoints[joint.custom_joint_index]->S = S;

  }

};

struct CustomEulerZYXJoint : public CustomJoint {
  CustomEulerZYXJoint () {
    mDoFCount = 3;
    S = MatrixNd::Zero (6,3);
    d_u = MatrixNd::Zero(mDoFCount,1);
  }

  virtual void jcalc (Model &model,
                      ModelData &model_data,
                      unsigned int joint_id,
                      const Math::VectorNd &q,
                      const Math::VectorNd &qdot)
  {
    double q0 = q[model.mJoints[joint_id].q_index];
    double q1 = q[model.mJoints[joint_id].q_index + 1];
    double q2 = q[model.mJoints[joint_id].q_index + 2];

    double s0 = sin (q0);
    double c0 = cos (q0);
    double s1 = sin (q1);
    double c1 = cos (q1);
    double s2 = sin (q2);
    double c2 = cos (q2);

    model_data.X_J[joint_id].E = Matrix3d(
                       c0 * c1,                s0 * c1,     -s1,
        c0 * s1 * s2 - s0 * c2, s0 * s1 * s2 + c0 * c2, c1 * s2,
        c0 * s1 * c2 + s0 * s2, s0 * s1 * c2 - c0 * s2, c1 * c2
        );

    S.setZero();
    S(0,0) = -s1;
    S(0,2) = 1.;

    S(1,0) = c1 * s2;
    S(1,1) = c2;

    S(2,0) = c1 * c2;
    S(2,1) = - s2;

    double qdot0 = qdot[model.mJoints[joint_id].q_index];
    double qdot1 = qdot[model.mJoints[joint_id].q_index + 1];
    double qdot2 = qdot[model.mJoints[joint_id].q_index + 2];

    model_data.v_J[joint_id] = S * Vector3d (qdot0, qdot1, qdot2);

    model_data.c_J[joint_id].set(
        -c1*qdot0*qdot1,
        -s1*s2*qdot0*qdot1 + c1*c2*qdot0*qdot2 - s2*qdot1*qdot2,
        -s1*c2*qdot0*qdot1 - c1*s2*qdot0*qdot2 - c2*qdot1*qdot2,
        0., 0., 0.
        );
  }

  virtual void jcalc_X_lambda_S ( Model &model,
                                  ModelData &model_data,
                                  unsigned int joint_id,
                                  const Math::VectorNd &q)
  {


      double q0 = q[model.mJoints[joint_id].q_index];
      double q1 = q[model.mJoints[joint_id].q_index + 1];
      double q2 = q[model.mJoints[joint_id].q_index + 2];

      double s0 = sin (q0);
      double c0 = cos (q0);
      double s1 = sin (q1);
      double c1 = cos (q1);
      double s2 = sin (q2);
      double c2 = cos (q2);


      model_data.X_lambda[joint_id] = SpatialTransformd (
          Matrix3d(
                           c0 * c1,                s0 * c1,     -s1,
            c0 * s1 * s2 - s0 * c2, s0 * s1 * s2 + c0 * c2, c1 * s2,
            c0 * s1 * c2 + s0 * s2, s0 * s1 * c2 - c0 * s2, c1 * c2
            ),
          Vector3d (0., 0., 0.)) * model.X_T[joint_id];

      S.setZero();
      S(0,0) = -s1;
      S(0,2) = 1.;

      S(1,0) = c1 * s2;
      S(1,1) = c2;

      S(2,0) = c1 * c2;
      S(2,1) = - s2;

      const Joint &joint = model.mJoints[joint_id];
      model.mCustomJoints[joint.custom_joint_index]->S = S;
  }
};

//==============================================================================
//Test Fixture
//==============================================================================


class CustomJointMultiBodyFixture : public ::testing::Test{
protected:
  virtual void SetUp () {

    reference_model.resize(NUMBER_OF_MODELS);
    custom_model.resize(NUMBER_OF_MODELS);

    body.resize(NUMBER_OF_MODELS);
    custom_joint.resize(NUMBER_OF_MODELS);

    reference_body_id.resize(NUMBER_OF_MODELS);
    custom_body_id.resize(NUMBER_OF_MODELS);

    for(int i=0; i < NUMBER_OF_MODELS; ++i){
      body.at(i).resize(3);
      custom_joint.at(i).resize(1);

      reference_body_id.at(i).resize(3);
      custom_body_id.at(i).resize(3);

    }

    q.resize(NUMBER_OF_MODELS);
    qdot.resize(NUMBER_OF_MODELS);
    qddot.resize(NUMBER_OF_MODELS);
    tau.resize(NUMBER_OF_MODELS);

    //========================================================
    //Test Model 1: Rx - multidof3 - custom
    //========================================================

    custom_rx_joint1 = CustomJointTypeRevoluteX();

    Matrix3d inertia1 = Matrix3d::Identity(3,3);

    Body body11 = Body (1., Vector3d (1.1, 1.2, 1.3), inertia1);
    Body body12 = Body (2., Vector3d (2.1, 2.2, 2.3), inertia1);
    Body body13 = Body (3., Vector3d (3.1, 3.2, 3.3), inertia1);

    ModelData reference1_data;
    Model reference1(reference1_data);

    ModelData custom1_data;
    Model custom1(custom1_data);

    Vector3d r1 = Vector3d(0.78,-0.125,0.37);

    double th1 = M_PI/6.0;
    double sinTh1 = sin(th1);
    double cosTh1 = cos(th1);

    Matrix3d rm1  = Matrix3d( 1.0,  0.,   0.,
                  0., cosTh1, -sinTh1,
                  0., sinTh1,  cosTh1);

    Vector3d r2 = Vector3d(-0.178,0.2125,-0.937);

    double th2 = M_PI/2.15;
    double sinTh2 = sin(th2);
    double cosTh2 = cos(th2);

    Matrix3d rm2  = Matrix3d(  cosTh2,  0., sinTh2,
                     0.,  1.,   0.,
                  -sinTh2,  0., cosTh2);



    unsigned int reference_body_id10 =
       reference1.AddBody (reference1_data, 0,
                SpatialTransformd(),
                Joint(JointTypeRevoluteX),
                body11);

    unsigned int reference_body_id11 =
       reference1.AddBody (reference1_data, reference_body_id10,
                SpatialTransformd(rm1,r1),
                Joint(JointTypeEulerZYX),
                body12);

    unsigned int reference_body_id12 =
       reference1.AddBody (reference1_data, reference_body_id11,
                SpatialTransformd(rm2,r2),
                Joint(JointTypeRevoluteX),
                body13);


    unsigned int custom_body_id10 =
       custom1.AddBody (custom1_data,   0,
                SpatialTransformd(),
                Joint(JointTypeRevoluteX),
                body11);

    unsigned int custom_body_id11 =
       custom1.AddBody (custom1_data,   custom_body_id10,
                SpatialTransformd(rm1,r1),
                Joint(JointTypeEulerZYX),
                body12);

    unsigned int custom_body_id12 =
       custom1.AddBodyCustomJoint (custom1_data,  custom_body_id11,
                SpatialTransformd(rm2,r2),
                &custom_rx_joint1,
                body13);

    VectorNd q1         = VectorNd::Zero (reference1.q_size);
    VectorNd qdot1      = VectorNd::Zero (reference1.qdot_size);
    VectorNd qddot1     = VectorNd::Zero (reference1.qdot_size);
    VectorNd tau1       = VectorNd::Zero (reference1.qdot_size);

    int idx = 0;
    reference_model.at(idx) = (reference1);
    custom_model.at(idx) = (custom1);

    reference_body_id.at(idx).at(0) = (reference_body_id10);
    reference_body_id.at(idx).at(1) = (reference_body_id11);
    reference_body_id.at(idx).at(2) = (reference_body_id12);

    custom_body_id.at(idx).at(0) = (custom_body_id10);
    custom_body_id.at(idx).at(1) = (custom_body_id11);
    custom_body_id.at(idx).at(2) = (custom_body_id12);

    body.at(idx).at(0) = (body11);
    body.at(idx).at(1) = (body12);
    body.at(idx).at(2) = (body13);
    custom_joint.at(idx).at(0) = (&custom_rx_joint1);

    q.at(idx)     = (q1);
    qdot.at(idx)  = (qdot1);
    qddot.at(idx)   = (qddot1);
    tau.at(idx)   = (tau1);



    //========================================================
    //Test Model 2: Rx - custom - multidof3
    //========================================================

    ModelData reference2_data;
    Model reference2(reference2_data);

    ModelData custom2_data;
    ModelData custom2(custom2_data);

    unsigned int reference_body_id20 =
       reference2.AddBody (0,
                SpatialTransformd(),
                Joint(JointTypeRevoluteX),
                body11);

    unsigned int reference_body_id21 =
       reference2.AddBody (reference_body_id20,
                SpatialTransformd(rm2,r2),
                Joint(JointTypeRevoluteX),
                body12);

    unsigned int reference_body_id22 =
       reference2.AddBody (reference_body_id21,
                SpatialTransformd(rm1,r1),
                Joint(JointTypeEulerZYX),
                body13);

    unsigned int custom_body_id20 =
       custom2.AddBody (  0,
                SpatialTransformd(),
                Joint(JointTypeRevoluteX),
                body11);

    unsigned int custom_body_id21 =
       custom2.AddBodyCustomJoint (  custom_body_id20,
                SpatialTransformd(rm2,r2),
                &custom_rx_joint1,
                body12);

    unsigned int custom_body_id22 =
       custom2.AddBody (  custom_body_id21,
                SpatialTransformd(rm1,r1),
                Joint(JointTypeEulerZYX),
                body13);



    VectorNd q2         = VectorNd::Zero (reference2.q_size);
    VectorNd qdot2      = VectorNd::Zero (reference2.qdot_size);
    VectorNd qddot2     = VectorNd::Zero (reference2.qdot_size);
    VectorNd tau2       = VectorNd::Zero (reference2.qdot_size);

    idx = 1;
    reference_model.at(idx) = (reference2);
    custom_model.at(idx) = (custom2);

    reference_body_id.at(idx).at(0) = (reference_body_id20);
    reference_body_id.at(idx).at(1) = (reference_body_id21);
    reference_body_id.at(idx).at(2) = (reference_body_id22);

    custom_body_id.at(idx).at(0) = (custom_body_id20);
    custom_body_id.at(idx).at(1) = (custom_body_id21);
    custom_body_id.at(idx).at(2) = (custom_body_id22);

    body.at(idx).at(0) = (body11);
    body.at(idx).at(1) = (body12);
    body.at(idx).at(2) = (body13);
    custom_joint.at(idx).at(0) = (&custom_rx_joint1);


    q.at(idx)     = (q2);
    qdot.at(idx)  = (qdot2);
    qddot.at(idx) = (qddot2);
    tau.at(idx)   = (tau2);

    //========================================================
    //Test Model 3: custom - multidof3 - Rx
    //========================================================

    Model reference3, custom3;

    unsigned int reference_body_id30 =
       reference3.AddBody (0,
                SpatialTransformd(),
                Joint(JointTypeRevoluteX),
                body11);

    unsigned int reference_body_id31 =
       reference3.AddBody (reference_body_id30,
                SpatialTransformd(rm1,r1),
                Joint(JointTypeEulerZYX),
                body12);

    unsigned int reference_body_id32 =
       reference3.AddBody (reference_body_id31,
                SpatialTransformd(rm2,r2),
                Joint(JointTypeRevoluteX),
                body13);

    unsigned int custom_body_id30 =
       custom3.AddBodyCustomJoint (  0,
                SpatialTransformd(),
                &custom_rx_joint1,
                body11);

    unsigned int custom_body_id31 =
       custom3.AddBody (  custom_body_id30,
                SpatialTransformd(rm1,r1),
                Joint(JointTypeEulerZYX),
                body12);

    unsigned int custom_body_id32 =
       custom3.AddBody (  custom_body_id31,
                SpatialTransformd(rm2,r2),
                Joint(JointTypeRevoluteX),
                body13);

    VectorNd q3     = VectorNd::Zero (reference3.q_size);
    VectorNd qdot3  = VectorNd::Zero (reference3.qdot_size);
    VectorNd qddot3   = VectorNd::Zero (reference3.qdot_size);
    VectorNd tau3   = VectorNd::Zero (reference3.qdot_size);

    idx = 2;
    reference_model.at(idx) = (reference3);
    custom_model.at(idx) = (custom3);

    reference_body_id.at(idx).at(0) = (reference_body_id30);
    reference_body_id.at(idx).at(1) = (reference_body_id31);
    reference_body_id.at(idx).at(2) = (reference_body_id32);

    custom_body_id.at(idx).at(0) = (custom_body_id30);
    custom_body_id.at(idx).at(1) = (custom_body_id31);
    custom_body_id.at(idx).at(2) = (custom_body_id32);

    body.at(idx).at(0) = (body11);
    body.at(idx).at(1) = (body12);
    body.at(idx).at(2) = (body13);
    custom_joint.at(idx).at(0) = (&custom_rx_joint1);

    q.at(idx)       = (q3);
    qdot.at(idx)    = (qdot3);
    qddot.at(idx)   = (qddot3);
    tau.at(idx)     = (tau3);

  }


  virtual void TearDown () {
  }

  vector < Model > reference_model;
  vector < Model > custom_model;

  vector < vector < Body >  > body;
  vector < vector < CustomJoint* > > custom_joint;

  vector < vector< unsigned int > > reference_body_id;
  vector < vector< unsigned int > > custom_body_id;

  vector < VectorNd > q;
  vector < VectorNd > qdot;
  vector < VectorNd > qddot;
  vector < VectorNd > tau;

  CustomJointTypeRevoluteX custom_rx_joint1;
  CustomJointTypeRevoluteX custom_rx_joint2;
  CustomJointTypeRevoluteX custom_rx_joint3;
};


//==============================================================================
//
// Tests
//  UpdateKinematicsCustom
//  Jacobians
//  InverseDynamics
//  CompositeRigidBodyAlgorithm
//  ForwardDynamics
//  CalcMInvTimestau
//  ForwardDynamicsContactsKokkevis
//
//==============================================================================

TEST_F ( CustomJointMultiBodyFixture, UpdateKinematics ) {

  VectorNd test;

  for(int idx =0; idx < NUMBER_OF_MODELS; ++idx){

    int dof = reference_model.at(idx).dof_count;
    for (unsigned int i = 0; i < dof ; i++) {
      q.at(idx)[i]      = i * 0.1;
      qdot.at(idx)[i]   = i * 0.15;
      qddot.at(idx)[i]  = i * 0.17;
    }

    UpdateKinematics (reference_model.at(idx),
                      q.at(idx),
                      qdot.at(idx),
                      qddot.at(idx));

    UpdateKinematics (custom_model.at(idx),
                      q.at(idx),
                      qdot.at(idx),
                      qddot.at(idx));


    Matrix3d Eref = reference_model.at(idx).X_base[
                      reference_body_id.at(idx).at(NUMBER_OF_BODIES-1)
                      ].E;
    Matrix3d Ecus = custom_model.at(idx).X_base[
                      custom_body_id.at(idx).at(NUMBER_OF_BODIES-1)
                      ].E;

    Matrix3d Eerr = Eref-Ecus;

    EXPECT_TRUE(EIGEN_MATRIX_NEAR (
      reference_model.at(idx).X_base[
        reference_body_id.at(idx).at(NUMBER_OF_BODIES-1)
        ].E,
      custom_model.at(idx).X_base[
        custom_body_id.at(idx).at(NUMBER_OF_BODIES-1)
        ].E,

      TEST_PREC));

    EXPECT_TRUE(EIGEN_MATRIX_NEAR (
      reference_model.at(idx).model_data.v[
        reference_body_id.at(idx).at(NUMBER_OF_BODIES-1)
        ],
      custom_model.at(idx).model_data.v[
        custom_body_id.at(idx).at(NUMBER_OF_BODIES-1)
        ],
      TEST_PREC));

    EXPECT_TRUE(EIGEN_MATRIX_NEAR (
      reference_model.at(idx).model_data.a[
        reference_body_id.at(idx).at(NUMBER_OF_BODIES-1)
        ],
      custom_model.at(idx).model_data.a[
        custom_body_id.at(idx).at(NUMBER_OF_BODIES-1)
        ],
      TEST_PREC));
  }
}

TEST_F (CustomJointMultiBodyFixture, UpdateKinematicsCustom) {

  for(int idx =0; idx < NUMBER_OF_MODELS; ++idx){
    int dof = reference_model.at(idx).dof_count;
    for (unsigned int i = 0; i < dof; i++) {
      q.at(idx)[i]        = i * 9.133758561390194e-01;
      qdot.at(idx)[i]     = i * 6.323592462254095e-01;
      qddot.at(idx)[i]    = i * 9.754040499940952e-02;
    }

    UpdateKinematicsCustom (reference_model.at(idx),
                            &q.at(idx), NULL, NULL);
    UpdateKinematicsCustom (custom_model.at(idx),
                            &q.at(idx), NULL, NULL);


    EXPECT_TRUE(EIGEN_MATRIX_NEAR (
      reference_model.at(idx).X_base[
        reference_body_id.at(idx).at(NUMBER_OF_BODIES-1)
        ].E,
      custom_model.at(idx).X_base[
        custom_body_id.at(idx).at(NUMBER_OF_BODIES-1)
        ].E,
      TEST_PREC));


    //velocity
    UpdateKinematicsCustom (reference_model.at(idx),
                            &q.at(idx),
                            &qdot.at(idx), NULL);
    UpdateKinematicsCustom (custom_model.at(idx),
                            &q.at(idx),
                            &qdot.at(idx), NULL);

    EXPECT_TRUE(EIGEN_MATRIX_NEAR (
          reference_model.at(idx).model_data.v[
            reference_body_id.at(idx).at(NUMBER_OF_BODIES-1)
            ],
          custom_model.at(idx).model_data.v[
            custom_body_id.at(idx).at(NUMBER_OF_BODIES-1)
            ],
          TEST_PREC));


    //All
    UpdateKinematicsCustom (reference_model.at(idx),
                            &q.at(idx),
                            &qdot.at(idx),
                            &qddot.at(idx));

    UpdateKinematicsCustom (custom_model.at(idx),
                            &q.at(idx),
                            &qdot.at(idx),
                            &qddot.at(idx));

    EXPECT_TRUE(EIGEN_MATRIX_NEAR (
      reference_model.at(idx).model_data.a[
        reference_body_id.at(idx).at(NUMBER_OF_BODIES-1)
        ],
      custom_model.at(idx).model_data.a[
        custom_body_id.at(idx).at(NUMBER_OF_BODIES-1)
        ],
      TEST_PREC));
  }


}

TEST_F (CustomJointMultiBodyFixture, Jacobians) {

  for(int idx = 0; idx < NUMBER_OF_MODELS; ++idx){
    int dof = reference_model.at(idx).dof_count;

    for (unsigned int i = 0; i < dof; i++) {
      q.at(idx)[i]        = i * 9.133758561390194e-01;
      qdot.at(idx)[i]     = i * 6.323592462254095e-01;
      qddot.at(idx)[i]    = i * 9.754040499940952e-02;
    }

    //position
    UpdateKinematics (reference_model.at(idx),
                      q.at(idx),
                      qdot.at(idx),
                      qddot.at(idx));

    UpdateKinematics (custom_model.at(idx),
                      q.at(idx),
                      qdot.at(idx),
                      qddot.at(idx));

    //Check the Spatial Jacobian
    MatrixNd Gref = MatrixNd(
          MatrixNd::Zero(6,reference_model.at(idx).dof_count));

    MatrixNd Gcus = MatrixNd(
          MatrixNd::Zero(6,reference_model.at(idx).dof_count));

    CalcBodySpatialJacobian ( reference_model.at(idx),
                              q.at(idx),
                              reference_body_id.at(idx).at(NUMBER_OF_BODIES-1),
                              Gref);

    CalcBodySpatialJacobian ( custom_model.at(idx),
                              q.at(idx),
                              custom_body_id.at(idx).at(NUMBER_OF_BODIES-1),
                              Gcus);

    for(int i=0; i<6;++i){
      for(int j=0; j<dof;++j){
        EXPECT_NEAR (
          Gref(i,j),
          Gcus(i,j),
          TEST_PREC);
      }
    }

    //Check the 6d point Jacobian
    Vector3d point_position (1.1, 1.2, 2.1);


    CalcPointJacobian6D (reference_model.at(idx),
                         q.at(idx),
                         reference_body_id.at(idx).at(NUMBER_OF_BODIES-1),
                         point_position,
                         Gref);

    CalcPointJacobian6D (custom_model.at(idx),
                         q.at(idx),
                         custom_body_id.at(idx).at(NUMBER_OF_BODIES-1),
                         point_position,
                         Gcus);

    for(int i=0; i<6;++i){
      for(int j=0; j<dof;++j){
        EXPECT_NEAR (
          Gref(i,j),
          Gcus(i,j),
          TEST_PREC);
      }
    }

    //Check the 3d point Jacobian
    MatrixNd GrefPt = MatrixNd::Constant (3,
                      reference_model.at(idx).dof_count,
                      0.);
    MatrixNd GcusPt = MatrixNd::Constant (3,
                      reference_model.at(idx).dof_count,
                      0.);

    CalcPointJacobian (reference_model.at(idx),
                       q.at(idx),
                       reference_body_id.at(idx).at(NUMBER_OF_BODIES-1),
                       point_position,
                       GrefPt);

    CalcPointJacobian (custom_model.at(idx),
                       q.at(idx),
                       custom_body_id.at(idx).at(NUMBER_OF_BODIES-1),
                       point_position,
                       GcusPt);

    for(int i=0; i<3;++i){
      for(int j=0; j<dof;++j){
        EXPECT_NEAR (
          GrefPt(i,j),
          GcusPt(i,j),
          TEST_PREC);
      }
    }
  }

}

TEST_F (CustomJointMultiBodyFixture, InverseDynamics) {

  for(int idx =0; idx < NUMBER_OF_MODELS; ++idx){

    int dof = reference_model.at(idx).dof_count;

    for (unsigned int i = 0; i < dof; i++) {
      q.at(idx)[i]        = i * 9.133758561390194e-01;
      qdot.at(idx)[i]     = i * 6.323592462254095e-01;
      qddot.at(idx)[i]    = i * 9.754040499940952e-02;
    }

    //position
    VectorNd tauRef = VectorNd::Zero (reference_model.at(idx).qdot_size);
    VectorNd tauCus = VectorNd::Zero (reference_model.at(idx).qdot_size);

    InverseDynamics(reference_model.at(idx),
                    q.at(idx),
                    qdot.at(idx),
                    qddot.at(idx),
                    tauRef);

    InverseDynamics(custom_model.at(idx),
                    q.at(idx),
                    qdot.at(idx),
                    qddot.at(idx),
                    tauCus);

    VectorNd tauErr = tauRef-tauCus;

    EXPECT_TRUE(EIGEN_MATRIX_NEAR (
      tauRef,
      tauCus,
      TEST_PREC));
  }

}

TEST_F (CustomJointMultiBodyFixture, CompositeRigidBodyAlgorithm) {

  for(int idx =0; idx < NUMBER_OF_MODELS; ++idx){

    int dof = reference_model.at(idx).dof_count;

    for (unsigned int i = 0; i < dof; i++) {
      q.at(idx)[i]    = (i+0.1) * 9.133758561390194e-01;
      qdot.at(idx)[i] = (i+0.1) * 6.323592462254095e-01;
      tau.at(idx)[i]  = (i+0.1) * 9.754040499940952e-02;
    }

    MatrixNd h_ref = MatrixNd::Constant (dof, dof, 0.);
    VectorNd c_ref = VectorNd::Constant (dof, 0.);
    VectorNd qddot_zero_ref = VectorNd::Constant (dof, 0.);
    VectorNd qddot_crba_ref = VectorNd::Constant (dof, 0.);

    MatrixNd h_cus = MatrixNd::Constant (dof, dof, 0.);
    VectorNd c_cus = VectorNd::Constant (dof, 0.);
    VectorNd qddot_zero_cus = VectorNd::Constant (dof, 0.);
    VectorNd qddot_crba_cus = VectorNd::Constant (dof, 0.);

    VectorNd qddotRef = VectorNd::Zero (dof);
    VectorNd qddotCus = VectorNd::Zero (dof);

    //Ref
    ForwardDynamics(reference_model.at(idx),
                    q.at(idx),
                    qdot.at(idx),
                    tau.at(idx),
                    qddotRef);

    CompositeRigidBodyAlgorithm (reference_model.at(idx),
                                 q.at(idx),
                                 h_ref);

    InverseDynamics (reference_model.at(idx),
                     q.at(idx),
                     qdot.at(idx),
                     qddot_zero_ref,
                     c_ref);

    LinSolveGaussElimPivot (h_ref,
                            c_ref * -1. + tau.at(idx),
                            qddot_crba_ref);

    //Custom
    ForwardDynamics(custom_model.at(idx),
                    q.at(idx),
                    qdot.at(idx),
                    tau.at(idx),
                    qddotCus);

    CompositeRigidBodyAlgorithm (custom_model.at(idx),
                                 q.at(idx),
                                 h_cus);

    InverseDynamics (custom_model.at(idx),
                     q.at(idx),
                     qdot.at(idx),
                     qddot_zero_cus,
                     c_cus);

    LinSolveGaussElimPivot (h_cus,
                            c_cus * -1. + tau.at(idx),
                            qddot_crba_cus);

    EXPECT_TRUE(EIGEN_MATRIX_NEAR(qddot_crba_ref,
                      qddot_crba_cus,
                      TEST_PREC));
  }
}

TEST_F (CustomJointMultiBodyFixture, ForwardDynamics) {

  for(int idx =0; idx < NUMBER_OF_MODELS; ++idx){

    int dof = reference_model.at(idx).dof_count;

    for (unsigned int i = 0; i < dof; i++) {
      q.at(idx)[i]       = (i+0.1) * 9.133758561390194e-01;
      qdot.at(idx)[i]    = (i+0.1) * 6.323592462254095e-01;
      qddot.at(idx)[i]   = (i+0.1) * 2.323592499940952e-01;
      tau.at(idx)[i]     = (i+0.1) * 9.754040499940952e-02;
    }


    VectorNd qddotRef = VectorNd::Zero (reference_model.at(idx).qdot_size);
    VectorNd qddotCus = VectorNd::Zero (reference_model.at(idx).qdot_size);

    ForwardDynamics(reference_model.at(idx),
                    q.at(idx),
                    qdot.at(idx),
                    tau.at(idx),
                    qddotRef);

    ForwardDynamics(custom_model.at(idx),
                    q.at(idx),
                    qdot.at(idx),
                    tau.at(idx),
                    qddotCus);

    EXPECT_TRUE(EIGEN_MATRIX_NEAR (
      qddotRef,
      qddotCus,
      TEST_PREC));
  }

}

TEST_F (CustomJointMultiBodyFixture, CalcMInvTimestau) {

  for(int idx =0; idx < NUMBER_OF_MODELS; ++idx){

    int dof = reference_model.at(idx).dof_count;

    for (unsigned int i = 0; i < dof; i++) {
      q.at(idx)[i]    = (i+0.1) * 9.133758561390194e-01;
      tau.at(idx)[i]  = (i+0.1) * 9.754040499940952e-02;

    }

    //reference
    VectorNd qddot_minv_ref = VectorNd::Zero(dof);

    CalcMInvTimesTau(reference_model.at(idx),
                     q.at(idx),
                     tau.at(idx),
                     qddot_minv_ref,
                     true);

    //custom
    VectorNd qddot_minv_cus = VectorNd::Zero(dof);
    CalcMInvTimesTau(custom_model.at(idx),
                     q.at(idx),
                     tau.at(idx),
                     qddot_minv_cus,
                     true);

    //check.
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(qddot_minv_ref,
                      qddot_minv_cus,
                      TEST_PREC));
  }

}

TEST_F (CustomJointMultiBodyFixture, ForwardDynamicsContactsKokkevis){

  for(int idx =0; idx < NUMBER_OF_MODELS; ++idx){

    int dof = reference_model.at(idx).dof_count;

    //Adding a 1 constraint to a system with 1 dof is
    //a no-no
    if(dof > 1){

      for (unsigned int i = 0; i < dof; i++) {
        q.at(idx)[i]      = (i+0.1) * 9.133758561390194e-01;
        qdot.at(idx)[i]   = (i+0.1) * 6.323592462254095e-01;
        tau.at(idx)[i]    = (i+0.1) * 9.754040499940952e-02;
      }

      VectorNd qddot_ref = VectorNd::Zero(dof);
      VectorNd qddot_cus = VectorNd::Zero(dof);

      VectorNd qdot_plus_ref = VectorNd::Zero(dof);
      VectorNd qdot_plus_cus = VectorNd::Zero(dof);

      Vector3d contact_point ( 0., 1., 0.);

      ConstraintSet constraint_set_ref;
      ConstraintSet constraint_set_cus;


      //Reference
      constraint_set_ref.AddConstraint(
                          reference_body_id.at(idx).at(NUMBER_OF_BODIES-1),
                          contact_point,
                          Vector3d (1., 0., 0.),
                          "ground_x");

      constraint_set_ref.AddConstraint(
                          reference_body_id.at(idx).at(NUMBER_OF_BODIES-1),
                          contact_point,
                          Vector3d (0., 1., 0.),
                          "ground_y");

      constraint_set_ref.Bind (reference_model.at(idx));

      //Custom
      constraint_set_cus.AddConstraint(
                          custom_body_id.at(idx).at(NUMBER_OF_BODIES-1),
                          contact_point,
                          Vector3d (1., 0., 0.),
                          "ground_x");

      constraint_set_cus.AddConstraint(
                          custom_body_id.at(idx).at(NUMBER_OF_BODIES-1),
                          contact_point,
                          Vector3d (0., 1., 0.),
                          "ground_y");

      constraint_set_cus.Bind (custom_model.at(idx));

      ComputeContactImpulsesDirect(reference_model.at(idx),
                                   q.at(idx),
                                   qdot.at(idx),
                                   constraint_set_ref,
                                   qdot_plus_ref);

      ForwardDynamicsContactsKokkevis (reference_model.at(idx),
                                       q.at(idx),
                                       qdot_plus_ref,
                                       tau.at(idx),
                                       constraint_set_ref,
                                       qddot_ref);

      ComputeContactImpulsesDirect(custom_model.at(idx),
                                   q.at(idx),
                                   qdot.at(idx),
                                   constraint_set_cus,
                                   qdot_plus_cus);

      ForwardDynamicsContactsKokkevis (custom_model.at(idx),
                                       q.at(idx),
                                       qdot_plus_cus,
                                       tau.at(idx),
                                       constraint_set_cus,
                                       qddot_cus);

      VectorNd qdot_plus_error = qdot_plus_ref - qdot_plus_cus;
      VectorNd qddot_error = qddot_ref - qddot_cus;

      EXPECT_TRUE(EIGEN_MATRIX_NEAR (qdot_plus_ref,
                         qdot_plus_cus,
                         TEST_PREC));

      EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_ref,
                         qddot_cus,
                         TEST_PREC));
    }
  }

}

//
//Completed?
// x  : implement test for UpdateKinematicsCustom
// x  : implement test for Jacobians
// x  : implement test for InverseDynamics
// x  : implement test for CompositeRigidBodyAlgorithm
// x  : implement test for ForwardDynamics
// x  : implement test for CalcMInvTimestau
// x  : implement test for ForwardDynamicsContactsKokkevis

int main( int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
