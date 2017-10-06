/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <iostream>
#include <limits>
#include <assert.h>
#include <string.h>

#include "rbdl/rbdl_mathutils.h"
#include "rbdl/Logging.h"

#include "rbdl/Model.h"
#include "rbdl/Joint.h"
#include "rbdl/Body.h"
#include "rbdl/Dynamics.h"
#include "rbdl/Kinematics.h"

namespace RigidBodyDynamics {

using namespace Math;

RBDL_DLLAPI void NonlinearEffects ( 
    Model &model,
    ModelDatad &model_data,
    const VectorNd &Q,
    const VectorNd &QDot,
    VectorNd &Tau) {
  LOG << "-------- " << __func__ << " --------" << std::endl;

  SpatialVectord spatial_gravity (0., 0., 0., -model.gravity[0], -model.gravity[1], -model.gravity[2]);

  // Reset the velocity of the root body
  model_data.v[0].setZero();
  model_data.a[0] = spatial_gravity;

  for (unsigned int i = 1; i < model.mJointUpdateOrder.size(); i++) {
    jcalc (model, model_data, model.mJointUpdateOrder[i], Q, QDot);
  }

  for (unsigned int i = 1; i < model.mBodies.size(); i++) {
    if (model.lambda[i] == 0) {
      model_data.v[i] = model_data.v_J[i];
      model_data.a[i] = model_data.X_lambda[i].apply(spatial_gravity);
    }	else {
      model_data.v[i] = model_data.X_lambda[i].apply(model_data.v[model.lambda[i]]) + model_data.v_J[i];
      model_data.c[i] = model_data.c_J[i] + crossm(model_data.v[i],model_data.v_J[i]);
      model_data.a[i] = model_data.X_lambda[i].apply(model_data.a[model.lambda[i]]) + model_data.c[i];
    }

    if (!model.mBodies[i].mIsVirtual) {
      model_data.f[i] = model_data.I[i] * model_data.a[i] + crossf(model_data.v[i],model_data.I[i] * model_data.v[i]);
    } else {
      model_data.f[i].setZero();
    }
  }

  for (unsigned int i = model.mBodies.size() - 1; i > 0; i--) {
    //if(model.mJoints[i].mJointType != JointTypeCustom){
      if (model.mJoints[i].mDoFCount == 1) {
        Tau[model.mJoints[i].q_index]
          = model_data.S[i].dot(model_data.f[i]);
      } else if (model.mJoints[i].mDoFCount == 3) {
        Tau.block<3,1>(model.mJoints[i].q_index, 0)
          = model_data.multdof3_S[i].transpose() * model_data.f[i];
      }
//    } else if(model.mJoints[i].mJointType == JointTypeCustom) {
//      unsigned int k = model.mJoints[i].custom_joint_index;
//      Tau.block(model.mJoints[i].q_index,0,
//          model.mCustomJoints[k]->mDoFCount, 1)
//        = model.mCustomJoints[k]->S.transpose() * model_data.f[i];
//    }

    if (model.lambda[i] != 0) {
      model_data.f[model.lambda[i]] = model_data.f[model.lambda[i]] + model_data.X_lambda[i].applyTranspose(model_data.f[i]);
    }
  }
}

RBDL_DLLAPI void ForwardDynamicsLagrangian (
    Model &model,
    ModelDatad &model_data,
    const VectorNd &Q,
    const VectorNd &QDot,
    const VectorNd &Tau,
    VectorNd &QDDot,
    Math::LinearSolver linear_solver,
    std::vector<SpatialVectord> *f_ext,
    Math::MatrixNd *H,
    Math::VectorNd *C) {
  LOG << "-------- " << __func__ << " --------" << std::endl;

  bool free_H = false;
  bool free_C = false;

  if (H == NULL) {
    H = new MatrixNd (MatrixNd::Zero(model.dof_count, model.dof_count));
    free_H = true;
  }

  if (C == NULL) {
    C = new VectorNd (VectorNd::Zero(model.dof_count));
    free_C = true;
  }

  // we set QDDot to zero to compute C properly with the InverseDynamics
  // method.
  QDDot.setZero();

  InverseDynamics (model, model_data, Q, QDot, QDDot, (*C), f_ext);
  CompositeRigidBodyAlgorithm (model, model_data, Q, *H, false);

  LOG << "A = " << std::endl << *H << std::endl;
  LOG << "b = " << std::endl << *C * -1. + Tau << std::endl;

#ifndef RBDL_USE_SIMPLE_MATH
  switch (linear_solver) {
    case (LinearSolverPartialPivLU) :
      QDDot = H->partialPivLu().solve (*C * -1. + Tau);
      break;
    case (LinearSolverColPivHouseholderQR) :
      QDDot = H->colPivHouseholderQr().solve (*C * -1. + Tau);
      break;
    case (LinearSolverHouseholderQR) :
      QDDot = H->householderQr().solve (*C * -1. + Tau);
      break;
    case (LinearSolverLLT) :
      QDDot = H->llt().solve (*C * -1. + Tau);
      break;
    default:
      LOG << "Error: Invalid linear solver: " << linear_solver << std::endl;
      assert (0);
      break;
  }
#else
  bool solve_successful = LinSolveGaussElimPivot (*H, *C * -1. + Tau, QDDot);
  assert (solve_successful);
#endif

  if (free_C) {
    delete C;
  }

  if (free_H) {
    delete H;
  }

  LOG << "x = " << QDDot << std::endl;
}

RBDL_DLLAPI void CalcMInvTimesTau ( Model &model,
    ModelDatad &model_data,
    const VectorNd &Q,
    const VectorNd &Tau,
    VectorNd &QDDot,
    bool update_kinematics) {

  LOG << "Q          = " << Q.transpose() << std::endl;
  LOG << "---" << std::endl;

  // Reset the velocity of the root body
  model_data.v[0].setZero();
  model_data.a[0].setZero();

  if (update_kinematics) {
    for (unsigned int i = 1; i < model.mBodies.size(); i++) {
      jcalc_X_lambda_S (model, model_data, model.mJointUpdateOrder[i], Q);

      model_data.v_J[i].setZero();
      model_data.v[i].setZero();
      model_data.c_J[i].setZero();
      model_data.pA[i].setZero();
      model_data.I[i].setSpatialMatrix(model_data.IA[i]);
    }
  }

  for (unsigned int i = 1; i < model.mBodies.size(); i++) {
    model_data.pA[i].setZero();
  }

  // ClearLogOutput();

  if (update_kinematics) {
    // Compute Articulate Body Inertias
    for (unsigned int i = model.mBodies.size() - 1; i > 0; i--) {
      unsigned int q_index = model.mJoints[i].q_index;

      if (model.mJoints[i].mDoFCount == 1){
          //&& model.mJoints[i].mJointType != JointTypeCustom) {
        model_data.U[i] = model_data.IA[i] * model_data.S[i];
        model_data.d[i] = model_data.S[i].dot(model_data.U[i]);
        //      LOG << "u[" << i << "] = " << model.u[i] << std::endl;
        unsigned int lambda = model.lambda[i];

        if (lambda != 0) {
          SpatialMatrixd Ia = model_data.IA[i] -
            model_data.U[i] * (model_data.U[i] / model_data.d[i]).transpose();
#ifdef EIGEN_CORE_H
          model_data.IA[lambda].noalias() += model_data.X_lambda[i].toMatrixTranspose()
            * Ia
            * model_data.X_lambda[i].toMatrix();
#else
          model.IA[lambda] += model_data.X_lambda[i].toMatrixTranspose()
            * Ia
            * model_data.X_lambda[i].toMatrix();
#endif
        }
      } else if (model.mJoints[i].mDoFCount == 3){
        //  && model.mJoints[i].mJointType != JointTypeCustom) {

        model_data.multdof3_U[i] = model_data.IA[i] * model_data.multdof3_S[i];

#ifdef EIGEN_CORE_H
        model_data.multdof3_Dinv[i] =
          (model_data.multdof3_S[i].transpose()*model_data.multdof3_U[i]).inverse().eval();
#else
        model.multdof3_Dinv[i] = 
          (model_data.multdof3_S[i].transpose() * model_data.multdof3_U[i]).inverse();
#endif
        //      LOG << "mCustomJoints[kI]->u[" << i << "] = "
        //<< model.mCustomJoints[kI]->u[i].transpose() << std::endl;

        unsigned int lambda = model.lambda[i];

        if (lambda != 0) {
          SpatialMatrixd Ia = model_data.IA[i]
            - ( model_data.multdof3_U[i]
                * model_data.multdof3_Dinv[i]
                * model_data.multdof3_U[i].transpose());
#ifdef EIGEN_CORE_H
          model_data.IA[lambda].noalias() +=
            model_data.X_lambda[i].toMatrixTranspose()
            * Ia
            * model_data.X_lambda[i].toMatrix();
#else
          model.IA[lambda] +=
            model_data.X_lambda[i].toMatrixTranspose()
            * Ia * model_data.X_lambda[i].toMatrix();
#endif
        }
      }/* else if (model.mJoints[i].mJointType == JointTypeCustom) {
        unsigned int kI     = model.mJoints[i].custom_joint_index;
        unsigned int dofI   = model.mCustomJoints[kI]->mDoFCount;
        model.mCustomJoints[kI]->U = model.IA[i] * model.mCustomJoints[kI]->S;

#ifdef EIGEN_CORE_H
        model.mCustomJoints[kI]->Dinv = (model.mCustomJoints[kI]->S.transpose()
            * model.mCustomJoints[kI]->U
            ).inverse().eval();
#else
        model.mCustomJoints[kI]->Dinv[i]=(model.mCustomJoints[kI]->S.transpose()
            * model.mCustomJoints[kI]->U
            ).inverse();
#endif
        //      LOG << "mCustomJoints[kI]->u[" << i << "] = "
        //<< model.mCustomJoints[kI]->u.transpose() << std::endl;
        unsigned int lambda = model.lambda[i];

        if (lambda != 0) {
          SpatialMatrixd Ia = model.IA[i]
            - ( model.mCustomJoints[kI]->U
                * model.mCustomJoints[kI]->Dinv
                * model.mCustomJoints[kI]->U.transpose());
#ifdef EIGEN_CORE_H
          model.IA[lambda].noalias() += model_data.X_lambda[i].toMatrixTranspose()
            * Ia
            * model_data.X_lambda[i].toMatrix();
#else
          model.IA[lambda] += model_data.X_lambda[i].toMatrixTranspose()
            * Ia * model_data.X_lambda[i].toMatrix();
#endif
        }

      }
      */
    }
  }

  // compute articulated bias forces
  for (unsigned int i = model.mBodies.size() - 1; i > 0; i--) {
    unsigned int q_index = model.mJoints[i].q_index;

    if (model.mJoints[i].mDoFCount == 1){
        //&& model.mJoints[i].mJointType != JointTypeCustom) {

      model_data.u[i] = Tau[q_index] - model_data.S[i].dot(model_data.pA[i]);
      // LOG << "u[" << i << "] = " << model.u[i] << std::endl;
      unsigned int lambda = model.lambda[i];
      if (lambda != 0) {
        SpatialVectord pa = model_data.pA[i] + model_data.U[i] * model_data.u[i] / model_data.d[i];

#ifdef EIGEN_CORE_H
        model_data.pA[lambda].noalias() += model_data.X_lambda[i].applyTranspose(pa);
#else
        model.pA[lambda] += model_data.X_lambda[i].applyTranspose(pa);
#endif
        LOG << "pA[" << lambda << "] = "
          << model_data.pA[lambda].transpose() << std::endl;
      }
    } else if (model.mJoints[i].mDoFCount == 3){
        //&& model.mJoints[i].mJointType != JointTypeCustom) {

      Vector3d tau_temp ( Tau[q_index],
          Tau[q_index + 1],
          Tau[q_index + 2]);
      model_data.multdof3_u[i] = tau_temp
        - model_data.multdof3_S[i].transpose()*model_data.pA[i];
      //      LOG << "multdof3_u[" << i << "] = "
      // << model.multdof3_u[i].transpose() << std::endl;
      unsigned int lambda = model.lambda[i];

      if (lambda != 0) {
        SpatialVectord pa = model_data.pA[i]
          + model_data.multdof3_U[i]
          * model_data.multdof3_Dinv[i]
          * model_data.multdof3_u[i];

#ifdef EIGEN_CORE_H
        model_data.pA[lambda].noalias() +=
          model_data.X_lambda[i].applyTranspose(pa);
#else
        model.pA[lambda] += model_data.X_lambda[i].applyTranspose(pa);
#endif
        LOG << "pA[" << lambda << "] = "
          << model_data.pA[lambda].transpose() << std::endl;
      }
    } /*else if (model.mJoints[i].mJointType == JointTypeCustom) {
      unsigned int kI     = model.mJoints[i].custom_joint_index;
      unsigned int dofI   = model.mCustomJoints[kI]->mDoFCount;
      VectorNd tau_temp   = VectorNd::Zero(dofI);

      for(int z=0; z<dofI;++z){
        tau_temp(z) = Tau[q_index+z];
      }
      model.mCustomJoints[kI]->u = 
        tau_temp - ( model.mCustomJoints[kI]->S.transpose()* model.pA[i]);
      //      LOG << "mCustomJoints[kI]->u"
      // << model.mCustomJoints[kI]->u.transpose() << std::endl;
      unsigned int lambda = model.lambda[i];

      if (lambda != 0) {
        SpatialVectord pa = model.pA[i]
          + (   model.mCustomJoints[kI]->U
              * model.mCustomJoints[kI]->Dinv
              * model.mCustomJoints[kI]->u);

#ifdef EIGEN_CORE_H
        model.pA[lambda].noalias() +=
          model_data.X_lambda[i].applyTranspose(pa);
#else
        model.pA[lambda] += model_data.X_lambda[i].applyTranspose(pa);
#endif
        LOG << "pA[" << lambda << "] = "
          << model.pA[lambda].transpose() << std::endl;
      }
    }*/
  }

  //  ClearLogOutput();

  for (unsigned int i = 1; i < model.mBodies.size(); i++) {
    unsigned int q_index = model.mJoints[i].q_index;
    unsigned int lambda = model.lambda[i];
    SpatialTransformd X_lambda = model_data.X_lambda[i];

    model_data.a[i] = X_lambda.apply(model_data.a[lambda]) + model_data.c[i];
    LOG << "a'[" << i << "] = " << model_data.a[i].transpose() << std::endl;

    if (model.mJoints[i].mDoFCount == 1){
     //   && model.mJoints[i].mJointType != JointTypeCustom) {
      QDDot[q_index] = (1./model_data.d[i])*(model_data.u[i]-model_data.U[i].dot(model_data.a[i]));
      model_data.a[i]     = model_data.a[i] + model_data.S[i] * QDDot[q_index];
    } else if (model.mJoints[i].mDoFCount == 3){
     //   && model.mJoints[i].mJointType != JointTypeCustom) {
      Vector3d qdd_temp = 
        model_data.multdof3_Dinv[i] * (model_data.multdof3_u[i]
            - model_data.multdof3_U[i].transpose()*model_data.a[i]);

      QDDot[q_index]      = qdd_temp[0];
      QDDot[q_index + 1]  = qdd_temp[1];
      QDDot[q_index + 2]  = qdd_temp[2];
      model_data.a[i]          = model_data.a[i] + model_data.multdof3_S[i] * qdd_temp;
    }/* else if (model.mJoints[i].mJointType == JointTypeCustom) {
      unsigned int kI     = model.mJoints[i].custom_joint_index;
      unsigned int dofI   = model.mCustomJoints[kI]->mDoFCount;

      VectorNd qdd_temp = model.mCustomJoints[kI]->Dinv
        * (  model.mCustomJoints[kI]->u 
            - model.mCustomJoints[kI]->U.transpose() * model_data.a[i]);

      for(int z=0; z<dofI;++z){
        QDDot[q_index+z]      = qdd_temp[z];
      }

      model_data.a[i] =    model_data.a[i]
        + model.mCustomJoints[kI]->S * qdd_temp;
    }*/
  }

  LOG << "QDDot = " << QDDot.transpose() << std::endl;
}

} /* namespace RigidBodyDynamics */
