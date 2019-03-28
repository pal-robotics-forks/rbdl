/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef RBDL_DYNAMICS_H
#define RBDL_DYNAMICS_H

#include <assert.h>
#include <iostream>

#include "rbdl/rbdl_math.h"
#include "rbdl/rbdl_mathutils.h"

#include "rbdl/Logging.h"
#include <rbdl/ModelData.h>
#include "rbdl/Joint.h"

namespace RigidBodyDynamics
{
struct Model;

/** \page dynamics_page Dynamics
 *
 * All functions related to kinematics are specified in the \ref
 * dynamics_group "Dynamics Module".
 *
 * \defgroup dynamics_group Dynamics
 * @{
 */

/** \brief Computes inverse dynamics with the Newton-Euler Algorithm
 *
 * This function computes the generalized forces from given generalized
 * states, velocities, and accelerations:
 *   \f$ \tau = M(q) \ddot{q} + N(q, \dot{q}) \f$
 *
 * \param model rigid body model
 * \param Q     state vector of the internal joints
 * \param QDot  velocity vector of the internal joints
 * \param QDDot accelerations of the internals joints
 * \param Tau   actuations of the internal joints (output)
 * \param f_ext External forces acting on the body in base coordinates (optional, defaults
 * to NULL)
 */

template <typename T>
void InverseDynamics(const Model &model, ModelData<T> &model_data, const Math::VectorN<T> &Q,
                     const Math::VectorN<T> &QDot, const Math::VectorN<T> &QDDot,
                     Math::VectorN<T> &Tau, std::vector<Math::SpatialVector<T> > *f_ext = NULL)
{
  // LOG << "-------- " << __func__ << " --------" << std::endl;

  // Reset the velocity of the root body
  model_data.v[0].setZero();
  model_data.a[0].set(0., 0., 0., -model.gravity[0], -model.gravity[1], -model.gravity[2]);

  for (unsigned int i = 1; i < model.mBodies.size(); i++)
  {
    unsigned int q_index = model.mJoints[i].q_index;
    unsigned int lambda = model.lambda[i];

    jcalc(model, model_data, i, Q, QDot);

    model_data.v[i] = model_data.X_lambda[i].apply(model_data.v[lambda]) + model_data.v_J[i];
    model_data.c[i] = model_data.c_J[i] + crossm(model_data.v[i], model_data.v_J[i]);

    // if(model.mJoints[i].mJointType != JointTypeCustom){
    if (model.mJoints[i].mDoFCount == 1)
    {
      model_data.a[i] = model_data.X_lambda[i].apply(model_data.a[lambda]) +
                        model_data.c[i] + model_data.S[i] * QDDot[q_index];
    }
    else if (model.mJoints[i].mDoFCount == 3)
    {
      model_data.a[i] =
          model_data.X_lambda[i].apply(model_data.a[lambda]) + model_data.c[i] +
          model_data.multdof3_S[i] *
              Vector3<T>(QDDot[q_index], QDDot[q_index + 1], QDDot[q_index + 2]);
    }
    //    }
    //    else if(model.mJoints[i].mJointType == JointTypeCustom){
    //      unsigned int k = model.mJoints[i].custom_joint_index;
    //      VectorNd customJointQDDot(model.mCustomJoints[k]->mDoFCount);
    //      for(int z=0; z<model.mCustomJoints[k]->mDoFCount; ++z){
    //        customJointQDDot[z] = QDDot[q_index+z];
    //      }
    //      model_data.a[i] =  model_data.X_lambda[i].apply(model_data.a[lambda])
    //        + model_data.c[i]
    //        + model.mCustomJoints[k]->S * customJointQDDot;
    //    }

    if (!model.mBodies[i].mIsVirtual)
    {
      model_data.f[i] = model_data.I[i] * model_data.a[i] +
                        Math::crossf(model_data.v[i], model_data.I[i] * model_data.v[i]);
    }
    else
    {
      model_data.f[i].setZero();
    }
  }

  if (f_ext != NULL)
  {
    for (unsigned int i = 1; i < model.mBodies.size(); i++)
    {
      unsigned int lambda = model.lambda[i];
      model_data.X_base[i] = model_data.X_lambda[i] * model_data.X_base[lambda];
      model_data.f[i] -= model_data.X_base[i].toMatrixAdjoint() * (*f_ext)[i];
    }
  }

  for (unsigned int i = model.mBodies.size() - 1; i > 0; i--)
  {
    // if(model.mJoints[i].mJointType != JointTypeCustom){
    if (model.mJoints[i].mDoFCount == 1)
    {
      Tau[model.mJoints[i].q_index] = model_data.S[i].dot(model_data.f[i]);
    }
    else if (model.mJoints[i].mDoFCount == 3)
    {
      Tau.template block<3, 1>(model.mJoints[i].q_index, 0) =
          model_data.multdof3_S[i].transpose() * model_data.f[i];
    }
    //    } else if (model.mJoints[i].mJointType == JointTypeCustom) {
    //      unsigned int k = model.mJoints[i].custom_joint_index;
    //      Tau.block(model.mJoints[i].q_index,0,
    //          model.mCustomJoints[k]->mDoFCount, 1)
    //        = model.mCustomJoints[k]->S.transpose() * model.f[i];
    //    }

    if (model.lambda[i] != 0)
    {
      model_data.f[model.lambda[i]] = model_data.f[model.lambda[i]] +
                                      model_data.X_lambda[i].applyTranspose(model_data.f[i]);
    }
  }
}

void InverseDynamics(Model &model, const Math::VectorNd &Q, const Math::VectorNd &QDot,
                     const Math::VectorNd &QDDot, Math::VectorNd &Tau,
                     std::vector<Math::SpatialVectord> *f_ext = NULL);


/** \brief Computes the coriolis forces
 *
 * This function computes the generalized forces from given generalized
 * states, velocities, and accelerations:
 *   \f$ \tau = M(q) \ddot{q} + N(q, \dot{q}) \f$
 *
 * \param model rigid body model
 * \param Q     state vector of the internal joints
 * \param QDot  velocity vector of the internal joints
 * \param Tau   actuations of the internal joints (output)
 */
RBDL_DLLAPI void NonlinearEffects(Model &model, ModelDatad &model_data, const Math::VectorNd &Q,
                                  const Math::VectorNd &QDot, Math::VectorNd &Tau);

/** \brief Computes forward dynamics with the Articulated Body Algorithm
 *
 * This function computes the generalized accelerations from given
 * generalized states, velocities and forces:
 *   \f$ \ddot{q} = M(q)^{-1} ( -N(q, \dot{q}) + \tau)\f$
 * It does this by using the recursive Articulated Body Algorithm that runs
 * in \f$O(n_{dof})\f$ with \f$n_{dof}\f$ being the number of joints.
 *
 * \param model rigid body model
 * \param Q     state vector of the internal joints
 * \param QDot  velocity vector of the internal joints
 * \param Tau   actuations of the internal joints
 * \param QDDot accelerations of the internal joints (output)
 * \param f_ext External forces acting on the body in base coordinates (optional, defaults
 * to NULL)
 */
template <typename T>
RBDL_DLLAPI void ForwardDynamics(const Model &model, ModelData<T> &model_data,
                                 const Math::VectorN<T> &Q, const Math::VectorN<T> &QDot,
                                 const Math::VectorN<T> &Tau, Math::VectorN<T> &QDDot,
                                 std::vector<SpatialVector<T> > *f_ext = NULL)
{
  // LOG << "-------- " << __func__ << " --------" << std::endl;

  SpatialVector<T> spatial_gravity(0., 0., 0., model.gravity[0], model.gravity[1],
                                   model.gravity[2]);

  unsigned int i = 0;

  //  LOG << "Q          = " << Q.transpose() << std::endl;
  //  LOG << "QDot       = " << QDot.transpose() << std::endl;
  //  LOG << "Tau        = " << Tau.transpose() << std::endl;
  //  LOG << "---" << std::endl;

  // Reset the velocity of the root body
  model_data.v[0].setZero();

  for (i = 1; i < model.mBodies.size(); i++)
  {
    unsigned int lambda = model.lambda[i];

    jcalc(model, model_data, i, Q, QDot);

    if (lambda != 0)
      model_data.X_base[i] = model_data.X_lambda[i] * model_data.X_base[lambda];
    else
      model_data.X_base[i] = model_data.X_lambda[i];

    model_data.v[i] = model_data.X_lambda[i].apply(model_data.v[lambda]) + model_data.v_J[i];

    /*
       LOG << "X_J (" << i << "):" << std::endl << X_J << std::endl;
       LOG << "v_J (" << i << "):" << std::endl << v_J << std::endl;
       LOG << "v_lambda" << i << ":" << std::endl << model_data.v.at(lambda) << std::endl;
       LOG << "X_base (" << i << "):" << std::endl << model_data.X_base[i] << std::endl;
       LOG << "X_lambda (" << i << "):" << std::endl << model_data.X_lambda[i] <<
       std::endl;
       LOG << "SpatialVelocity (" << i << "): " << model_data.v[i] << std::endl;
       */
    model_data.c[i] = model_data.c_J[i] + crossm(model_data.v[i], model_data.v_J[i]);
    model_data.I[i].setSpatialMatrix(model_data.IA[i]);

    model_data.pA[i] = crossf(model_data.v[i], model_data.I[i] * model_data.v[i]);

    if (f_ext != NULL)
    {
      // LOG << "External force (" << i << ") = " <<
      // model_data.X_base[i].toMatrixAdjoint() * (*f_ext)[i] << std::endl;
      model_data.pA[i] -= model_data.X_base[i].toMatrixAdjoint() * (*f_ext)[i];
    }
  }

  // ClearLogOutput();

  // LOG << "--- first loop ---" << std::endl;

  for (i = model.mBodies.size() - 1; i > 0; i--)
  {
    unsigned int q_index = model.mJoints[i].q_index;

    if (model.mJoints[i].mDoFCount == 1)
    {
      //&& model.mJoints[i].mJointType != JointTypeCustom) {

      model_data.U[i] = model_data.IA[i] * model_data.S[i];
      model_data.d[i] = model_data.S[i].dot(model_data.U[i]);
      model_data.u[i] = Tau[q_index] - model_data.S[i].dot(model_data.pA[i]);
      //      LOG << "u[" << i << "] = " << model_data.U[i] << std::endl;

      unsigned int lambda = model.lambda[i];
      if (lambda != 0)
      {
        SpatialMatrix<T> Ia =
            model_data.IA[i] - model_data.U[i] * (model_data.U[i] / model_data.d[i]).transpose();

        SpatialVector<T> pa = model_data.pA[i] + Ia * model_data.c[i] +
                              model_data.U[i] * model_data.u[i] / model_data.d[i];

#ifdef EIGEN_CORE_H
        model_data.IA[lambda].noalias() += model_data.X_lambda[i].toMatrixTranspose() *
                                           Ia * model_data.X_lambda[i].toMatrix();
        model_data.pA[lambda].noalias() += model_data.X_lambda[i].applyTranspose(pa);
#else
        model_data.IA[lambda] += model_data.X_lambda[i].toMatrixTranspose() * Ia *
                                 model_data.X_lambda[i].toMatrix();

        model_data.pA[lambda] += model_data.X_lambda[i].applyTranspose(pa);
#endif
        // LOG << "pA[" << lambda << "] = "
        //     << model_data.pA[lambda].transpose() << std::endl;
      }
    }
    else if (model.mJoints[i].mDoFCount == 3)
    {
      //&& model.mJoints[i].mJointType != JointTypeCustom) {
      model_data.multdof3_U[i] = model_data.IA[i] * model_data.multdof3_S[i];
#ifdef EIGEN_CORE_H
      model_data.multdof3_Dinv[i] =
          (model_data.multdof3_S[i].transpose() * model_data.multdof3_U[i]).inverse().eval();
#else
      model_data.multdof3_Dinv[i] =
          (model_data.multdof3_S[i].template transpose() * model_data.multdof3_u[i]).inverse();
#endif
      Vector3<T> tau_temp(Tau[q_index], Tau[q_index + 1], Tau[q_index + 2]);
      model_data.multdof3_u[i] =
          tau_temp - model_data.multdof3_S[i].transpose() * model_data.pA[i];

      // LOG << "multdof3_u[" << i << "] = "
      //                      << model_data.multdof3_u[i].transpose() << std::endl;
      unsigned int lambda = model.lambda[i];
      if (lambda != 0)
      {
        SpatialMatrix<T> Ia = model_data.IA[i] -
                              model_data.multdof3_U[i] * model_data.multdof3_Dinv[i] *
                                  model_data.multdof3_U[i].transpose();
        SpatialVector<T> pa = model_data.pA[i] + Ia * model_data.c[i] +
                              model_data.multdof3_U[i] * model_data.multdof3_Dinv[i] *
                                  model_data.multdof3_u[i];
#ifdef EIGEN_CORE_H
        model_data.IA[lambda].noalias() += model_data.X_lambda[i].toMatrixTranspose() *
                                           Ia * model_data.X_lambda[i].toMatrix();

        model_data.pA[lambda].noalias() += model_data.X_lambda[i].applyTranspose(pa);
#else
        model_data.IA[lambda] += model_data.X_lambda[i].toMatrixTranspose() * Ia *
                                 model_data.X_lambda[i].toMatrix();

        model_data.pA[lambda] += model_data.X_lambda[i].applyTranspose(pa);
#endif
        //        LOG << "pA[" << lambda << "] = "
        //            << model_data.pA[lambda].transpose()
        //            << std::endl;
      }
    } /* else if (model.mJoints[i].mJointType == JointTypeCustom) {
       unsigned int kI   = model.mJoints[i].custom_joint_index;
       unsigned int dofI = model.mCustomJoints[kI]->mDoFCount;
       model.mCustomJoints[kI]->U =
         model_data.IA[i] * model.mCustomJoints[kI]->S;

 #ifdef EIGEN_CORE_H
       model.mCustomJoints[kI]->Dinv
         = (model.mCustomJoints[kI]->S.transpose()
             * model.mCustomJoints[kI]->U).inverse().eval();
 #else
       model.mCustomJoints[kI]->Dinv
         = (model.mCustomJoints[kI]->S.transpose()
             * model.mCustomJoints[kI]->U).inverse();
 #endif
       VectorNd tau_temp(dofI);
       for(int z=0;z<dofI;++z){
         tau_temp(z) = Tau[q_index+z];
       }
       model.mCustomJoints[kI]->u = tau_temp
         - model.mCustomJoints[kI]->S.transpose() * model_data.pA[i];

       //      LOG << "multdof3_u[" << i << "] = "
       //      << model_data.multdof3_u[i].transpose() << std::endl;
       unsigned int lambda = model.lambda[i];
       if (lambda != 0) {
         SpatialMatrixd Ia = model_data.IA[i]
           - (model.mCustomJoints[kI]->U
               * model.mCustomJoints[kI]->Dinv
               * model.mCustomJoints[kI]->U.transpose());
         SpatialVectord pa =  model_data.pA[i]
           + Ia * model_data.c[i]
           + (model.mCustomJoints[kI]->U
               * model.mCustomJoints[kI]->Dinv
               * model.mCustomJoints[kI]->u);

 #ifdef EIGEN_CORE_H
         model_data.IA[lambda].noalias() += model_data.X_lambda[i].toMatrixTranspose()
           * Ia
           * model_data.X_lambda[i].toMatrix();
         model_data.pA[lambda].noalias() += model_data.X_lambda[i].applyTranspose(pa);
 #else
         model_data.IA[lambda] += model_data.X_lambda[i].toMatrixTranspose()
           * Ia
           * model_data.X_lambda[i].toMatrix();
         model_data.pA[lambda] += model_data.X_lambda[i].applyTranspose(pa);
 #endif
         LOG << "pA[" << lambda << "] = "
           << model_data.pA[lambda].transpose()
           << std::endl;
       }
     }
     */
  }

  //  ClearLogOutput();

  model_data.a[0] = spatial_gravity * -1.;

  for (i = 1; i < model.mBodies.size(); i++)
  {
    unsigned int q_index = model.mJoints[i].q_index;
    unsigned int lambda = model.lambda[i];
    SpatialTransform<T> X_lambda = model_data.X_lambda[i];

    model_data.a[i] = X_lambda.apply(model_data.a[lambda]) + model_data.c[i];
    // LOG << "a'[" << i << "] = " << model_data.a[i].transpose() << std::endl;

    if (model.mJoints[i].mDoFCount == 1)
    {
      //&& model.mJoints[i].mJointType != JointTypeCustom) {
      QDDot[q_index] = (T(1.) / model_data.d[i]) *
                       (model_data.u[i] - model_data.U[i].dot(model_data.a[i]));
      model_data.a[i] = model_data.a[i] + model_data.S[i] * QDDot[q_index];
    }
    else if (model.mJoints[i].mDoFCount == 3)
    {
      //&& model.mJoints[i].mJointType != JointTypeCustom) {
      Vector3<T> qdd_temp =
          model_data.multdof3_Dinv[i] *
          (model_data.multdof3_u[i] -
           (model_data.multdof3_U[i].transpose() * model_data.a[i]));
      QDDot[q_index] = qdd_temp[0];
      QDDot[q_index + 1] = qdd_temp[1];
      QDDot[q_index + 2] = qdd_temp[2];
      model_data.a[i] = model_data.a[i] + model_data.multdof3_S[i] * qdd_temp;
    } /* else if (model.mJoints[i].mJointType == JointTypeCustom) {
       unsigned int kI = model.mJoints[i].custom_joint_index;
       unsigned int dofI=model.mCustomJoints[kI]->mDoFCount;

       VectorNd qdd_temp = model.mCustomJoints[kI]->Dinv
         * (  model.mCustomJoints[kI]->u
             - model.mCustomJoints[kI]->U.transpose()
             * model_data.a[i]);

       for(int z=0; z<dofI; ++z){
         QDDot[q_index+z] = qdd_temp[z];
       }

       model_data.a[i] = model_data.a[i]
         + model.mCustomJoints[kI]->S * qdd_temp;
     } */
  }

  //  LOG << "QDDot = " << QDDot.transpose() << std::endl;
}

void ForwardDynamics(Model &model, const Math::VectorNd &Q, const Math::VectorNd &QDot,
                     const Math::VectorNd &Tau, Math::VectorNd &QDDot,
                     std::vector<SpatialVectord> *f_ext = NULL);


/** \brief Computes forward dynamics by building and solving the full Lagrangian equation
 *
 * This method builds and solves the linear system
 * \f[ 	H \ddot{q} = -C + \tau	\f]
 * for \f$\ddot{q}\f$ where \f$H\f$ is the joint space inertia matrix
 * computed with the CompositeRigidBodyAlgorithm(), \f$C\f$ the bias
 * force (sometimes called "non-linear effects").
 *
 * \param model rigid body model
 * \param Q     state vector of the internal joints
 * \param QDot  velocity vector of the internal joints
 * \param Tau   actuations of the internal joints
 * \param QDDot accelerations of the internal joints (output)
 * \param linear_solver specification which method should be used for solving the linear
 * system
 * \param f_ext External forces acting on the body in base coordinates (optional, defaults
 * to NULL)
 * \param H     preallocated workspace area for the joint space inertia matrix of size
 * dof_count x dof_count (optional, defaults to NULL and allocates temporary matrix)
 * \param C     preallocated workspace area for the right hand side vector of size
 * dof_count x 1 (optional, defaults to NULL and allocates temporary vector)
 */

template <typename T>
void CompositeRigidBodyAlgorithm(const Model &model, ModelData<T> &model_data,
                                 const VectorN<T> &Q, MatrixN<T> &H,
                                 bool update_kinematics = true)
{
  // LOG << "-------- " << __func__ << " --------" << std::endl;

  H.setZero();
  assert(H.rows() == model.dof_count && H.cols() == model.dof_count);

  for (unsigned int i = 1; i < model.mBodies.size(); i++)
  {
    if (update_kinematics)
    {
      jcalc_X_lambda_S(model, model_data, i, Q);
    }
    model_data.Ic[i] = model_data.I[i];
  }

  for (unsigned int i = model.mBodies.size() - 1; i > 0; i--)
  {
    if (model.lambda[i] != 0)
    {
      model_data.Ic[model.lambda[i]] = model_data.Ic[model.lambda[i]] +
                                       model_data.X_lambda[i].applyTranspose(model_data.Ic[i]);
    }

    unsigned int dof_index_i = model.mJoints[i].q_index;

    if (model.mJoints[i].mDoFCount == 1)
    {
      //&& model.mJoints[i].mJointType != JointTypeCustom) {

      SpatialVector<T> F = model_data.Ic[i] * model_data.S[i];
      H(dof_index_i, dof_index_i) = model_data.S[i].dot(F);

      unsigned int j = i;
      unsigned int dof_index_j = dof_index_i;

      while (model.lambda[j] != 0)
      {
        F = model_data.X_lambda[j].applyTranspose(F);
        j = model.lambda[j];
        dof_index_j = model.mJoints[j].q_index;

        // if(model.mJoints[j].mJointType != JointTypeCustom) {
        if (model.mJoints[j].mDoFCount == 1)
        {
          H(dof_index_i, dof_index_j) = F.dot(model_data.S[j]);
          H(dof_index_j, dof_index_i) = H(dof_index_i, dof_index_j);
        }
        else if (model.mJoints[j].mDoFCount == 3)
        {
          Vector3<T> H_temp2 = (F.transpose() * model_data.multdof3_S[j]).transpose();
          //          LOG << F.transpose() << std::endl
          //              << model_data.multdof3_S[j] << std::endl;
          //          LOG << H_temp2.transpose() << std::endl;

          H.template block<1, 3>(dof_index_i, dof_index_j) = H_temp2.transpose();
          H.template block<3, 1>(dof_index_j, dof_index_i) = H_temp2;
        }
        //        } else if (model.mJoints[j].mJointType == JointTypeCustom){
        //          unsigned int k      = model.mJoints[j].custom_joint_index;
        //          unsigned int dof    = model.mCustomJoints[k]->mDoFCount;
        //          VectorNd H_temp2    =
        //            (F.transpose() * model.mCustomJoints[k]->S).transpose();

        //          LOG << F.transpose()
        //            << std::endl
        //            << model.mCustomJoints[j]->S << std::endl;

        //          LOG << H_temp2.transpose() << std::endl;

        //          H.block(dof_index_i,dof_index_j,1,dof) = H_temp2.transpose();
        //          H.block(dof_index_j,dof_index_i,dof,1) = H_temp2;
        //        }
      }
    }
    else if (model.mJoints[i].mDoFCount == 3)
    {
      //    && model.mJoints[i].mJointType != JointTypeCustom) {
      Matrix63<T> F_63 = model_data.Ic[i].toMatrix() * model_data.multdof3_S[i];
      H.template block<3, 3>(dof_index_i, dof_index_i) =
          model_data.multdof3_S[i].transpose() * F_63;

      unsigned int j = i;
      unsigned int dof_index_j = dof_index_i;

      while (model.lambda[j] != 0)
      {
        F_63 = model_data.X_lambda[j].toMatrixTranspose() * (F_63);
        j = model.lambda[j];
        dof_index_j = model.mJoints[j].q_index;

        //  if(model.mJoints[j].mJointType != JointTypeCustom){
        if (model.mJoints[j].mDoFCount == 1)
        {
          Vector3<T> H_temp2 = F_63.transpose() * (model_data.S[j]);

          H.template block<3, 1>(dof_index_i, dof_index_j) = H_temp2;
          H.template block<1, 3>(dof_index_j, dof_index_i) = H_temp2.transpose();
        }
        else if (model.mJoints[j].mDoFCount == 3)
        {
          Matrix3<T> H_temp2 = F_63.transpose() * (model_data.multdof3_S[j]);

          H.template block<3, 3>(dof_index_i, dof_index_j) = H_temp2;
          H.template block<3, 3>(dof_index_j, dof_index_i) = H_temp2.transpose();
        }
        //        } else if (model.mJoints[j].mJointType == JointTypeCustom){
        //          unsigned int k = model.mJoints[j].custom_joint_index;
        //          unsigned int dof = model.mCustomJoints[k]->mDoFCount;

        //          MatrixNd H_temp2 = F_63.transpose() * (model.mCustomJoints[k]->S);

        //          H.block(dof_index_i,dof_index_j,3,dof) = H_temp2;
        //          H.block(dof_index_j,dof_index_i,dof,3) = H_temp2.transpose();
        //        }
      }
    }
    //  else if (model.mJoints[i].mJointType == JointTypeCustom) {
    //      unsigned int kI = model.mJoints[i].custom_joint_index;
    //      unsigned int dofI = model.mCustomJoints[kI]->mDoFCount;

    //      MatrixNd F_Nd = model_data.Ic[i].toMatrix()
    //        * model.mCustomJoints[kI]->S;

    //      H.block(dof_index_i, dof_index_i,dofI,dofI)
    //        = model.mCustomJoints[kI]->S.transpose() * F_Nd;

    //      unsigned int j = i;
    //      unsigned int dof_index_j = dof_index_i;

    //      while (model.lambda[j] != 0) {
    //        F_Nd = model_data.X_lambda[j].toMatrixTranspose() * (F_Nd);
    //        j = model.lambda[j];
    //        dof_index_j = model.mJoints[j].q_index;

    //        if(model.mJoints[j].mJointType != JointTypeCustom){
    //          if (model.mJoints[j].mDoFCount == 1) {
    //            MatrixNd H_temp2 = F_Nd.transpose() * (model_data.S[j]);
    //            H.block(   dof_index_i,  dof_index_j,
    //                H_temp2.rows(),H_temp2.cols()) = H_temp2;
    //            H.block(dof_index_j,dof_index_i,
    //                H_temp2.cols(),H_temp2.rows()) = H_temp2.transpose();
    //          } else if (model.mJoints[j].mDoFCount == 3) {
    //            MatrixNd H_temp2 = F_Nd.transpose() * (model_data.multdof3_S[j]);
    //            H.block(dof_index_i,   dof_index_j,
    //                H_temp2.rows(),H_temp2.cols()) = H_temp2;
    //            H.block(dof_index_j,   dof_index_i,
    //                H_temp2.cols(),H_temp2.rows()) = H_temp2.transpose();
    //          }
    //        } else if (model.mJoints[j].mJointType == JointTypeCustom){
    //          unsigned int k   = model.mJoints[j].custom_joint_index;
    //          unsigned int dof = model.mCustomJoints[k]->mDoFCount;

    //          MatrixNd H_temp2 = F_Nd.transpose() * (model.mCustomJoints[k]->S);

    //          H.block(dof_index_i,dof_index_j,3,dof) = H_temp2;
    //          H.block(dof_index_j,dof_index_i,dof,3) = H_temp2.transpose();
    //        }
    //      }
    //    }
  }
}

void CompositeRigidBodyAlgorithm(Model &model, const VectorNd &Q, MatrixNd &H,
                                 bool update_kinematics = true);

RBDL_DLLAPI void ForwardDynamicsLagrangian(
    Model &model, ModelDatad &model_data, const VectorNd &Q, const VectorNd &QDot,
    const VectorNd &Tau, VectorNd &QDDot,
    Math::LinearSolver linear_solver = Math::LinearSolverColPivHouseholderQR,
    std::vector<Math::SpatialVectord> *f_ext = NULL, Math::MatrixNd *H = NULL,
    Math::VectorNd *C = NULL);

/** \brief Computes the effect of multiplying the inverse of the joint
 * space inertia matrix with a vector in linear time.
 *
 * \param model rigid body model
 * \param Q     state vector of the generalized positions
 * \param Tau   the vector that should be multiplied with the inverse of
 *              the joint space inertia matrix
 * \param QDDot vector where the result will be stored
 * \param update_kinematics whether the kinematics should be updated (safer, but at a
 * higher computational cost)
 *
 * This function uses a reduced version of the Articulated %Body Algorithm
 * to compute
 *
 *   \f$ \ddot{q} = M(q)^{-1} ( -N(q, \dot{q}) + \tau)\f$
 *
 * in \f$O(n_{\textit{dof}}\f$) time.
 *
 * \note When calling this function repeatedly for the same values of Q make sure
 * to set the last parameter to false as this avoids expensive
 * recomputations of transformations and articulated body inertias.
 */
RBDL_DLLAPI void CalcMInvTimesTau(Model &model, ModelDatad &model_data,
                                  const Math::VectorNd &Q, const Math::VectorNd &Tau,
                                  Math::VectorNd &QDDot, bool update_kinematics = true);

/** @} */
}

/* RBDL_DYNAMICS_H */
#endif
