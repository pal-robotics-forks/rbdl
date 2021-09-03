#ifndef RBDL_ENERGY_H
#define RBDL_ENERGY_H

#include <assert.h>
#include <iostream>

#include "rbdl/rbdl_math.h"
#include "rbdl/rbdl_mathutils.h"
#include "rbdl/Model.h"
#include <rbdl/Kinematics.h>
#include "rbdl/Logging.h"

namespace RigidBodyDynamics
{
using namespace Math;

double CalcTotalMass(const Model &model);

double CalcPartialMass(const Model &model, const std::string &body_1, const std::string &body_2);

template <typename T>
Vector3<T> CalcCOM(const Model &model, ModelData<T> &model_data, const VectorN<T> &Q,
                   bool update_kinematics = true)
{
  // update the Kinematics if necessary
  if (update_kinematics)
  {
    UpdateKinematicsCustom<T>(model, model_data, &Q, NULL, NULL);
  }
  Vector3<T> com;
  com.setZero();
  T total_mass = 0;

  for (size_t i = 1; i < model.mBodies.size(); ++i)
  {
    int body_id = i;  // rbdl_model_.GetBodyId(link_names_[i].c_str());
    Vector3<T> link_com;
    link_com = CalcBodyToBaseCoordinates<T>(model, model_data, Q, body_id,
                                            model.mBodies[body_id].mCenterOfMass.cast<T>(), false);

    com += model.mBodies[body_id].mMass * link_com;
    total_mass += model.mBodies[body_id].mMass;
  }
  return com / total_mass;
}

Math::Vector3d CalcCOM(Model &model, const Math::VectorNd &Q, bool update_kinematics = true);

Math::Vector3d CalcCOMVelocity(const Model &model, ModelDatad &model_data,
                               const Math::VectorNd &Q, const Math::VectorNd &QDot,
                               bool update_kinematics = true);

Math::Vector3d CalcCOMVelocity(Model &model, const Math::VectorNd &Q,
                               const Math::VectorNd &QDot, bool update_kinematics = true);

Math::Vector3d CalcCOMAcceleartion(const Model &model, ModelDatad &model_data,
                                   const Math::VectorNd &Q, const Math::VectorNd &QDot,
                                   const Math::VectorNd &QDDot, bool update_kinematics = true);

Math::Vector3d CalcCOMAcceleartion(Model &model, const Math::VectorNd &Q,
                                   const Math::VectorNd &QDot, const Math::VectorNd &QDDot,
                                   bool update_kinematics = true);

Math::Vector3d CalcCOMAccelerationBias(const Model &model, ModelDatad &model_data,
                                       const Math::VectorNd &Q, const Math::VectorNd &QDot,
                                       bool update_kinematics = true);

Math::Vector3d CalcCOMAccelerationBias(Model &model, const Math::VectorNd &Q,
                                       const Math::VectorNd &QDot,
                                       bool update_kinematics = true);

void CalcCOMJacobian(const Model &model, ModelDatad &model_data, const Math::VectorNd &Q,
                     Math::MatrixNd &COMJ, bool update_kinematics = true);

void CalcCOMJacobian(Model &model, const Math::VectorNd &Q, Math::MatrixNd &COMJ,
                     bool update_kinematics = true);

void CalcCOMJacobian_inefficient(const Model &model, ModelDatad &model_data,
                                 const Math::VectorNd &Q, Math::MatrixNd &COMJ,
                                 bool update_kinematics = true);


//  void CCM_CCRBI_Jacobian_com(Model &model,
//                              const Math::VectorNd &Q,
//                              Math::MatrixNd &AG,
//                              bool update_kinematics = true);

/// @todo: Document
Math::SpatialVectord CalcEnergy_inefficient(const Model &model, ModelDatad &model_data,
                                            const Math::VectorNd &Q, const Math::VectorNd &QDot,
                                            bool update_kinematics = true,
                                            unsigned int method = 0);

Math::Matrix3d calcGlobalInertiaTensorFromCOM(Model &model, const Math::VectorNd &Q, bool update_kinematics = true);
}

#endif
