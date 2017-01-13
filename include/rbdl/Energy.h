#ifndef RBDL_ENERGY_H
#define RBDL_ENERGY_H

#include <assert.h>
#include <iostream>

#include "rbdl/rbdl_math.h"
#include "rbdl/rbdl_mathutils.h"
#include "rbdl/Model.h"

#include "rbdl/Logging.h"

namespace RigidBodyDynamics {

  double CalcTotalMass(Model &model);

  Math::Vector3d CalcCOM(Model &model,
                         const Math::VectorNd &Q,
                         bool update_kinematics = true);


  Math::Vector3d CalcCOMVelocity(Model &model,
                                 const Math::VectorNd &Q,
                                 const Math::VectorNd &QDot,
                                 bool update_kinematics = true);

  Math::Vector3d CalcCOMAcceleartion(Model &model,
                                     const Math::VectorNd &Q,
                                     const Math::VectorNd &QDot,
                                     const Math::VectorNd &QDDot,
                                     bool update_kinematics = true);

  Math::Vector3d CalcCOMAccelerationBias(Model &model,
                                         const Math::VectorNd &Q,
                                         const Math::VectorNd &QDot,
                                         bool update_kinematics = true);

  void CalcCOMJacobian (
      Model &model,
      const Math::VectorNd &Q,
      Math::MatrixNd &COMJ,
      bool update_kinematics = true);


  void CalcCOMJacobian_inefficient (
      Model &model,
      const Math::VectorNd &Q,
      Math::MatrixNd &COMJ,
      bool update_kinematics = true
      );


//  void CCM_CCRBI_Jacobian_com(Model &model,
//                              const Math::VectorNd &Q,
//                              Math::MatrixNd &AG,
//                              bool update_kinematics = true);

  /// @todo: Document
  Math::SpatialVectord CalcEnergy_inefficient (
      Model &model,
      const Math::VectorNd &Q,
      const Math::VectorNd &QDot,
      bool update_kinematics = true,
      unsigned int method = 0
      );

}

#endif
