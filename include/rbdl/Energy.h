#ifndef _ENERGY_H
#define _ENERGY_H

#include <rbdl/rbdl_math.h>
#include <assert.h>
#include <iostream>
#include "rbdl/Logging.h"

namespace RigidBodyDynamics {

  /// @todo: Document
  Math::SpatialVector CalcEnergy_ineficient (
      Model &model,
      const Math::VectorNd &Q,
      const Math::VectorNd &QDot,
      bool update_kinematics = false,
      unsigned int method = 0
      );

  void Energy_matlab(
      Model &model,
      const Math::VectorNd &Q,
      const Math::VectorNd &QDot,
      Math::SpatialMatrix &Itot,
      Math::SpatialVector htot,
      double &total_KE,
      Math::Vector3d &vcm,
      bool update_kinematics = false
      );

  void CCM_CCRBI(Model &model,
                 const Math::VectorNd &Q,
                 const Math::VectorNd &QDot,
                 Math::SpatialMatrix &I,
                 bool update_kinematics = false);


  void CCM_CCRBI(Model &model,
                 const Math::VectorNd &Q,
                 const Math::VectorNd &QDot,
                 Math::SpatialMatrix &I,
                 Math::MatrixNd &AG,
                 Math::SpatialVector &h,
                 bool update_kinematics = false);

  void CCM_CCRBI_Jacobian_com(Model &model,
                              const Math::VectorNd &Q,
                              Math::MatrixNd &AG,
                              bool update_kinematics = false);



  /*
  void systemMomentumMatrix(Model &model,
                            const Math::VectorNd &Q,
                            const Math::VectorNd &QDot,
                            Math::MatrixNd &A,
                            Math::SpatialVector &h,
                            bool update_kinematics = false);

  void massMatrixFromCCRBI(Model &model,
                           const Math::VectorNd &Q,
                           const Math::VectorNd &QDot,
                           Math::MatrixNd &H );
  */
}


#endif

