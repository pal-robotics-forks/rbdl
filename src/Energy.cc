#include <iostream>
#include <limits>
#include <cstring>
#include <assert.h>

#include "rbdl/rbdl_mathutils.h"
#include "rbdl/Logging.h"

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"
#include "rbdl/Energy.h"

/// @todo change all the Eigen math types to rbdl math types

namespace RigidBodyDynamics {

  using namespace Math;
  
  double CalcTotalMass(Model &model){
    double totalMass = 0.;
    for(size_t i=0; i<model.mBodies.size(); ++i){
      totalMass += model.mBodies[i].mMass;
    }
    return totalMass;
  }

  Math::SpatialVector CalcMoment(
      Model &model,
      const Math::VectorNd &Q,
      const Math::VectorNd &QDot,
      bool update_kinematics) {
    
    if (update_kinematics) {
      UpdateKinematicsCustom(model, &Q, &QDot, NULL);
    }

    Math::SpatialVector spatial_moment;
    spatial_moment.setZero();

    for(unsigned int i=1; i<model.mBodies.size(); ++i){
      unsigned int body_id = i;//rbdl_model_.GetBodyId(link_names_[i].c_str());

        /// X'*I*X*(X^-1*v)
        spatial_moment += model.X_base[body_id].toMatrix().transpose()*
            model.mBodies[body_id].mSpatialInertia*
            model.X_base[body_id].toMatrix()*
            model.X_base[body_id].toMatrix().inverse()*model.v[body_id];

    }
    
    return spatial_moment;
  }
  
  Math::SpatialVector CalcMomentCOM (
      Model &model,
      const Math::VectorNd &Q,
      const Math::VectorNd &QDot,
      bool update_kinematics
      ){

     Math::Vector3d com = CalCOM(model, Q, update_kinematics);
     //return  Math::Xtrans(com).applyAdjoint(CalcMoment(model, Q, QDot, update_kinematics));
     return  Math::Xtrans(com).toMatrix().transpose()*CalcMoment(model, Q, QDot, update_kinematics);
   }
   
  Math::SpatialVector CalcEnergy_ineficient(Model &model,
      const Math::VectorNd &Q,
      const Math::VectorNd &QDot,
      bool update_kinematics,
      unsigned int method
      ) {
    
    if (update_kinematics) {
      //Compute energy & COM & mass of the robot;
      RigidBodyDynamics::UpdateKinematicsCustom(model, &Q, &QDot, NULL);
    }

    SpatialVector spatial_moment = RigidBodyDynamics::Math::SpatialVectorZero;

    /// @todo: Is there a problem with the floating base when calculation COM, is the base_link COM getting stuck at the origin
    /// of the inertial frame?

    /// @warning We start with body 1 becouse the first link is the root world link which seems to have 1 kg mass
    for(size_t i=1; i<model.mBodies.size(); ++i){
      unsigned int body_id = i;//rbdl_model_.GetBodyId(link_names_[i].c_str());

      if(method == 0){
        /// X'*I*X*(X^-1*v)
        spatial_moment += model.X_base[body_id].toMatrix().transpose()*
            model.mBodies[body_id].mSpatialInertia*
            model.X_base[body_id].toMatrix()*
            model.X_base[body_id].toMatrix().inverse()*model.v[body_id];
      }
      else if(method == 1){
        //Optimization using spatial operators, in the book there is a further optimization to reduce operations
        //    spatial_moment += model.X_base[body_id].apply(model.mBodies[body_id].mSpatialInertia*
        //        (model.X_base[body_id].applyInverse(model.v[body_id])));

        spatial_moment += model.X_base[body_id].toMatrix().transpose()*
            model.mBodies[body_id].mSpatialInertia*
            model.X_base[body_id].toMatrix()*
            (model.X_base[body_id].applyInverse(model.v[body_id]));

      }
      else{
        //Spatial momentum vectors have the same properties as spatial force vectors
        // Optimization using spatial operators, in the book there is a further optimization to reduce operations
        spatial_moment += model.X_base[body_id].applyTranspose(model.mBodies[body_id].mSpatialInertia*(model.v[body_id]));
      }

      //std::cerr<<"Spatial moment: "<<i<<"  is: "<<spatial_moment.transpose()<<std::endl;
    }
    
    return spatial_moment;
  }
  
  // EnerMo  calculate energy, momentum and related quantities
  // EnerMo(robot,q,qd)  returns a structure containing the fields KE, PE,
  // htot, Itot, mass, cm and vcm.  These fields contain the kinetic and
  // potential energies of the whole system, the total spatial momentum, the
  // total spatial inertia, total mass, position of centre of mass, and the
  // linear velocity of centre of mass, respectively.  Vector quantities are
  // expressed in base coordinates.  PE is defined to be zero when cm is
  // zero.
  void Energy_matlab(
      Model &model,
      const VectorNd &Q,
      const VectorNd &QDot,
      SpatialMatrix &Itot,
      SpatialVector htot,
      double &total_KE,
      Eigen::Vector3d &vcm,
      bool update_kinematics = false
      ) {
    
    if (update_kinematics) {
      UpdateKinematicsCustom(model, &Q, &QDot, NULL);
    }

    Itot.setZero();
    htot.setZero();

    std::vector<SpatialMatrix> Ic(model.mBodies.size());
    std::vector<SpatialVector> hc(model.mBodies.size());
    Eigen::VectorXd KE(model.mBodies.size());

    for(unsigned int i=1; i<model.mBodies.size(); ++i){
      Ic[i] = model.mBodies[i].mSpatialInertia;
      hc[i] = Ic[i] * model.v[i];
      KE[i] = 0.5 * model.v[i].transpose() * hc[i];
    }

    for (int i = model.mBodies.size() -1; i >= 1; --i){
      unsigned int lambda = model.lambda[i];
      if( lambda != 0){
        /// @todo optimize using spatial inertia operators
        Ic[lambda] = Ic[lambda] + model.X_base[i].toMatrix().transpose()*Ic[i]*model.X_base[i].toMatrix();
        hc[lambda] = hc[lambda] + model.X_base[i].toMatrix().transpose()*hc[i];
      }
      else{
        Itot = Itot + model.X_base[i].toMatrix().transpose()*Ic[i]*model.X_base[i].toMatrix();
        htot = htot + model.X_base[i].toMatrix().transpose()*hc[i];
      }
    }

    Eigen::Vector3d g = model.gravity;
    Eigen::Vector3d h = htot.segment(3,3); //3D linear momentum

    double mass; Eigen::Vector3d cm; SpatialRigidBodyInertia mcI;
    mcI.createFromMatrix(Itot);
    mass = mcI.m;
    cm = mcI.h/mcI.m;

    total_KE = KE.sum(); //std::accumulate
    //double total_PE = mass * (cm.dot(g));
    vcm = h / mass; // COM velocity
  }

  void CCM_CCRBI(Model &model,
                 const VectorNd &Q,
                 const VectorNd &QDot,
                 SpatialMatrix &I,
                 bool update_kinematics){

    if (update_kinematics) {
      UpdateKinematicsCustom(model, &Q, &QDot, NULL);
    }

    std::vector<SpatialMatrix> Ic(model.mBodies.size());
    std::vector<SpatialVector> hc(model.mBodies.size());
    Eigen::VectorXd KE(model.mBodies.size());

    for(unsigned int i=1; i<model.mBodies.size(); ++i){
      Ic[i] = model.mBodies[i].mSpatialInertia;
    }
    for (int i = model.mBodies.size() -1; i >= 1; --i){
      unsigned int lambda = model.lambda[i];
      if( lambda != 0){
        /// @todo optimize using spatial inertia operators
        Ic[lambda] = Ic[lambda] + model.X_lambda[i].toMatrix().transpose()*Ic[i]*model.X_lambda[i].toMatrix();
      }
    }
    I = Ic[1];

  }

  void CCM_CCRBI(Model &model,
                 const VectorNd &Q,
                 const VectorNd &QDot,
                 SpatialMatrix &I,
                 MatrixNd &AG,
                 SpatialVector &h,
                 bool update_kinematics){

    assert(AG.cols() == model.dof_count);
    assert(AG.rows() == 6);

    if (update_kinematics) {
      UpdateKinematicsCustom(model, &Q, &QDot, NULL);
    }

    std::vector<SpatialMatrix> Ic(model.mBodies.size());
    h.setZero();

    for(unsigned int i=1; i<model.mBodies.size(); ++i){
      Ic[i] = model.mBodies[i].mSpatialInertia;
    }
    for (int i = model.mBodies.size() -1; i >= 1; --i){
      unsigned int lambda = model.lambda[i];
      if( lambda != 0){
        /// @todo optimize using spatial inertia operators
        Ic[lambda] = Ic[lambda] + model.X_lambda[i].toMatrix().transpose()*Ic[i]*model.X_lambda[i].toMatrix();
      }
    }
    I = Ic[1];

    for(unsigned int i=1; i<model.mBodies.size(); ++i){
      SpatialVector AGi = model.X_base[i].toMatrix().transpose()*Ic[i]* model.mJoints[i].mJointAxes[0];
      for(unsigned int j=0; j<6; ++j){
        AG(j, i-1) = AGi(j);
      }
      h = h + AGi*QDot(i-1);
    }

  }
  
  /**
    Using a methodology similar to the COM of mass jacobian
  */
  void CCM_CCRBI_Jacobian_com(Model &model,
                              const VectorNd &Q,
                              MatrixNd &AG,
                              bool update_kinematics){
    /// @todo: Reuse computations of the screws and think to formulate everything directly at the com coordinate frame
    assert(AG.cols() == model.dof_count);
    assert(AG.rows() == 6);
    AG.setZero();
    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematicsCustom (model, &Q, NULL, NULL);
    }

    Vector3d com = CalCOM(model, Q, update_kinematics);

    Eigen::MatrixXd jacobian_link(6, model.dof_count);
    for (unsigned int i = 1; i < model.mBodies.size(); ++i) {
      CalcPoseJacobian(model, Q, i, model.mBodies[i].mCenterOfMass, jacobian_link, false);
      /// @todo is Ic initialized?
      SpatialMatrix Icom = model.X_base[i].toMatrix().transpose()*model.mBodies[i].mSpatialInertia*model.X_base[i].toMatrix();//Inertia at the base mass frame
      AG += Icom*jacobian_link;
    }
    AG = Xtrans (com).toMatrix().transpose()*AG; //Change to COM coordinates
  }

  /*
  void systemMomentumMatrix(Model &model,
                            const VectorNd &Q,
                            const VectorNd &QDot,
                            Eigen::MatrixXd &A,
                            SpatialVector &h,
                            bool update_kinematics){

    if (update_kinematics) {
      UpdateKinematicsCustom(model, &Q, &QDot, NULL);
    }

    Eigen::MatrixXd J(6*(model.mBodies.size() - 1), model.dof_count);
    Eigen::MatrixXd I(6*(model.mBodies.size() - 1), 6);

    unsigned int row_index = 0;
    Eigen::MatrixXd Jac_temp(6, model.dof_count);
    for(unsigned int i=1; i<model.mBodies.size() - 1; ++i){
      Eigen::Vector3d zero;
      zero.setZero();
      CalcPoseJacobian(model, Q, i, zero, Jac_temp, false);
      J.block(row_index, 0, 6, model.dof_count) = Jac_temp;
      row_index += 6;
    }

    row_index;
    for(unsigned int i=1; i<model.mBodies.size() - 1; ++i){
      I.block(row_index, 0, 6, 6) = model.X_base[i].toMatrix().transpose()*
          model.mBodies[i].mSpatialInertia*
          model.X_base[i].toMatrix();
      row_index += 6;
    }

    h = I*J*QDot;
    A = I*J;
  }
*/

  /*
  void massMatrixFromCCRBI(Model &model,
                           const VectorNd &Q,
                           const VectorNd &QDot,
                           Math::MatrixNd &H){

    Eigen::MatrixXd J(6*(model.mBodies.size() - 1), model.dof_count);
    unsigned int row_index = 0;
    Eigen::MatrixXd Jac_temp(6, model.dof_count);
    for(unsigned int i=1; i<model.mBodies.size() - 1; ++i){
      Eigen::Vector3d zero;
      zero.setZero();
      CalcPoseJacobian(model, Q, i, zero, Jac_temp, false);
      J.block(row_index, 0, 6, model.dof_count) = Jac_temp;
      row_index += 6;
    }

    SpatialMatrix I;
    CCM_CCRBI(model, Q, QDot, I, false);

     H = J.transpose()*I*J;
  }
  */

}
