#include <rbdl/Energy.h>
#include <rbdl/Kinematics.h>

namespace RigidBodyDynamics {

  using namespace Math;

  double CalcTotalMass(Model &model){
    double totalMass = 0.;
    for(size_t i=0; i<model.mBodies.size(); ++i){
      totalMass += model.mBodies[i].mMass;
    }
    return totalMass;
  }

  // We have to acumulate the spatial transforms, it seems that the cross product of acumulationg the
  // com displacements is a wrong assumption
  void CalcAcumulatedMass(Model &model,
                          ModelDatad &model_data,
                          const VectorNd &Q){
    assert(model_data.acumulated_mass.size() == model.mBodies.size());

    for(unsigned int i = 1; i<model.mBodies.size(); ++i){
      Vector3d comi = CalcBodyToBaseCoordinates(model, model_data, Q, i, model.mBodies[i].mCenterOfMass, false);
      model_data.acumulated_mass[i] = model.mBodies[i].mMass*Xtrans_mat(comi);
    }

    for (unsigned int i = model.mBodies.size() - 1; i > 0; i--) {
      unsigned int lambda = model.lambda[i];
      model_data.acumulated_mass[lambda] += model_data.acumulated_mass[i];
    }
  }


  Math::Vector3d CalcCOM(Model &model,
                         ModelDatad &model_data,
                         const Math::VectorNd &Q,
                         bool update_kinematics){
    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematicsCustom<double>(model, model_data, &Q, NULL, NULL);
    }
    Eigen::Vector3d com;
    com.setZero();
    double total_mass = 0;

    for(unsigned int i=1; i<model.mBodies.size(); ++i){
      int body_id = i;//rbdl_model_.GetBodyId(link_names_[i].c_str());
      Eigen::Vector3d link_com;
      link_com = CalcBodyToBaseCoordinates(model, model_data, Q,
                                           body_id, model.mBodies[body_id].mCenterOfMass, false);

      com += model.mBodies[body_id].mMass*link_com;
      total_mass +=  model.mBodies[body_id].mMass;
    }
    return  com/total_mass;

  }


  Math::Vector3d CalcCOMVelocity(Model &model,
                                 ModelDatad &model_data,
                                 const Math::VectorNd &Q,
                                 const Math::VectorNd &QDot,
                                 bool update_kinematics){

    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematicsCustom<double>(model, model_data, &Q, &QDot, NULL);
    }
    Eigen::Vector3d com_vel;
    com_vel.setZero();
    double total_mass = 0;

    for(unsigned int i=1; i<model.mBodies.size(); ++i){
      int body_id = i;//rbdl_model_.GetBodyId(link_names_[i].c_str());
      Eigen::Vector3d link_com_vel;
      link_com_vel = CalcPointVelocity(model, model_data, Q, QDot,
                                       body_id, model.mBodies[body_id].mCenterOfMass, false);

      com_vel += model.mBodies[body_id].mMass*link_com_vel;
      total_mass +=  model.mBodies[body_id].mMass;
    }
    return  com_vel/total_mass;

  }

  Math::Vector3d CalcCOMAcceleartion(Model &model,
                                     ModelDatad &model_data,
                                     const Math::VectorNd &Q,
                                     const Math::VectorNd &QDot,
                                     const Math::VectorNd &QDDot,
                                     bool update_kinematics){
    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematics(model, model_data, Q, QDot, QDDot);
    }
    Eigen::Vector3d com_acc;
    com_acc.setZero();
    double total_mass = 0;

    for(unsigned int i=1; i<model.mBodies.size(); ++i){
      int body_id = i;
      Eigen::Vector3d link_com_acc;
      link_com_acc = CalcPointAcceleration(model, model_data, Q, QDot, QDDot,
                                           body_id, model.mBodies[body_id].mCenterOfMass, false);

      com_acc += model.mBodies[body_id].mMass*link_com_acc;
      total_mass +=  model.mBodies[body_id].mMass;
    }
    return  com_acc/total_mass;
  }

  Math::Vector3d CalcCOMAccelerationBias(Model &model,
                                         ModelDatad &model_data,
                                         const Math::VectorNd &Q,
                                         const Math::VectorNd &QDot,
                                         bool update_kinematics){
    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematicsCustom<double>(model, model_data, &Q, &QDot, NULL);
    }
    Eigen::Vector3d com_acc_bias;
    com_acc_bias.setZero();
    double total_mass = 0;

    for(unsigned int i=1; i<model.mBodies.size(); ++i){
      int body_id = i;
      Eigen::Vector3d link_com_acc_bias;
      link_com_acc_bias = CalcPointAccelerationBias(model, model_data, Q, QDot,
                                                    body_id, model.mBodies[body_id].mCenterOfMass, false);

      com_acc_bias += model.mBodies[body_id].mMass*link_com_acc_bias;
      total_mass +=  model.mBodies[body_id].mMass;
    }
    return  com_acc_bias/total_mass;
  }



  void CalcCOMJacobian (
      Model &model,
      ModelDatad &model_data,
      const Math::VectorNd &Q,
      Math::MatrixNd &COMJ,
      bool update_kinematics){

    LOG << "-------- " << __func__ << " --------" << std::endl;

    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematicsCustom<double>(model, model_data, &Q, NULL, NULL);
    }

    double total_mass = 0;
    for (unsigned int j = 1; j < model.mBodies.size(); j++) {
      double mass_link = model.mBodies[j].mMass;
      total_mass += mass_link;
    }


    CalcAcumulatedMass(model, model_data, Q);


    for (unsigned int j = 1; j < model.mBodies.size(); j++) {
      unsigned int q_index = model.mJoints[j].q_index;

      //if(model.mJoints[j].mJointType != JointTypeCustom){
        if (model.mJoints[j].mDoFCount == 1) {
          COMJ.block(0,q_index, 3, 1) =
              (model_data.acumulated_mass[j]
               *(model_data.X_base[j].inverse().toMatrix()
                 * model_data.S[j])).block(3,0,3,1);
        } else if (model.mJoints[j].mDoFCount == 3) {
          COMJ.block(0, q_index, 3, 3) =
              (model_data.acumulated_mass[j]
                 * (model_data.X_base[j].inverse()).toMatrix()
               * model_data.multdof3_S[j]).block(3,0,3,3);
        }
//      }
//      else{
//        unsigned int k = model.mJoints[j].custom_joint_index;

//        COMJ.block(0, q_index, 3, model.mCustomJoints[k]->mDoFCount) =
//            (model_data.acumulated_mass[j]
//               * (model_data.X_base[j].inverse()).toMatrix()
//             * model.mCustomJoints[k]->S).block(
//              3,0,3,model.mCustomJoints[k]->mDoFCount);
//      }

    }

    COMJ = COMJ/total_mass;

  }

  void CalcCOMJacobian_inefficient (
      Model &model,
      ModelDatad &model_data,
      const VectorNd &Q,
      MatrixNd &COMJ,
      bool update_kinematics
      ) {

    assert (COMJ.rows() == 3 && COMJ.cols() == model.dof_count );

    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematicsCustom<double>(model, model_data, &Q, NULL, NULL);
    }

    double total_mass = 0;
    for (unsigned int j = 1; j < model.mBodies.size(); j++) {
      double mass_link = model.mBodies[j].mMass;
      total_mass += mass_link;
    }

    COMJ.setZero();
    Vector3d zero;
    zero.setZero();
    Eigen::MatrixXd jacobian_link(3, model.dof_count);
    for (size_t j = 1; j < model.mBodies.size(); ++j) {
      jacobian_link.setZero();
      double mass_link = model.mBodies[j].mMass;
      /// @todo: make optinional to not especify extra displacement of the tip
      CalcPointJacobian(model, model_data, Q, j,
                        model.mBodies[j].mCenterOfMass, jacobian_link, false);
      COMJ += mass_link*jacobian_link;
    }
    COMJ = COMJ/total_mass;
  }

  //  void CCM_CCRBI_Jacobian_com(Model &model,
  //                              const Math::VectorNd &Q,
  //                              Math::MatrixNd &AG,
  //                              bool update_kinematics){

  //  }

  /// @todo: Document
  SpatialVectord CalcEnergy_inefficient(
      Model &model,
      ModelDatad &model_data,
      const Math::VectorNd &Q,
      const Math::VectorNd &QDot,
      bool update_kinematics,
      unsigned int method
      ){

    if (update_kinematics) {
      //Compute energy & COM & mass of the robot;
      RigidBodyDynamics::UpdateKinematicsCustom<double>(model, model_data, &Q, &QDot, NULL);
    }

    SpatialVectord spatial_moment = RigidBodyDynamics::Math::SpatialVectord::Zero();

    /// @todo: Is there a problem with the floating base when calculation COM, is the base_link COM getting stuck at the origin
    /// of the inertial frame?

    /// @warning We start with body 1 becouse the first link is the root world link which seems to have 1 kg mass
    for(size_t i=1; i<model.mBodies.size(); ++i){
      unsigned int body_id = i;//rbdl_model_.GetBodyId(link_names_[i].c_str());

      if(method == 0){
        /// X'*I*X*(X^-1*v)
        spatial_moment += model_data.X_base[body_id].toMatrix().transpose()*
            model.I[body_id].toMatrix()*
            model_data.X_base[body_id].toMatrix()*
            model_data.X_base[body_id].inverse().toMatrix()*model_data.v[body_id];
      }
      else if(method == 1){
        //Optimization using spatial operators, in the book there is a further optimization to reduce operations
        //    spatial_moment += model_data.X_base[body_id].apply(model.mBodies[body_id].mSpatialInertia*
        //        (model_data.X_base[body_id].applyInverse(model_data.v[body_id])));

        spatial_moment += model_data.X_base[body_id].toMatrix().transpose()*
            model.I[body_id].toMatrix()*
            model_data.X_base[body_id].toMatrix()*
            (model_data.X_base[body_id].inverse().apply(model_data.v[body_id]));

      }
      else{
        //Spatial momentum vectors have the same properties as spatial force vectors
        // Optimization using spatial operators, in the book there is a further optimization to reduce operations
        spatial_moment += model_data.X_base[body_id].applyTranspose(model.I[body_id].toMatrix()*(model_data.v[body_id]));
      }

      //std::cerr<<"Spatial moment: "<<i<<"  is: "<<spatial_moment.transpose()<<std::endl;
    }

    return spatial_moment;

  }


}
