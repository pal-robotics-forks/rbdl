#ifndef _RBDL_UTILS_
#define _RBDL_UTILS_

#include <Eigen/Dense>
#include <exception>

inline void createGeneralizedVector(const Eigen::Vector3d &floatingBasePosition, const Eigen::Quaterniond &floatingBaseOrientation,
                                    const Eigen::VectorXd &jointStates, Eigen::VectorXd &state){

  if(state.rows() != jointStates.rows() + 7){
    throw std::runtime_error("Missmatch between vectors creating generalized vector");
  }

  state.segment(0, 3) = floatingBasePosition;
  state(3) = floatingBaseOrientation.x();
  state(4) = floatingBaseOrientation.y();
  state(5) = floatingBaseOrientation.z();
  state(state.rows() - 1) = floatingBaseOrientation.w();

  state.segment(6, jointStates.rows()) = jointStates;
}

inline void createGeneralizedVelocityVector(const Eigen::Vector3d &floatingBaseLinearVelocity, const Eigen::Vector3d &floatingBaseAngularVelocity,
                                    const Eigen::VectorXd &jointStatesVelocity, Eigen::VectorXd &state){

  if(state.rows() != jointStatesVelocity.rows() + 6){
    throw std::runtime_error("Missmatch between vectors creating generalized vector");
  }

  state.segment(0, 3) = floatingBaseLinearVelocity;
  state.segment(3, 3) = floatingBaseAngularVelocity;
  state.segment(6, jointStatesVelocity.rows()) = jointStatesVelocity;
}

inline void createGeneralizedAccelerationVector(const Eigen::Vector3d &floatingBaseLinearAcceleration, const Eigen::Vector3d &floatingBaseAngularAcceleration,
                                    const Eigen::VectorXd &jointStatesAcceleration, Eigen::VectorXd &state){

  if(state.rows() != jointStatesAcceleration.rows() + 6){
    throw std::runtime_error("Missmatch between vectors creating generalized vector");
  }

  state.segment(0, 3) = floatingBaseLinearAcceleration;
  state.segment(3, 3) = floatingBaseAngularAcceleration;
  state.segment(6, jointStatesAcceleration.rows()) = jointStatesAcceleration;
}




#endif
