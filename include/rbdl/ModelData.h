#ifndef _MODEL_DATA_
#define _MODEL_DATA_

#include <rbdl/rbdl_math.h>

/** \brief Namespace for all structures of the RigidBodyDynamics library
*/
namespace RigidBodyDynamics {

template <class T>
struct FixedBodyData{
  /// \brief Transforms spatial quantities expressed for the parent to the
  // fixed body.
  Math::SpatialTransform<T> mBaseTransform;
};

using FixedBodyDataD = FixedBodyData<double>;

template <typename T>
class ModelData{

public:

  // State information
  /// \brief The spatial velocity of the bodies
  std::vector<Math::SpatialVector<T> > v;
  /// \brief The spatial acceleration of the bodies
  std::vector<Math::SpatialVector<T> > a;
  /// \brief The spatial bias acceleration of the bodies
  std::vector<Math::SpatialVector<T> > a_bias;

  /// \brief For computing COM jacobian efficiently
  std::vector<Math::SpatialMatrix<T> > acumulated_mass;

  /// \brief Transformation from the base to bodies reference frame
  std::vector<Math::SpatialTransform<T> > X_base;

  // Joint state variables
  std::vector<Math::SpatialTransform<T> > X_J;
  std::vector<Math::SpatialVector<T> > v_J;
  std::vector<Math::SpatialVector<T> > c_J;

  /// \brief The joint axis for joint i
  std::vector<Math::SpatialVector<T> > S;

  // Special variables for joints with 3 degrees of freedom
  /// \brief Motion subspace for joints with 3 degrees of freedom
  std::vector<Math::Matrix63<T> > multdof3_S;

  // Bodies

  /** \brief Transformation from the parent body to the current body
       * \f[
       *	X_{\lambda(i)} = {}^{i} X_{\lambda(i)}
       * \f]
       */
  std::vector<Math::SpatialTransform<T> > X_lambda;

  // Dynamics variables
  /// \brief The velocity dependent spatial acceleration
  std::vector<Math::SpatialVector<T> > c;


  std::vector<FixedBodyData<T> > mFixedBodiesData;

  /// \brief Internal forces on the body (used only InverseDynamics())
  std::vector<Math::SpatialVector<T> > f;


  /// \brief The spatial inertia of the bodies
  std::vector<Math::SpatialMatrix<T> > IA;
  /// \brief The spatial bias force
  std::vector<Math::SpatialVector<T> > pA;
  /// \brief Temporary variable U_i (RBDA p. 130)
  std::vector<Math::SpatialVector<T> > U;
  /// \brief Temporary variable D_i (RBDA p. 130)
  Math::VectorN<T> d;
  /// \brief Temporary variable u (RBDA p. 130)
  Math::VectorN<T> u;

  std::vector<Math::SpatialRigidBodyInertia<T > > I;

  // Special variables for joints with 3 degrees of freedom
  /// \brief Motion subspace for joints with 3 degrees of freedom
  std::vector<Math::Matrix63<T> > multdof3_U;
  std::vector<Math::Matrix3<T> > multdof3_Dinv;
  std::vector<Math::Vector3<T> > multdof3_u;
  std::vector<unsigned int> multdof3_w_index;

  /// \brief The spatial inertia of body i (used only in
  ///  CompositeRigidBodyAlgorithm())
  std::vector<Math::SpatialRigidBodyInertia<T> > Ic;
  std::vector<Math::SpatialVector<T> > hc;

};

using ModelDatad = ModelData<double>;

}

#endif
