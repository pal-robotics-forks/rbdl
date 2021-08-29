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

  template <typename C>
  FixedBodyData<C> cast() const{

    FixedBodyData<C> casted;
    casted.mBaseTransform = mBaseTransform.template cast<C>();

    return casted;
  }
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

//  std::vector<FixedBodyData<T> > mFixedBodiesData;

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
  std::vector<Math::SpatialVector<T> > hdotc;

  template <typename C>
  ModelData<C> cast() const{

    ModelData<C> casted;

    for(size_t i=0; i<v.size(); ++i){
      casted.v.push_back(v[i]. template cast<C>());
    }

    for(size_t i=0; i<a.size(); ++i){
      casted.a.push_back(a[i]. template cast<C>());
    }

    for(size_t i=0; i<a_bias.size(); ++i){
      casted.a_bias.push_back(a_bias[i]. template cast<C>());
    }

    for(size_t i=0; i<acumulated_mass.size(); ++i){
      casted.acumulated_mass.push_back(acumulated_mass[i]. template cast<C>());
    }

    for(size_t i=0; i<X_base.size(); ++i){
      casted.X_base.push_back(X_base[i]. template cast<C>());
    }

    for(size_t i=0; i<X_J.size(); ++i){
      casted.X_J.push_back(X_J[i]. template cast<C>());
    }

    for(size_t i=0 ; i<v_J.size(); ++i){
      casted.v_J.push_back(v_J[i]. template cast<C>());
    }

    for(size_t i=0; i<c_J.size(); ++i){
      casted.c_J.push_back(c_J[i]. template cast<C>());
    }

    for(size_t i=0; i<S.size(); ++i){
      casted.S.push_back(S[i]. template cast<C>());
    }

    for(size_t i=0; i<multdof3_S.size(); ++i){
      casted.multdof3_S.push_back(multdof3_S[i]. template cast<C>());
    }

    // Bodies
    for(size_t i=0; i<X_lambda.size(); ++i){
      casted.X_lambda.push_back(X_lambda[i]. template cast<C>());
    }

    for(size_t i=0; i<c.size(); ++i){
      casted.c.push_back(c[i]. template cast<C>());
    }

//    /// @todo shoul this be in model data??
//    for(size_t i=0; i<mFixedBodiesData.size(); ++i){
//      casted.mFixedBodiesData.push_back(mFixedBodiesData[i]. template cast<C>());
//    }

    for(size_t i=0; i<f.size(); ++i){
      casted.f.push_back(f[i]. template cast<C>());
    }

    for(size_t i=0; i<IA.size(); ++i){
      casted.IA.push_back(IA[i]. template cast<C>());
    }

    for(size_t i=0; i<pA.size(); ++i){
      casted.pA.push_back(pA[i]. template cast<C>());
    }

    for(size_t i=0; i<U.size(); ++i){
      casted.U.push_back(U[i]. template cast<C>());
    }

    for(size_t i=0; i<d.size(); ++i){
      casted.d = d. template cast<C>();
    }

    for(size_t i=0; i<u.size(); ++i){
      casted.u = u. template cast<C>();
    }

    for(size_t i=0; i<I.size(); ++i){
      casted.I.push_back(I[i]. template cast<C>());
    }

    for(size_t i=0; i<multdof3_U.size(); ++i){
      casted.multdof3_U.push_back(multdof3_U[i]. template cast<C>());
    }

    for(size_t i=0; i<multdof3_Dinv.size(); ++i){
      casted.multdof3_Dinv.push_back(multdof3_Dinv[i]. template cast<C>());
    }

    for(size_t i=0; i<multdof3_u.size(); ++i){
      casted.multdof3_u.push_back(multdof3_u[i]. template cast<C>());
    }

    /// @todo this should go in the fix model data?
    for(size_t i=0; i<multdof3_w_index.size(); ++i){
      casted.multdof3_w_index = multdof3_w_index;
    }

    for(size_t i=0; i<Ic.size(); ++i){
      casted.Ic.push_back(Ic[i]. template cast<C>());
    }

    for(size_t i=0; i<hc.size(); ++i){
      casted.hc.push_back(hc[i]. template cast<C>());
    }

    return casted;
  }

};

using ModelDatad = ModelData<double>;

}

#endif
