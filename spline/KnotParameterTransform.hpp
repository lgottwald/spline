#ifndef _SPLINE_KNOT_PARAMETER_TRANSFORM_HPP_
#define _SPLINE_KNOT_PARAMETER_TRANSFORM_HPP_

#include <cmath>
#include "BSplineCurve.hpp"

namespace spline
{


/** Pass floating point values by value but
 *  other types by reference
 */
template<typename T>
using real_arg_t = typename std::conditional<std::is_floating_point<T>::value, T, T &>::type;

/**
 * Inverse the logarithmic transformation of parameters
 * to get the knots. If all parameters are zero the knots will
 * be equally spaced. The first knot is always a
 * and the last knot is b.
 */
template<class REAL, int DEGREE>
REAL inverse_parameter_log_transform( const REAL *params, real_arg_t<const REAL> a, real_arg_t<const REAL> b, BSplineCurve<DEGREE, REAL> &curve )
{
   //use sopace for coefficients as workspace
   auto &workspace = curve.getCoefficients();
   std::size_t nparams = workspace.size() - DEGREE - 1;

   //compute exp of params and store intermediate results in knot vector
   for( std::size_t i = 0; i < nparams; ++i )
   {
      curve.getKnot( i ) = exp( params[i] );
   }

   REAL Z = curve.getKnot( nparams - 1 );

   for( int i = nparams - 2; i >= 0 ; --i )
   {
      Z = curve.getKnot( i ) + curve.getKnot( i ) * Z;
   }

   workspace[0] = ( b - a ) / ( Z + 1 );

   for( std::size_t i = 1; i < nparams + 1; ++i )
   {
      workspace[i] = workspace[i - 1] * curve.getKnot( i - 1 );
   }

   //now compute the knots
   curve.getKnot( 0 ) = a;
   REAL dist = (b-a);
   for( std::size_t i = 1; i < nparams + 2; ++i )
   {
      curve.getKnot( i ) = curve.getKnot( i - 1 ) + workspace[i - 1];
      REAL d =  curve.getKnot( i ) - curve.getKnot( i - 1 );
      dist = dist > d ? d : dist;
   }
   return dist;
}

/**
 * Transform knot values into corresponding parameters.
 */
template<class REAL, int DEGREE, class VEC>
void parameter_log_transform( const BSplineCurve<DEGREE, REAL> &curve, VEC &params )
{
   std::size_t n = curve.getNumKnots() - 2;
   params.resize( n + 1 );

   for( std::size_t i = 0; i < params.size(); ++i )
   {
      params[i] = curve.getKnot( i + 1 ) - curve.getKnot( i );
   }

   for( std::size_t i = 0; i < n; ++i )
   {
      params[i] = log( params[i + 1] / params[i] );
   }

   params.resize( n );
}

}//spline

#endif