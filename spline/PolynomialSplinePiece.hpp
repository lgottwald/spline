#ifndef _SPLINE_POLYNOMIAL_PIECE_HPP_
#define _SPLINE_POLYNOMIAL_PIECE_HPP_

#include "DegreeTraits.hpp"
#include "SplineIntervalTraits.hpp"
#include "FloatingPointTraits.hpp"

namespace spline
{

template<typename Impl>
class Polynomial;

template<typename Impl>
class UnivariateSplineCurve;

template<typename Impl>
class PolynomialPiece;

template<typename Impl>
struct FloatingPointTraits<PolynomialPiece<Impl>>
{
   using type = floating_point<Impl>;
};

template<typename Impl>
struct DegreeTraits<PolynomialPiece<Impl>>
{
   constexpr static int value = get_degree<Impl>();
};

/**
 * Implement the polynomial base class for
 * polynomial pieces of an UnivariateSplineCurve.
 */
template<typename Impl>
class PolynomialPiece : public Polynomial<PolynomialPiece<Impl>>
{
public:

   /**
    * Evaluate the polynomial piece at x value.
    */
   template < typename T >
   T operator()( const T &x ) const
   {
      return spline.template evaluate<0>( x, interval );
   }

   /**
    * Evaluate the polynomial pieces derivative at x value.
    */
   template < int D = 1, typename T >
   T derivative( const T &x ) const
   {
      return spline.template evaluate<D>( x, interval );
   }

   /**
    * Create a polynomial piece given an UnivariateSplineCurve and an interval.
    */
   PolynomialPiece( interval_type<Impl> interval, const UnivariateSplineCurve<Impl> &spline ) :
      interval( interval ), spline( spline )
   {
   }

private:
   interval_type<Impl> interval;
   const UnivariateSplineCurve<Impl> &spline;
};


}

#endif
