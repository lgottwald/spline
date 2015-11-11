#ifndef _SPLINE_DERIVATIVE_HPP_
#define _SPLINE_DERIVATIVE_HPP_

#include <type_traits>
#include "DegreeTraits.hpp"
#include "FloatingPointTraits.hpp"
#include "SplineIntervalTraits.hpp"

namespace spline
{

template < typename Implementation >
class UnivariateSplineCurve;

template < typename Implementation >
class Polynomial;

/**
 * Class representing the derivative of a spline or polynomial function.
 */
template<typename Impl, int NDERIV, bool spline = std::is_base_of<UnivariateSplineCurve<Impl>, Impl>::value, bool poly = std::is_base_of<Polynomial<Impl>, Impl>::value >
class Derivative;

/**
 * Specialization of DegreeTraits to get the degree of the derivative.
 */
template<typename Impl, int NDERIV, bool spline, bool poly>
struct DegreeTraits<Derivative<Impl, NDERIV, spline, poly>>
{
   constexpr static int value = get_degree<Impl>() - NDERIV;
};

/**
 * Specialization of FloatingPointTraits to get the floating point type
 * of the derivative.
 */
template<typename Impl, int NDERIV, bool spline, bool poly>
struct FloatingPointTraits<Derivative<Impl, NDERIV, spline, poly>>
{
   using type = floating_point<Impl>;
};

/**
 * Specialization of IntervalTypeTraits to get the interval type
 * for derivatives of spline curves.
 */
template<typename Impl, int NDERIV>
struct IntervalTypeTraits<Derivative<Impl, NDERIV, true, false>>
{
   using type = interval_type<Impl>;
};

/**
 * Specialization for Derivative of UnivariateSplineCurve.
 */
template<typename Impl, int NDERIV>
class Derivative<Impl, NDERIV, true, false> : public UnivariateSplineCurve<Derivative<Impl, NDERIV, true, false>>
{
public:
   using real_t = floating_point< Derivative<Impl, NDERIV, true, false> >;
   using interval_t = interval_type< Derivative<Impl, NDERIV, true, false> >;

   // Just forward calls
   real_t getInfimumImpl( const interval_t interval ) const
   {
      return f.getInfimum( interval );
   }

   real_t getSupremumImpl( const interval_t interval ) const
   {
      return f.getSupremum( interval );
   }

   real_t getInfimumImpl() const
   {
      return f.getInfimum();
   }

   real_t getSupremumImpl() const
   {
      return f.getSupremum();
   }

   std::size_t numIntervalsImpl() const
   {
      return f.numIntervals();
   }

   interval_t findIntervalImpl( const real_t x, const interval_t interval ) const
   {
      return f.findInterval( x, interval );
   }

   template < int D, typename T >
   T evaluate( const T &x, const interval_t interval ) const
   {
      return f.template evaluate < D + NDERIV > ( x, interval );
   }

   Derivative( const UnivariateSplineCurve<Impl> &f ) : f( f )
   {
   }

private:
   const UnivariateSplineCurve<Impl> &f;
};

/**
 * Specialization for Derivative of polynomial.
 */
template<typename Impl, int NDERIV>
class Derivative<Impl, NDERIV, false, true> : public Polynomial<Derivative<Impl, NDERIV, false, true>>
{
public:
   template<typename T>
   T operator()( const T &x ) const
   {
      return poly.template derivative<NDERIV>( x );
   }

   template<int D = 1, typename T>
   T derivative( const T &x ) const
   {
      return poly.template derivative < D + NDERIV > ( x );
   }

   Derivative( const Polynomial<Impl> &poly ) : poly( poly ) {}

private:
   const Polynomial<Impl> &poly;
};

/**
 * Function to conveniently create the derivative
 * of a polynomial.
 */
template<int NDERIV = 1, typename Impl>
Derivative<Impl, NDERIV> differentiate(
   const Polynomial<Impl> &func )
{
   return Derivative<Impl, NDERIV>( func );
}

/**
 * Function to conveniently create the derivative
 * of an UnivariateSplineCurve.
 */
template<int NDERIV = 1, typename Impl>
Derivative<Impl, NDERIV> differentiate(
   const UnivariateSplineCurve<Impl> &func )
{
   return Derivative<Impl, NDERIV>( func );
}



}

#endif
