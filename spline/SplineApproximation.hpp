#ifndef _SPLINE_APPROXIMATION_HPP_
#define _SPLINE_APPROXIMATION_HPP_

#include <vector>
#include <cmath>
#include <cassert>
#include <omp.h>
#include <simd/alloc.hpp>
#include <simd/pack.hpp>
#include "BSplineCurve.hpp"
#include "PiecewisePolynomial.hpp"
#include <cpplsq/gn_sbfgs_min.hpp>
#include "KnotParameterTransform.hpp"

namespace spline
{

/**
 * Implementation of the parameter log transformation
 * for use with the least squares minimizer.
 *
 */
template<typename REAL, typename Func>
struct ParameterLogTransformator
{
   ParameterLogTransformator( const Func &f, REAL a, REAL b, REAL minkntdist ) : f( f ), a( a ), b( b ), minkntdist(minkntdist) {}

   /**
    * Compute the knots from the transformed parameters
    * and compute the variation diminishing approximation
    * coefficients using the stored function.
    * Returns a pointer to the BSplineCurve representing the result.
    */
   const BSplineCurve< 3, cpplsq::SingleDiff<REAL> > *operator()( const cpplsq::SingleDiff<REAL> *params )
   {
      auto kntdist = inverse_parameter_log_transform( params, a, b, sdcurve );
      if(kntdist.getValue() < minkntdist )
         return nullptr;
      sdcurve.approximateVariationDiminishing( f.get() );
      return &sdcurve;
   }

   /**
    * Compute the knots from the transformed parameters
    * and compute the variation diminishing approximation
    * coefficients using the stored function.
    * Returns a pointer to the BSplineCurve representing the result.
    */
   const BSplineCurve< 3, cpplsq::MultiDiff<REAL> > *operator()( const cpplsq::MultiDiff<REAL> *params )
   {
      auto kntdist = inverse_parameter_log_transform( params, a, b, mdcurve );
      if(kntdist.getValue() < minkntdist )
         return nullptr;
      mdcurve.approximateVariationDiminishing( f.get() );
      return &mdcurve;
   }

   /**
    * Allocate the space for N paramters. Allocation should not
    * occur before since the MultiDiff context will only be initialized
    * when this function is called.
    */
   void num_parameters( std::size_t N )
   {
      mdcurve.setKnotInterval( a, b, N + 2 );
      sdcurve.setKnotInterval( a, b, N + 2 );
   }

private:
   std::reference_wrapper<const Func> f;
   REAL a, b;
   REAL minkntdist;
   BSplineCurve< 3, cpplsq::MultiDiff<REAL> > mdcurve;
   BSplineCurve< 3, cpplsq::SingleDiff<REAL> > sdcurve;
};

/**
 * A residual for approximating a spline curve.
 * Works on with the BSplineCurve
 * returned by the ParameterLogTransformator.
 */
template<typename REAL>
struct SplineResidual
{
   /**
    * Create residual evaluating curve at x and expecting y.
    * The error is then curve(x)-y but is normalized by
    * the given value y_norm.
    */
   SplineResidual( REAL x, REAL y, REAL y_norm ) : x( x ), y( y ), y_norm( y_norm ) {}

   /**
    * Computes the error for the given curve.
    */
   template<typename T>
   T operator()( const BSplineCurve<3, T> *curvePtr )
   {
      if(!curvePtr)
         return std::numeric_limits<REAL>::quiet_NaN();
      auto &curve = *curvePtr;
      return ( y - curve( T( x ) ) ) / y_norm;
   }

private:
   REAL x;
   REAL y;
   REAL y_norm;
};


/**
 * Approximate the given function by a
 * cubic BSplineCurve.
 * \param f       The function to approximate
 * \param x       vector of x values for creating the residuals
 *                and the interval in which the function is approximated.
 * \param nknots  Number of knots used for the spline.
 * \param epsilon Value of epsilon used for the approximation.
 *
 */
template <typename REAL, typename Function>
BSplineCurve<3, REAL> ApproximateFunction(
   const Function &f,
   const std::vector<REAL> &x,
   const std::vector<REAL> &y,
   int nknots,
   REAL delta = 1e-2,
   REAL epsilon = 1e-6,
   REAL minkntdist = 1e-6
)
{
   const std::size_t nresiduals = x.size();
   const REAL a = x.front();
   const REAL b = x.back();

   std::vector<SplineResidual<REAL>> r;
   r.reserve( nresiduals );

   for( std::size_t i = 0; i < nresiduals; ++i )
   {
      r.emplace_back( x[i], y[i], y[i] + delta );
   }

   simd::aligned_vector<REAL> params( nknots - 2, 0.0 );
   cpplsq::gn_sbfgs_min<cpplsq::Silent>( epsilon, params, r, ParameterLogTransformator<REAL, Function>( f, a, b, minkntdist ) );

   BSplineCurve<3, REAL> curve( a, b, nknots );
   inverse_parameter_log_transform( params.data(), a, b, curve );
   curve.approximateVariationDiminishing( f );
   return curve;
}


/**
 * Approximate piecewise linear function by cubic spline. Try
 * to use as few knots as possible to achieve the desired error
 * rates.
 *
 * \param f                   piecewise linear function
 * \param a                   left endpoint of interval to approximate the function
 * \param b                   right endpoint of interval to approximate the function
 * \param max_rel_err         maximum relative error to nonzero values of piecewise linear function
 * \param delta               delta to use for computing the relative error.
 */
template <typename REAL, typename IMPL>
BSplineCurve<3, REAL> ApproximatePiecewiseLinear(
  const UnivariateSplineCurve<IMPL> &f,
  REAL a,
  REAL b,
  REAL &max_rel_err = 0.05,
  REAL delta = 1e-2,
  REAL epsilon = 1e-6,
  REAL minkntdist = 1e-6
)
{
   static_assert(get_degree<IMPL>() == 1, "Given function should be a spline with degree one");
   using interval_t = interval_type< IMPL >;
   int nknts = 5;
   REAL maxerr;
   // For piecewise linear functions the cubic spline will have the largest error at
   // the breakpoints. So use only the breakpoints for measuring the error.
   std::vector<REAL> x( f.numIntervals() + 1 );
   x.front() = a;

   for( std::size_t i = 1; i < x.size() - 1; ++i )
   {
      x[i] = f.getInfimum( i );
   }

   x.back() = b;
   std::vector<REAL> y;
   y.reserve( x.size() );
   std::transform( x.begin(), x.end(), std::back_inserter( y ), [&](REAL x) { return f(x); } );
   auto y_minmax = std::minmax_element( y.begin(), y.end() );
   delta  = -*y_minmax.first + delta;

   while( 1 )
   {
      BSplineCurve<3, REAL> curve = ApproximateFunction( f, x, y, nknts, delta, epsilon, minkntdist );
      maxerr = 0;
      for( interval_t i = 1; i < interval_t(f.numIntervals()); ++i )
      {
         REAL x = f.getInfimum( i );
         REAL fx = f( x );
         REAL err = ( fx - curve( x ) ) / ( fx + delta );
         maxerr = std::max(maxerr, std::abs(err) );
      }
      if( maxerr > max_rel_err )
      {
         goto refit;
      }

      max_rel_err = maxerr;

      return curve;
   refit:
      ++nknts;

      if( nknts > 300 )
      {
         max_rel_err = maxerr;
	 return curve;
      }
   }
}

} // namespace spline

#endif

