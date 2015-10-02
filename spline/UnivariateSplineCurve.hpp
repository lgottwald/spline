#ifndef _UNIVARIATE_SPLINE_CURVE_HPP_
#define _UNIVARIATE_SPLINE_CURVE_HPP_

#include "Polynomial.hpp"
#include "SplineIntervalTraits.hpp"
#include "DegreeTraits.hpp"
#include "FloatingPointTraits.hpp"
#include <simd/pack.hpp>
#include <type_traits>
#include <vector>
#include <algorithm>
#include "Derivative.hpp"
#include "PolynomialSplinePiece.hpp"
#include <cpplsq/AutoDiff.hpp>

namespace spline
{

/**
 * Crtp base class for univariate spline curves.
 */
template < typename Implementation >
class UnivariateSplineCurve
{
   static_assert( has_degree<Implementation>(), "DegreeTraits required for implementation of UnivariateSplineCurve" );
   static_assert( floating_point_type_defined<Implementation>(), "FloatingPointTraits required for implementation of UnivariateSplineCurve" );
   static_assert( interval_type_defined<Implementation>(), "IntervalTypeTraits required for implementation of UnivariateSplineCurve" );
public:
   UnivariateSplineCurve() : interval_( 0 )
   {
   }

   //=========== Interface that needs to be impemented in derived class ==============
   using interval_t = interval_type<Implementation>;
   using real_t = floating_point<Implementation>;
   constexpr static int DEGREE = get_degree<Implementation>();

   /**
    * Get the infimum of given interval.
    */
   real_t getInfimum( const interval_t interval ) const
   {
      return getImpl().getInfimumImpl( interval );
   }

   /**
    * Get the supremum of given interval.
    */
   real_t getSupremum( const interval_t interval ) const
   {
      return getImpl().getSupremumImpl( interval );
   }

   /**
    * Get the number of intervals.
    */
   std::size_t numIntervals() const
   {
      return getImpl().numIntervalsImpl();
   }

   /**
    * Find the interval for a given x value starting the search
    * with the given hint.
    */
   template < typename T >
   interval_t findInterval( const T &x, const interval_t hint ) const
   {
      return getImpl().template findIntervalDispatch<T>( x, hint );
   }

   /**
    * Evaluate the spline curve at value x in the given interval.
    *
    * \param x           x value for evaluation
    * \param interval    interval to evaluate the curve in
    */
   template < int D = 0, typename T >
   cpplsq::ValueType<T> evaluate( const T &x, const interval_t interval ) const
   {
      return getImpl().template evaluate<D>( x, interval );
   }

   //============ Dispatch call of find interval for simd::pack and adouble types =============

   /**
    * Statically dispatch the call to find the interval.
    * For simd::pack type the first value of the pack
    * is used to find the interval.
    *
    * \param x    The value to find the interval for.
    * \param hint hint for interval to start search.
    */
   template < typename T >
   interval_t findIntervalDispatch(
      const typename std::enable_if<std::is_same< T, simd::pack<real_t> >::value, T >::type &x,
      const interval_t hint
   ) const
   {
      return getImpl().findIntervalImpl( x[0], hint );
   }

   /**
   * Statically dispatch the call to find the interval.
   * If the type is the same as the curves type for reals
   * then just forward.
   *
   * \param x    The value to find the interval for.
   * \param hint hint for interval to start search.
   */
   template < typename T,
            typename std::enable_if<std::is_same<T, real_t>::value, int >::type = 0 >
   interval_t findIntervalDispatch(
      const T &x,
      const interval_t hint ) const
   {
      return getImpl().findIntervalImpl( x, hint );
   }

   /**
    * Statically dispatch the call to find the interval.
    * For AutoDiff types used their value.
    *
    * \param x    The value to find the interval for.
    * \param hint hint for interval to start search.
    */
   template < typename T,
            typename std::enable_if < cpplsq::is_diff_type<T>() && !std::is_same<T, real_t>::value, int >::type = 0 >
   interval_t findIntervalDispatch(
      const T &x,
      const interval_t hint ) const
   {
      return getImpl().findIntervalImpl( x.getValue(), hint );
   }

   /**
    * Find interval and store result to guide the next interval search.
    *
    * \param x    The value to find the interval for.
    * \param hint hint for interval to start search.
    */
   template < typename T >
   interval_t findInterval( const T &x ) const
   {
      auto i = findInterval( x, interval_ );
      interval_ = i;
      return i;
   }

   //============ evaluation is forwarded to implementation and interval is cached ==============

   /**
    * Evaluate the univariate spline curve at the given x value.
    * The interval of the x value is searched and cached to guide
    * future searches.
    * NOT guaranteed to be threadsafe because of this caching
    * though it should not fail on most platforms where
    * concurrent read and write will not yield
    * garbage values.
    *
    * \param x    The value to evaluate this curve at.
    * \return     The value of the curve at the given x value. The type
    *             will be the type of the given x or if x is an AutoDiff
    *             expression template the type will be the corresponding
    *             type to store the expressions value.
    */
   template < typename T >
   cpplsq::ValueType<T> operator()( const T &x ) const
   {
      auto i = findInterval( x, interval_ );
      interval_ = i;
      return getImpl().template evaluate<0>( x, i );
   }

   /**
   * Evaluate the univariate spline curves derivative at the given
   * x value. The interval of the x value is searched and cached
   * to guide future searches.
   * NOT guaranteed to be threadsafe because of this caching
   * though it should not fail on most platforms where
   * concurrent read and write will not yield
   * garbage values.
   *
   * \param x    The value to evaluate this curve at.
   * \tparam D   The number of the derivative to evaluate. E.g. D=1 evaluates the first derivative.
   * \return     The value of the curves derivative at the given x value.
   *             The type will be the type of the given x or if x is an AutoDiff
   *             expression template the type will be the corresponding
   *             type to store the expressions value.
   */
   template < int D = 1, typename T >
   cpplsq::ValueType<T> derivative( const T &x ) const
   {
      int i = findInterval( x, interval_ );
      interval_ = i;
      return getImpl().template evaluate<D>( x, i );
   }


   /**
    * Finds x values for the minimum and maximum values
    * within the given interval using the given tolerance.
    *
    * \param a          left endpoint of interval
    * \param b          right endpoint of interval
    * \param tolerance  The tolerance used to find the values.
    * \return A pair with the first value being the x value where the curve attains
    *         the minimum within the given interval and the second value
    *         being the x value where the curve attains its maximum within the interval.
    */
   std::pair<real_t, real_t> findMinMax( real_t a, real_t b, real_t tolerance ) const
   {
      std::pair<real_t, real_t> minmax;
      real_t minval, maxval;
      findMinMax( minval, maxval, minmax.first, minmax.second, a, b, tolerance );
      return minmax;
   }

   /**
   * Finds minimum and maximum values in the given interval using the given tolerance.
   *
   * \param a          left endpoint of interval
   * \param b          right endpoint of interval
   * \param tolerance  The tolerance used to find the min and max values.
   * \return A pair with the first value being the minimum value attained by the curve
   *         within the given interval and the second value being the maximum value.
   */
   std::pair<real_t, real_t> findMinMaxVal( real_t a, real_t b, real_t tolerance ) const
   {
      std::pair<real_t, real_t> minmaxval;
      real_t minarg, maxarg;
      findMinMax( minmaxval.first, minmaxval.second, minarg, maxarg, a, b, tolerance );
      return minmaxval;
   }

   /**
    * Finds x values for the minimum value within the
    * given interval using the given tolerance.
    *
    * \param a          left endpoint of interval
    * \param b          right endpoint of interval
    * \param tolerance  The tolerance used to find the values.
    * \return The x value where this curve attains the minimum within the given interval.
    */
   real_t findMin( real_t a, real_t b, real_t tolerance ) const
   {
      auto deriv = differentiate( *this );
      a = std::max(a, getInfimum(0));
      b = std::min(b, getSupremum(numIntervals()-1));
      interval_t l = findInterval( a, interval_ );
      interval_t r = findInterval( b, l );
      interval_ = r;
      real_t min = a;

      real_t minval = evaluate<0>( a, l );
      {
         real_t val = evaluate<0>( b, r );

         if( val < minval )
         {
            minval = val;
            min = b;
         }
      }
      ++r;

      for( interval_t i = l; i != r; ++i )
      {
         std::complex<real_t> localRoots[DEGREE - 1];
         auto pp = deriv.getPolynomialPiece( i );
         int nroots = pp.getRoots( localRoots, tolerance * 1e-2 );

         for( int rt = 0; rt < nroots; ++rt )
         {
            if( std::abs( localRoots[rt].imag() ) * 2 <= tolerance &&
                  localRoots[rt].real() >= std::max( a, getInfimum( i ) ) &&
                  localRoots[rt].real() <= std::min( b, getSupremum( i ) )
              )
            {
               real_t val = evaluate<0>( localRoots[rt].real(), i );

               if( val < minval )
               {
                  min = localRoots[rt].real();
                  minval = val;
               }

            }
         }
      }

      return min;
   }

   /**
    * Finds x values for the maximum value within the
    * given interval using the given tolerance.
    *
    * \param a          left endpoint of interval
    * \param b          right endpoint of interval
    * \param tolerance  The tolerance used to find the values.
    * \return The x value where this curve attains the maximum within the given interval.
    */
   real_t findMax( real_t a, real_t b, real_t tolerance ) const
   {
      auto deriv = differentiate( *this );
      a = std::max(a, getInfimum(0));
      b = std::min(b, getSupremum(numIntervals()-1));
      interval_t l = findInterval( a, interval_ );
      interval_t r = findInterval( b, l );
      interval_ = r;
      real_t max = a;

      real_t maxval = evaluate<0>( a, l );
      {
         real_t val = evaluate<0>( b, r );

         if( val > maxval )
         {
            maxval = val;
            max = b;
         }
      }
      ++r;

      for( interval_t i = l; i != r; ++i )
      {
         std::complex<real_t> localRoots[DEGREE - 1];
         auto pp = deriv.getPolynomialPiece( i );
         int nroots = pp.getRoots( localRoots, tolerance * 1e-2 );

         for( int rt = 0; rt < nroots; ++rt )
         {
            if( std::abs( localRoots[rt].imag() ) * 2 <= tolerance &&
                  localRoots[rt].real() >= std::max( a, getInfimum( i ) ) &&
                  localRoots[rt].real() <= std::min( b, getSupremum( i ) )
              )
            {
               real_t val = evaluate<0>( localRoots[rt].real(), i );

               if( val > maxval )
               {
                  max = localRoots[rt].real();
                  maxval = val;
               }

            }
         }
      }

      return max;
   }


   /**
    * Finds all real roots within the
    * given interval using the given tolerance.
    *
    * \param a          left endpoint of interval
    * \param b          right endpoint of interval
    * \param tolerance  The tolerance used to find the values.
    * \return Vector of all real roots of this curve within the given interval.
    */
   std::vector<real_t> getRealRoots( real_t a, real_t b, real_t tolerance )
   {
      std::vector<real_t> roots;

      interval_t l = findInterval( a, interval_ );
      interval_t r = findInterval( b, l );
      interval_ = r;
      ++r;

      for( interval_t i = l; i != r; ++i )
      {
         std::complex<real_t> localRoots[DEGREE];
         int nroots = getPolynomialPiece( i ).getRoots( localRoots, tolerance * 1e-2 );
         std::size_t size = roots.size();

         for( int rt = 0; rt < nroots; ++rt )
         {
            if( std::abs( localRoots[rt].imag() ) * 2 <= tolerance &&
                  localRoots[rt].real() >= std::max( a, getInfimum( i ) ) &&
                  localRoots[rt].real() <= std::min( b, getSupremum( i ) ) )
               roots.push_back( localRoots[rt].real() );
         }

         std::sort( roots.begin() + size, roots.end() );
      }

      return roots;
   }

   /**
    * Solves an equation with this curve on the left handside
    * and a given polynomial on the right hand side
    * within the given interval using the given tolerance.
    *
    * \param rhs        the right hand side polynomial
    * \param a          left endpoint of interval
    * \param b          right endpoint of interval
    * \param tolerance  The tolerance used to find the values.
    * \return Vector of all solutions to the equation within the given interval.
    */
   template < typename RHS_IMPL >
   std::vector<real_t> solveEquation( const Polynomial<RHS_IMPL> &rhs, real_t a, real_t b, real_t tolerance )
   {
      std::vector<real_t> sol;
      interval_t l = findInterval( a, interval_ );
      interval_t r = findInterval( b, l );
      interval_ = r;
      ++r;

      for( interval_t i = l; i != r; ++i )
      {
         std::complex<real_t> localSol[DEGREE];
         int nsols = getPolynomialPiece( i ).solveEquation( rhs, localSol, tolerance * 1e-2 );
         std::size_t size = sol.size();

         for( int s = 0; s < nsols; ++s )
         {
            if( std::abs( localSol[s].imag() ) * 2 <= tolerance &&
                  localSol[s].real() >= std::max( a, getInfimum( i ) ) &&
                  localSol[s].real() <= std::min( b, getSupremum( i ) )
              )
               sol.push_back( localSol[s].real() );
         }

         std::sort( sol.begin() + size, sol.end() );
      }

      return sol;
   }

   /**
    * Returns the polynomial representing this curve in the given interval.
    *
    * \param interval   The interval to return the polynomial for.
    * \return           The representing polynomial.
    */
   PolynomialPiece<Implementation> getPolynomialPiece( const interval_t interval ) const
   {
      return PolynomialPiece<Implementation>( interval, *this );
   }

private:

   void findMinMax( real_t &minval, real_t &maxval, real_t &minarg, real_t &maxarg, real_t a, real_t b, real_t tolerance ) const
   {
      auto deriv = differentiate( *this );
      a = std::max(a, getInfimum(0));
      b = std::min(b, getSupremum(numIntervals()-1));
      interval_t l = findInterval( a, interval_ );
      interval_t r = findInterval( b, l );
      interval_ = r;
      minarg = a;
      maxarg = b;
      minval = evaluate<0>( a, l );
      maxval = evaluate<0>( b, r );

      if( minval > maxval )
      {
         std::swap( minval, maxval );
         std::swap( minarg, maxarg );
      }

      ++r;

      for( interval_t i = l; i != r; ++i )
      {
         std::complex<real_t> localRoots[DEGREE - 1];
         auto pp = deriv.getPolynomialPiece( i );
         int nroots = pp.getRoots( localRoots, tolerance * 1e-2 );

         for( int rt = 0; rt < nroots; ++rt )
         {
            if( std::abs( localRoots[rt].imag() ) * 2  <= tolerance &&
                  localRoots[rt].real() >= std::max( a, getInfimum( i ) ) &&
                  localRoots[rt].real() <= std::min( b, getSupremum( i ) )
              )
            {
               real_t val = evaluate<0>( localRoots[rt].real(), i );

               if( val < minval )
               {
                  minarg = localRoots[rt].real();
                  minval = val;
               }

               if( val > maxval )
               {
                  maxarg = localRoots[rt].real();
                  maxval = val;
               }
            }
         }
      }
   }

   Implementation &getImpl()
   {
      return *static_cast<Implementation *>( this );
   }

   const Implementation &getImpl() const
   {
      return *static_cast<const Implementation *>( this );
   }

   mutable interval_t interval_;
};

} // namespace spline

#endif
