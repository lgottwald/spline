#ifndef _PIECEWISE_POLYNOMIAL_IMPL_HPP_
#define _PIECEWISE_POLYNOMIAL_IMPL_HPP_

#include "UnivariateSplineCurve.hpp"
#include "Polynomial.hpp"
#include "Infinity.hpp"
#include <type_traits>
#include <simd/alloc.hpp>
#include <simd/pack.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>

namespace spline
{
using boost::array;

template<int DEGREE, typename REAL>
class PiecewisePolynomial;

template<int DEGREE, typename REAL>
struct DegreeTraits<PiecewisePolynomial<DEGREE, REAL>>
{
   constexpr static int value = DEGREE;
};

template<int DEGREE, typename REAL>
struct FloatingPointTraits<PiecewisePolynomial<DEGREE, REAL>>
{
   using type = REAL;
};

template<int DEGREE, typename REAL>
struct IntervalTypeTraits<PiecewisePolynomial<DEGREE, REAL>>
{
   using type = std::size_t;
};

/**
 * Class implementing the UnivariateSplineCurve crtp base class
 * using the piecewise polynomial representation internally.
 * This yields faster evaluation but is less stable then the BSpline
 * representation.
 */
template<int DEGREE, typename REAL>
class PiecewisePolynomial : public UnivariateSplineCurve< PiecewisePolynomial<DEGREE, REAL> >
{
public:
   using interval_t = interval_type<PiecewisePolynomial<DEGREE, REAL>>;

   PiecewisePolynomial() {}

   /**
    * Create a PiecewisePolynomial from UnivariateSplineCurve
    * using taylor expansion to compute the  coefficients. TODO currently
    * implemented assuming degree equal to 3.
    */
   template<typename IMPL>
   PiecewisePolynomial( const UnivariateSplineCurve< IMPL > &curve ) : breakpoints( curve.numIntervals() + simd::pack_size<REAL>() ), coeffs( curve.numIntervals() )
   {
      /*
       * Now convert the SplineCurve to its piecewise polynomial representation
       * by computing the coefficients of the polynomial pieces using taylor expansion.
       */

      for( interval_t i = 0; i < coeffs.size(); ++i )
      {
         REAL a = curve.getInfimum( i );
         breakpoints[i] = a;
         REAL f0 = curve.template evaluate<0>( a, i );
         REAL f1 = curve.template evaluate<1>( a, i );
         REAL f2 = curve.template evaluate<2>( a, i );
         REAL f3 = curve.template evaluate<3>( a, i );

         REAL a2 = a * a;

         coeffs[i][0] = -a * a2 * f3 / 6 + a2 * f2 / 2 - a * f1 + f0;
         coeffs[i][1] = a2 * f3 / 2 - a * f2 + f1;
         coeffs[i][2] = -a * f3 / 2 + f2 / 2;
         coeffs[i][3] = f3 / 6;

      }

      //add sentinels for array of breakpoints to speedup the search for an interval
      breakpoints[0] = -infinity<REAL>();

      for( interval_t i = breakpoints.size() - simd::pack_size<REAL>(); i < breakpoints.size(); ++i )
      {
         breakpoints[i] = infinity<REAL>();
      }
   }

   template<typename Iter1, typename Iter2>
   PiecewisePolynomial( typename std::enable_if <
                        std::is_same< typename std::iterator_traits<Iter1>::value_type, REAL>::value, Iter1
                        >::type breakPointBegin,
                        Iter1 breakPointEnd,
                        typename std::enable_if <
                        std::is_same < typename std::iterator_traits<Iter2>::value_type, array < REAL, DEGREE + 1 > >::value, Iter2
                        >::type coeffBegin,
                        Iter2 coeffEnd ) :  breakpoints( breakPointBegin, breakPointEnd ), coeffs( coeffBegin, coeffEnd )
   {
      breakpoints[0] = -infinity<REAL>();

      breakpoints.back() = infinity<REAL>();

      for( int i = 0; i < simd::pack_size<REAL>() - 1; ++i )
      {
         breakpoints.push_back( infinity<REAL>() );
      }
   }

   template<typename Iter1, typename Iter2>
   PiecewisePolynomial( typename std::enable_if <
                        std::is_same< typename std::iterator_traits<Iter1>::value_type, REAL>::value, Iter1
                        >::type breakPointBegin,
                        Iter1 breakPointEnd,
                        typename std::enable_if <
                        std::is_same< typename std::iterator_traits<Iter2>::value_type, REAL >::value, Iter2
                        >::type breakValBegin,
                        Iter2 breakValEnd ) : breakpoints( ( breakPointEnd - breakPointBegin ) + simd::pack_size<REAL>() + 1 )
   {
      int nbreaks = ( breakPointEnd - breakPointBegin );

      breakpoints[0] = -infinity<REAL>();
      std::copy( breakPointBegin, breakPointEnd, &breakpoints[1] );

      for( std::size_t i = nbreaks + 1; i < breakpoints.size(); ++i )
         breakpoints[i] = infinity<REAL>();

      /* Piecewise linear should stay constant outside the given interval */
      coeffs.resize( nbreaks + 1 );
      coeffs[0][0] = *breakValBegin;

      for( int j = 1; j < DEGREE + 1; ++j )
         coeffs[0][j] = REAL( 0 );

      for( std::size_t i = 1; i < coeffs.size() - 1; ++i )
      {
         REAL x1 = breakpoints[i];
         REAL x2 = breakpoints[i + 1];
         REAL y1 = *( breakValBegin + i - 1 );
         REAL y2 = *( breakValBegin + i );
         coeffs[i][1] = ( y2 - y1 ) / ( x2 - x1 );
         coeffs[i][0] = y1 - coeffs[i][1] * x1;

         for( int j = 2; j < DEGREE + 1; ++j )
            coeffs[i][j] = REAL( 0 );
      }

      coeffs.back()[0] = *( breakValEnd - 1 );

      for( int j = 1; j < DEGREE + 1; ++j )
         coeffs.back()[j] = REAL( 0 );
   }



   REAL getInfimumImpl( const interval_t interval ) const
   {
      return breakpoints[interval];
   }

   REAL getSupremumImpl( const interval_t interval ) const
   {
      return breakpoints[interval + 1];
   }

   std::size_t numIntervalsImpl() const
   {
      return coeffs.size();
   }

   interval_t findIntervalImpl( const REAL x, const interval_t hint ) const
   {
      if( x >= breakpoints[hint + 1] )
      {
         //search forward
         simd::pack<REAL> key( x );
         std::size_t simdsize = simd::prev_size<REAL>( breakpoints.size() );
         interval_t i;

         for( i = simd::prev_size<REAL>( hint ); i < simdsize; i += simd::pack_size<REAL>() )
         {
            int mask = simd::movemask( key < simd::aligned_load( &breakpoints[i] ) );

            if( mask )
               return i - 1 + simd::pack_size<REAL>() - __builtin_popcount( mask );
         }

         while( 1 )
         {
            if( x < breakpoints[i + 1] )
               return i;

            ++i;
         }
      }
      else if( x < breakpoints[hint] )
      {
         //search backward
         simd::pack<REAL> key( x );

         for( interval_t i = simd::next_size<REAL>( hint ); 1; i -= simd::pack_size<REAL>() )
         {
            int mask = simd::movemask( key < simd::aligned_load( &breakpoints[i] ) );

            if( !( mask & 1 ) )
               return i - 1 + simd::pack_size<REAL>() - __builtin_popcount( mask );
         }
      }
      else
      {
         //hint is correct
         return hint;
      }
   }

   template<int D, typename T>
   cpplsq::ValueType<T> evaluate( const T &x, const interval_t interval ) const
   {
      return evaluatePolynomial<DEGREE, D>( coeffs[interval].data(), x );
   }

   std::vector < array < REAL, DEGREE + 1 > > & getCoefficients()
   {
      return coeffs;
   }

   const std::vector < array < REAL, DEGREE + 1 > > & getCoefficients() const
   {
      return coeffs;
   }

   simd::aligned_vector<REAL> &getBreakpoints()
   {
      return breakpoints;
   }

   const simd::aligned_vector<REAL> &getBreakpoints() const
   {
      return breakpoints;
   }

private:

   friend class boost::serialization::access;
   template<class Archive>
   void serialize( Archive &ar, const unsigned int version )
   {
      ar &breakpoints;
      ar &coeffs;
   }

   simd::aligned_vector<REAL> breakpoints;
   std::vector < array < REAL, DEGREE + 1 > > coeffs;
};


} // namespace spline

#endif
