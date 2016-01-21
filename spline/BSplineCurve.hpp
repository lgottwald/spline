#ifndef B_SPLINE_CURVE_HPP
#define B_SPLINE_CURVE_HPP

#include <algorithm>
#include <vector>
#include <cassert>
#include <simd/alloc.hpp>
#include "UnivariateSplineCurve.hpp"
#include <boost/serialization/vector.hpp>
#include "Infinity.hpp"
#include <cmath>
#include <type_traits>

namespace spline
{

template < int DEGREE, typename REAL >
class BSplineCurve;

/**
 * Specialization of degree traits
 * for BSplineCurve
 */
template < int DEG, typename REAL >
struct DegreeTraits<BSplineCurve<DEG, REAL>>
{
   constexpr static int value = DEG;
};

/**
 * Specialization of FloatingPointTraits
 * for BSplineCurve. Used to retrieve
 * the type used for floating point numbers
 * in a BSplineCurve.
 */
template < int DEG, typename REAL >
struct FloatingPointTraits<BSplineCurve<DEG, REAL>>
{
   using type = REAL;
};

/**
 * Specialization of IntervalTypeTraits
 * for BSplineCurve. Used to retrieve
 * the type that identifies an interval
 * of the BSplineCurve.
 */
template < int DEG, typename REAL >
struct IntervalTypeTraits<BSplineCurve<DEG, REAL>>
{
   using type = long;
};

/**
 * This class implements a univariate spline curve represented as
 * a linear combination of bsplines. Uses the CRTP base class
 * spline::UnivariateSplineCurve and implements the required functionality.
 */
template < int DEGREE, typename REAL >
class BSplineCurve : public UnivariateSplineCurve< BSplineCurve<DEGREE, REAL> >
{
public:
   using interval_t = interval_type<BSplineCurve<DEGREE, REAL>>;

   /**
    * Create a BSplineCurve with given knots and coefficients.
    *
    * \param knotsCont       container with the knots
    * \param coeffsCont      container with the coefficients
    */
   template<typename Cont1, typename Cont2>
   BSplineCurve( const Cont1 &knotsCont, const Cont2 &coeffsCont ) :
      UnivariateSplineCurve< BSplineCurve<DEGREE, REAL> >()
   {
      std::size_t n = simd::next_size<REAL>(knotsCont.size()+2*DEGREE + 2);
      knots.reserve(n);

      REAL a = 2 * (*knotsCont.begin()) - *(knotsCont.begin() + 1);
      REAL b = 2 * (*(knotsCont.end() - 1)) - *(knotsCont.end() - 2);

      knots.emplace_back(-infinity<REAL>());

      for( int i = 0; i < DEGREE; ++i )
         knots.emplace_back(a);

      for( auto & knt : knotsCont )
         knots.emplace_back(knt);

      for( int i = 0; i < DEGREE; ++i )
         knots.emplace_back(b);

      for( auto i = knots.size(); i <= n; ++i )
         knots.emplace_back(infinity<REAL>());

      coeffs.reserve(coeffsCont.size() + 2);
      coeffs.emplace_back(*(coeffsCont.begin()));

      for( auto & cf : coeffsCont )
         coeffs.emplace_back(cf);

      coeffs.emplace_back( *(coeffsCont.end() - 1) );
   }

   /**
    * Create a BSplineCurve with space for nknots many knots
    * defined on the given interval.
    *
    * \param a       left endpoint of interval
    * \param b       right endpoint of interval
    * \param nknots  number of  knots
    */
   BSplineCurve( REAL a, REAL b, std::size_t nknots )
   {
      setKnotInterval( a, b, nknots );
   }

   BSplineCurve() = default;
   BSplineCurve( BSplineCurve<DEGREE, REAL> && ) = default;
   BSplineCurve( const BSplineCurve<DEGREE, REAL> & ) = default;
   BSplineCurve<DEGREE, REAL> &operator=( const BSplineCurve<DEGREE, REAL> & ) = default;
   BSplineCurve<DEGREE, REAL> &operator=( BSplineCurve<DEGREE, REAL> && ) = default;

   /**
    * Initialize the knot vector for knots between a and b with nknots many inner knots
    * and allocate enough space for the coefficients.
    * \param a       left endpoint of interval
    * \param b       right endpoint of interval
    * \param nknots  number of knots
    */
   void setKnotInterval( REAL a, REAL b, std::size_t nknots )
   {
      //initialize knot vector with sentinels and double knots
      size_t n = simd::next_size<REAL>( nknots + DEGREE * 2 );
      knots.resize( n );
      knots[0] = -infinity<REAL>();

      for( auto i = 1; i < DEGREE; ++i )
         knots[i] = a;

      for( auto i = nknots + DEGREE; i < nknots + 2 * DEGREE - 1; ++i )
         knots[i] = b;

      for( auto i = nknots + 2 * DEGREE - 1; i < n; ++i )
         knots[i] = infinity<REAL>();

      coeffs.resize( nknots + DEGREE - 1 );
   }

   /**
    * Gets reference to i'th inner knot
    */
   REAL &getKnot( const std::size_t i )
   {
      return knots[i + DEGREE];
   }


   /**
    * Gets const reference to i'th inner knot
    */
   const REAL &getKnot( const std::size_t i ) const
   {
      return knots[i + DEGREE];
   }

   /**
    * Gets number of inner knots
    */
   std::size_t getNumKnots() const
   {
      return coeffs.size() - DEGREE + 1;
   }

   /**
    * Returns reference to vector of coefficients of this bspline curve.
    */
   std::vector<REAL> &getCoefficients()
   {
      return coeffs;
   }


   /**
    * Returns const reference to vector of coefficients of this bspline curve.
    */
   const std::vector<REAL> &getCoefficients() const
   {
      return coeffs;
   }

   /**
    * Choose coefficients for current knot vector
    * to yield a variation diminishing approximation
    * for the given function. This approximation
    * has generalizes many of the properties of piecewise
    * linear approximation to higher degree splines.
    * E.g. the error of the approximation will converge to
    * zero if the number of knots goes to infinity and it
    * preserves the bounds and the shape of the given function.
    */
   template<typename Function>
   void approximateVariationDiminishing( const Function &f )
   {
      for( std::size_t i = 0; i < coeffs.size(); ++i )
      {
         REAL c = 0.0;

         for( std::size_t j = i + 1; j <= i + DEGREE; ++j )
            c += knots[j];

         coeffs[i] = f( c / DEGREE );
      }
   }

   /**
    * Implementation of getInfimum. Returns the infimum of
    * the given interval.
    */
   REAL getInfimumImpl( const interval_t interval ) const
   {
      return knots[interval + DEGREE];
   }

   /**
    * Implementation of getSupremum. Returns the supremum of
    * the given interval.
    */
   REAL getSupremumImpl( const interval_t interval ) const
   {
      return knots[interval + DEGREE + 1];
   }

   /**
    * Implementation of getInfimum. Returns the infimum of
    * the first interval.
    */
   REAL getInfimumImpl() const
   {
      return knots[DEGREE];
   }

   /**
    * Implementation of getSupremum. Returns the supremum of
    * the last interval.
    */
   REAL getSupremumImpl() const
   {
      return knots[coeffs.size()];
   }

   /**
    * Implementation of numIntervals. Returns the number
    * of intervals.
    */
   std::size_t numIntervalsImpl() const
   {
      return coeffs.size() - DEGREE;
   }


   /**
    * Implementation of findInterval. Uses a linear
    * search starting from the hint and dispatches to call
    * an implementation using simd instruction in order to
    * compare multiple values at once if supported.
    */
   template<typename T>
   interval_t findIntervalImpl( const T &x, const interval_t hint ) const
   {
      interval_t i = findIntervalSelectImpl<T>( x, hint );
      assert( ( x >= getKnot( i ) && x < getKnot( i + 1 ) ) || // no extrapolation
              ( x >= getKnot( getNumKnots() - 1 ) && i == interval_t( getNumKnots() - 2 ) ) || // extrapolation of last interval
              ( x <= getKnot( 0 ) && i == 0 ) || // extrapolation of first interval
              x != x || getKnot( i ) != getKnot( i ) || getKnot( i + 1 ) != getKnot( i + 1 ) ); // nan
      return i;
   }

   /**
    * Linear search using simd instructions starting from the hint
    */
   template<typename T>
   interval_t findIntervalSelectImpl( const typename std::enable_if<simd::supported<REAL>::value, T>::type x, const interval_t hint ) const
   {
      interval_t idx = hint + DEGREE;

      if( x >= knots[idx + 1] )
      {
         //in this case search forward
         simd::pack<T> key( x );
         std::size_t simdsize = simd::prev_size<T>( knots.size() );
         std::size_t i;

         for( i = simd::prev_size<T>( idx ); i < simdsize; i += simd::pack_size<T>() )
         {
            int mask = simd::movemask( key < simd::aligned_load( &knots[i] ) );

            if( mask )
            {
               idx =  i - 1 + interval_t( simd::pack_size<REAL>() ) - __builtin_popcount( mask ) - DEGREE;
               return std::min( interval_t( numIntervalsImpl() - 1 ), idx );
            }
         }

         while( 1 )
         {
            if( x < knots[i + 1] )
            {
               idx = i - DEGREE;
               return std::min( interval_t( numIntervalsImpl() - 1 ), idx );
            }

            ++i;
         }
      }
      else if( x < knots[idx] )
      {
         //in this case search backward
         simd::pack<T> key( x );

         for( interval_t i = simd::next_size<T>( idx ); 1; i -= simd::pack_size<T>() )
         {
            int mask = simd::movemask( key < simd::aligned_load( &knots[i] ) );

            if( !( mask & 1 ) )
            {
               idx = i - 1 + interval_t( simd::pack_size<T>() ) - __builtin_popcount( mask ) - DEGREE;
               return std::max( interval_t( 0 ), idx );
            }
         }
      }
      else
      {
         //hint is correct
         return hint;
      }
   }

   template <typename T>
   interval_t findIntervalSelectImpl( const typename std::enable_if < !simd::supported<T>::value, T >::type &x, interval_t i ) const
   {
      i += DEGREE;

      while( i < interval_t( coeffs.size() - 1 ) && x >= knots[i + 1] )
      {
         ++i;
      }

      while( i > DEGREE && x < knots[i] )
      {
         --i;
      }

      i -= DEGREE;
      return i;
   }

   /**
    * Evaluate this BSplineCurve.
    *
    * \param x         x value for evaluating the spline
    * \param interval  interval for which spline will be evaluated
    * \tparam D        derivative to evaluate.
    */
   template<int D, typename T, typename TVAL = cpplsq::ValueType<T>>
   TVAL evaluate( const T& x, const interval_t interval ) const
   {

      if( DEGREE == 0 )
         return coeffs[interval];

      TVAL y( 0.0 );

      /* evaluate bsplines for interval and xvalue */
      TVAL basis[DEGREE + 1 - D];
      evaluateBSplines < DEGREE - D > ( x, interval, basis );

      /* if D>0 then differentiate coefficients and evaluate with the
         differentiated ones */
      if( D )
      {
         /* Set up an array of the coefficients of the basis splines that are non-zero */
         REAL dcoeffs[DEGREE + 1];
         std::copy( coeffs.begin() + interval, coeffs.begin() + interval + DEGREE + 1, dcoeffs );

         /* Differentiate the coefficients D times */
         for( int d = 0; d < D; ++d )
         {
            for( int j = 1; j < DEGREE + 1 - d; ++j )
            {
               dcoeffs[j - 1] = ( DEGREE - d ) *
                                ( dcoeffs[j] - dcoeffs[j - 1] ) /
                                ( knots[interval + DEGREE + j] - knots[interval + j + d] );
            }
         }

         for( int j = 0; j < DEGREE + 1 - D; ++j )
            y += dcoeffs[j] * basis[j];
      }
      else
      {
         /* evaluate with coefficients */
         for( int j = 0; j < DEGREE + 1; ++j )
            y += coeffs[interval + j] * basis[j];
      }

      return y;
   }

   /**
    * evaluates the DEG+1 many nonzero bsplines of degree DEG
    * and stores the results in the given array B
    */
   template<int DEG, typename T, typename TVAL = cpplsq::ValueType<T>>
   void evaluateBSplines( const T &x, const interval_t i, TVAL B[DEG + 1] ) const
   {
      B[0] = 1.0;
      interval_t index = i + DEGREE;

      for( int k = 1; k < DEG + 1; ++k )
      {
         B[k] = B[k - 1] * W( index, k, x );

         for( int j = k - 1; j >= 1; --j )
         {
            interval_t idx = index - k + j;
            B[j] = B[j - 1] * W( idx, k, x ) + B[j] * ( 1.0 - W( idx + 1, k, x ) );
         }

         B[0] = B[0] * ( 1.0 - W( index - k + 1, k, x ) );
      }
   }

   /**
    * evaluates the nonzero bsplines up to degree DEG
    * and stores the results in the given array B
    * where B[0][0] contains the nonzero bespline of degree 0 which is the constant 1
    * B[1][0] B[1][1] will contain the two nonzero bsplines of degree 1 and so on.
    */
   template<int DEG, typename T, typename TVAL = cpplsq::ValueType<T>>
   void evaluateAllBSplines( const T &x, const interval_t i, T B[DEG + 1][DEG + 1] ) const
   {
      B[0][0] = T( 1.0 );
      interval_t index = i + DEGREE;

      for( int k = 1; k < DEG + 1; ++k )
      {
         B[k][k] = B[k - 1][k - 1] * W( index, k, x );

         for( int j = k - 1; j >= 1; --j )
         {
            interval_t idx = index - k + j;
            B[k][j] = B[k - 1][j - 1] * W( idx, k, x ) + B[k - 1][j] * ( T( 1.0 ) - W( idx + 1, k, x ) );
         }

         B[k][0] = B[k - 1][0] * ( T( 1.0 ) - W( index - k + 1, k, x ) );
      }
   }
private:

   template<typename T, typename TVAL = cpplsq::ValueType<T>>
   TVAL W( const interval_t i, const int k, const T &x ) const
   {
      return ( x - knots[i] ) / ( knots[i + k] - knots[i] );
   }

   friend class boost::serialization::access;
   template<class Archive>
   void serialize( Archive &ar, const unsigned int version )
   {
      ar &knots;
      ar &coeffs;
   }

   simd::aligned_vector<REAL> knots;
   std::vector<REAL> coeffs;
};


} // namespace spline

#endif // B_SPLINE_CURVE_HPP

