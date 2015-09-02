#ifndef _SPLINE_POLYNOMIAL_HPP_
#define _SPLINE_POLYNOMIAL_HPP_

#include <complex>
#include <algorithm>
#include <type_traits>
#include <cpplsq/AutoDiff.hpp>
#include "DegreeTraits.hpp"
#include "FloatingPointTraits.hpp"

#define RESULT(...) decltype(__VA_ARGS__) { return (__VA_ARGS__); }

namespace spline
{

/**
 * A polynomial crtp base class
 */
template< class Implementation >
class Polynomial;

namespace internal
{

/**
 * Compile time factorial to compute highest order coefficient
 * from derivative
 */
constexpr int factorial( int N )
{
   return N > 1 ? N * factorial( N - 1 ) : 1;
}

/**
 * Class representing the difference of two polynomial
 * types.
 */
template < typename Impl1, typename Impl2 >
struct PolynomialDifference : public Polynomial< PolynomialDifference<Impl1, Impl2> >
{
   template<typename T>
   T operator()( const T &x ) const
   {
      return poly1( x ) - poly2( x );
   }

   template<int D, typename T>
   T derivative( const T &x ) const
   {
      return poly1.template derivative<D>( x ) - poly2.template derivative<D>( x );
   }

   PolynomialDifference(
      const Polynomial<Impl1> &poly1,
      const Polynomial<Impl2> &poly2
   ) : poly1( poly1 ), poly2( poly2 ) {}

private:
   const Polynomial<Impl1> &poly1;
   const Polynomial<Impl2> &poly2;
};

/**
 * Find root of degree one polynomial/line analytically.
 */
template < class IMPL,
         int NROOTS = get_degree<IMPL>(),
         int MIN_ITER = 4,
         int MAX_ITER = 100 >
typename std::enable_if< ( NROOTS == 1 ), int >::type
polynomial_roots(
   const Polynomial<IMPL> &poly,
   std::complex<floating_point<IMPL>> *roots,
   const floating_point<IMPL> tolerance )
{
   using REAL = floating_point<IMPL>;
   REAL Slope = poly.derivative( REAL( 0. ) );

   if( std::abs( Slope ) <= tolerance )
      return 0;

   roots[0] = -poly( REAL( 0. ) ) / Slope;
   return 1;
}


/**
 * Find root of degree two polynomial analytically.
 */
template < class IMPL,
         int NROOTS = get_degree<IMPL>(),
         int MIN_ITER = 4,
         int MAX_ITER = 100 >
typename std::enable_if< ( NROOTS == 2 ), int >::type
polynomial_roots(
   const Polynomial<IMPL> &poly,
   std::complex<floating_point<IMPL>> *roots,
   const floating_point<IMPL> tolerance )
{
   using REAL = floating_point<IMPL>;
   using COMPLEX = std::complex<REAL>;

   REAL QuadraticCoeff = poly.template derivative<2>( REAL( 0. ) ) / REAL( 2 );

   if( std::abs( QuadraticCoeff ) <= tolerance )
      return polynomial_roots < IMPL, NROOTS - 1, MIN_ITER, MAX_ITER > ( poly, roots, tolerance );

   REAL LinearCoeff = poly.derivative( REAL( 0 ) );
   REAL Constant = poly( REAL( 0 ) );
   COMPLEX q = REAL( -0.5 ) * ( LinearCoeff + std::copysign( REAL( 1.0 ), LinearCoeff ) *
                                std::sqrt( COMPLEX( LinearCoeff * LinearCoeff - 4 * QuadraticCoeff * Constant ) ) );
   roots[0] = q / QuadraticCoeff;
   roots[1] = Constant / q;
   return 2;
}

/**
 * Find root of higher degree polynomials iteratively
 * using Durand–Kerner method.
 */
template < class IMPL,
         int NROOTS = get_degree<IMPL>(),
         int MIN_ITER = 4,
         int MAX_ITER = 100 >
typename std::enable_if < ( NROOTS > 2 ), int >::type
polynomial_roots(
   const Polynomial<IMPL> &poly,
   std::complex<floating_point<IMPL>> *roots,
   const floating_point<IMPL> tolerance )
{
   using REAL = floating_point<IMPL>;
   using COMPLEX = std::complex<REAL>;
   REAL HighestOrderCoeff = poly.template derivative<NROOTS>( REAL( 0 ) ) / factorial( NROOTS );

   //If highest order coefficient is zero within tolerance call
   //function to find roots for lower degree polynomial.
   if( std::abs( HighestOrderCoeff ) <= tolerance )
      return polynomial_roots < IMPL, NROOTS - 1, MIN_ITER, MAX_ITER > ( poly, roots, tolerance );

   roots[0] = REAL( 1.0 );

   for( int i = 1; i < NROOTS; ++i )
      roots[i] = roots[i - 1] * COMPLEX( REAL( 0.4 ), REAL( 0.9 ) );

   //first perform MIN_ITER many iterations without checking for convergence
   //since MIN_ITER is known at compiletime these loop iterations are likely
   //to be unrolled by the compiler
   for( unsigned n = 0; n < MIN_ITER; ++n )
   {
      for( int i = 0; i < NROOTS; ++i )
      {
         COMPLEX q = HighestOrderCoeff;

         for( int j = 1; j < NROOTS; ++j )
            q *= ( roots[i] - roots[( i + j ) % NROOTS] );

         roots[i] -= poly( roots[i] ) / q;
      }
   }

   //Now keep iterating until convergence
   for( unsigned n = MIN_ITER; n < MAX_ITER; ++n )
   {
      bool toleranceReached = true;

      for( int i = 0; i < NROOTS; ++i )
      {
         COMPLEX q = HighestOrderCoeff;

         for( int j = 1; j < NROOTS; ++j )
            q *= ( roots[i] - roots[( i + j ) % NROOTS] );

         COMPLEX val = poly( roots[i] );

         roots[i] -= val / q;
         toleranceReached = toleranceReached && std::abs( val.real() ) + std::abs( val.imag() ) <= tolerance;
      }

      if( toleranceReached )
         break;
   }

   return NROOTS;
}


}

/**
 * Specialization of DegreeTraits for polynomial difference type.
 */
template <typename Impl1, typename Impl2>
struct DegreeTraits<internal::PolynomialDifference<Impl1, Impl2>>
{
   constexpr static int value = get_degree<Impl1>() > get_degree<Impl2>() ? get_degree<Impl1>() : get_degree<Impl2>();
};

/**
 * Specialization of FloatingPointTraits for polynomial difference type.
 */
template <typename Impl1, typename Impl2>
struct FloatingPointTraits<internal::PolynomialDifference<Impl1, Impl2>>
{
   using type = decltype( floating_point<Impl1>() - floating_point<Impl2>() );
};

/**
 * Compile time function to compute factor for evaluating
 * the derivative of polynomial.
 */
constexpr int deriv_factor( int deriv, int power )
{
   return deriv ? deriv_factor( deriv - 1, power - 1 ) * power : 1;
}

/**
 * Evaluate a polynomial by the given coefficients. If it is
 * differentiated degree times then the result is zero
 * and this sfinae overload will be used.
 */
template < int DEGREE, int DERIVATIVE = 0, int POWER = DERIVATIVE, typename REAL, typename T,
         typename std::enable_if < ( DERIVATIVE > DEGREE ), int >::type = 0 >
               REAL evaluatePolynomial( const REAL *coeffs, const T &x )
{
   return 0;
}

/**
 * Evaluate a polynomial by the given coefficients. This
 * sfinae overload will be used for evaluating the term
 * with highest power and no recursive call will be made.
 */
template < int DEGREE, int DERIVATIVE = 0, int POWER = DERIVATIVE, typename REAL, typename T,
typename std::enable_if < ( DERIVATIVE <= DEGREE && POWER == DEGREE ), int >::type = 0 >
REAL evaluatePolynomial( const REAL *coeffs, const T &x )
{
   return  deriv_factor( DERIVATIVE, POWER ) * coeffs[POWER] ;
}


template < int DEGREE, int DERIVATIVE = 0, int POWER = DERIVATIVE, typename REAL, typename T,
typename std::enable_if < ( DERIVATIVE <= DEGREE && POWER < DEGREE ), int >::type = 0 >
cpplsq::ValueType<T> evaluatePolynomial( const REAL *coeffs, const T &x )
{
   return ( deriv_factor( DERIVATIVE, POWER ) * coeffs[POWER] ) +
   x * evaluatePolynomial < DEGREE, DERIVATIVE, POWER + 1, REAL, T > ( coeffs, x );
}


template <typename IMPL1, typename IMPL2>
internal::PolynomialDifference<IMPL1, IMPL2> operator-(
   const Polynomial<IMPL1> &poly1,
   const Polynomial<IMPL2> &poly2
)
{
   return internal::PolynomialDifference<IMPL1, IMPL2>( poly1, poly2 );
}

template< class Implementation >
class Polynomial
{
   public:
   static_assert( has_degree<Implementation>(), "DegreeTraits required for implementation of Polynomial" );
   static_assert( floating_point_type_defined<Implementation>(), "FloatingPointTraits required for implementation of Polynomial" );
   using REAL = floating_point<Implementation>;
   // ========== Interface required for implementation ==========
   /**
    * Evaluate the this polynomial. Must be implemented by all
    * polynomial types deriving from the crtp Polynomial base class.
    *
    * \param x    x value to evaluate this polynomial
    * \return     the value of the polynomial at the given value.
    */
   template<typename T>
   T operator()( const T &x ) const
{
   return getImpl()( x );
}

/**
 * Evaluate derivatives of this polynomial. Must be implemented by all
 * polynomial types deriving from the crtp Polynomial base class.
 *
 * \param x    x value to evaluate this polynomial
 * \tparam D   the derivative to evaluate. E.g. D=1 evaluates the first derivative.
 * \return     the value of this polynomials derivative at the given value.
 */
template<int D = 1, typename T>
T derivative( const T &x ) const
{
   return getImpl().template derivative<D>( x );
}


// ============== functionality implemented here ==============

/**
 * Compute all (complex) roots of this polynomial
 * and store them in the given array. Returns the
 * number of roots found, i.e. the actual degree of
 * the polynomial which might be
 * lower than the degree known at compile time if the
 * highest order coefficient is zero. The degree is
 * equal to the number of roots found.
 *
 * \param roots      array to store the roots having size of at least the degree.
 * \param tolerance  epsilon/tolerance for computing the roots
 * \return number of complex roots found and thus the actual degree of the polynomial
 */
int getRoots(
   std::complex<REAL> *roots,
   REAL tolerance = 1e-12
) const
{
   return internal::polynomial_roots( *this, roots, tolerance );
}


/**
 * Solve a polynomial equation having this polynomial
 * at the left hand side and the given polynomial
 * at the right hand side. Stores the solutions in the
 * given array and returns the numbre of solutions
 * found. The given array size should be at least
 * the degree.
 *
 * \param rhs        Right hand side polynomial.
 * \param solution   Array to store solutions with size of at least the degree.
 * \param tolerance  epsilon/tolerance for computing the roots
 * \return number of complex solutions found.
 */
template<typename OtherImpl>
int solveEquation(
   const Polynomial<OtherImpl> &rhs,
   std::complex<REAL> *solution,
   const REAL tolerance = REAL( 1e-12 )
) const
{
   return ( ( *this ) - rhs ).getRoots( solution, tolerance );
}

/**
 * Compute the minimum of this polynomial in the given interval.
 *
 * \param a          Left endpoint of the interval.
 * \param b          Right endpoint of the interval.
 * \param tolerance  epsilon/tolerance for computing the roots
 * \return The x value in the given interval where this polynomial attains the lowest value
 */
REAL findMin( REAL a, REAL b, REAL tolerance = REAL( 1e-12 ) )
{
   auto cmp = [this]( const std::complex<REAL> &x, const std::complex<REAL> &y )
   {
      return ( *this )( x.real() ) < ( *this )( y.real() );
   };

   return findExtremePoint( a, b, std::move( cmp ), tolerance );
}

/**
 * Compute the maximum of this polynomial in the given interval.
 *
 * \param a          Left endpoint of the interval.
 * \param b          Right endpoint of the interval.
 * \param tolerance  epsilon/tolerance for computing the roots
 * \return The x value in the given interval where this polynomial attains the highest value
 */
REAL findMax( REAL a, REAL b, REAL tolerance = REAL( 1e-12 ) )
{
   auto cmp = [this]( const std::complex<REAL> &x, const std::complex<REAL> &y )
   {
      return ( *this )( x.real() ) > ( *this )( y.real() );
   };

   return findExtremePoint( a, b, cmp, tolerance );
}

private:
template<typename CMP>
REAL findExtremePoint( REAL a, REAL b, CMP cmp, REAL tolerance = REAL( 1e-12 ) )
{
   std::complex<REAL> roots[get_degree<Implementation>() - 1];
   int nroots = internal::polynomial_roots( differentiate( *this ), roots, tolerance );
   REAL m = cmp( a, b ) ? a : b;
   auto end = std::remove_if( roots, roots + nroots, [ = ]( const std::complex<REAL> &x )
   {
      return x.imag() > tolerance;
   } );
   REAL m2 = std::min_element( roots, end, cmp )->real();
   return cmp( m, m2 ) ? m : m2;
}

Implementation &getImpl()
{
   return *static_cast<Implementation *>( this );
}

const Implementation &getImpl() const
{
   return *static_cast<const Implementation *>( this );
}

};



} // namespace spline

#endif
