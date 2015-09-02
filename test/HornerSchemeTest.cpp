#include "catch/catch.hpp"
#include <cmath>
#include <spline/Polynomial.hpp>
#include <boost/array.hpp>

using boost::array;

template<int DIM>
double polynomial_evaluation_reference( const array < double, DIM + 1 > &coeffs, double x )
{
   double res = coeffs[0];

   for( int i = 1; i < DIM + 1; ++i )
      res += coeffs[i] * std::pow( x, i );

   return res;
}

TEST_CASE( "Horner scheme evaluates correctly", "[spline]" )
{
   array<double, 5> coeffs = {1, 2, 3, 4, 5};
   array<double, 4> dcoeffs =  {coeffs[1], 2 * coeffs[2], 3 * coeffs[3], 4 * coeffs[4]};
   array<double, 3> d2coeffs =  {dcoeffs[1], 2 * dcoeffs[2], 3 * dcoeffs[3]};
   double x = M_PI;
   SECTION( "evaluation works for polynomial itself" )
   {
      double value = spline::evaluatePolynomial<4>( coeffs.data(), x );
      double reference = polynomial_evaluation_reference<4>( coeffs, x );
      REQUIRE( value == Approx( reference ) );
   }

   SECTION( "evaluation works for first derivative" )
   {
      double value = spline::evaluatePolynomial<4, 1>( coeffs.data(), x );
      double reference = polynomial_evaluation_reference<3>( dcoeffs, x );

      REQUIRE( value == Approx( reference ) );
   }

   SECTION( "evaluation works for second derivative" )
   {
      double value = spline::evaluatePolynomial<4, 2>( coeffs.data(), x );
      double reference = polynomial_evaluation_reference<2>( d2coeffs, x );
      REQUIRE( value == Approx( reference ) );
   }
}
