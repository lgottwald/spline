#include "catch/catch.hpp"
#include <spline/Polynomial.hpp>
#include <spline/SimplePolynomial.hpp>


bool complex_less_pred( std::complex<double> const &x, std::complex<double> const &y )
{
   return x.real() < y.real();
}

TEST_CASE( "Polynomial roots are found correctly for degree 1 to 5", "[spline]" )
{
   double tolerance = 1e-12;
   std::complex<double> roots[4];

   SECTION( "test for degree 1 polynomial" )
   {
      //degree 3 polynomial that should be considered as degree 1 because the leading two coefficients are smaller than tolerance
      double coeffs[] = { -5, 3, 0.0, tolerance / 1.1};
      spline::SimplePolynomial<3, double> poly( coeffs );

      int nroots = poly.getRoots( roots, tolerance );
      REQUIRE( nroots == 1 );
      REQUIRE( roots[0].imag() == 0.0 );
      REQUIRE( roots[0].real() == 5.0 / 3 );

   }

   SECTION( "test for degree 2 polynomial" )
   {
      double coeffs[] = {2, 5, 3};
      spline::SimplePolynomial<2, double> poly( coeffs );
      int nroots = poly.getRoots( roots, tolerance );
      REQUIRE( nroots == 2 );
      REQUIRE( Approx( roots[0].imag() ) == 0 );
      REQUIRE( Approx( roots[1].imag() ) == 0 );
      std::sort( roots, roots + 2, complex_less_pred );
      REQUIRE( Approx( roots[0].real() ) == -1 );
      REQUIRE( Approx( roots[1].real() ) == -2.0 / 3.0 );
   }

   SECTION( "test for degree 3 polynomial" )
   {
      double coeffs[] = {2, -1, -2, 1};
      spline::SimplePolynomial<3, double> poly( coeffs );
      int nroots = poly.getRoots( roots, tolerance );
      REQUIRE( nroots == 3 );
      REQUIRE( Approx( roots[0].imag() ) == 0 );
      REQUIRE( Approx( roots[1].imag() ) == 0 );
      REQUIRE( Approx( roots[2].imag() ) == 0 );

      std::sort( roots, roots + 3, complex_less_pred );
      REQUIRE( Approx( roots[0].real() ) == -1 );
      REQUIRE( Approx( roots[1].real() ) == 1 );
      REQUIRE( Approx( roots[2].real() ) == 2 );
   }

   SECTION( "test for degree 4 polynomial" )
   {
      double coeffs[] = {15, -9, -11, 3, 2};
      spline::SimplePolynomial<4, double> poly( coeffs );
      int nroots = poly.getRoots( roots, tolerance );

      REQUIRE( nroots == 4 );
      REQUIRE( Approx( roots[0].imag() ) == 0 );
      REQUIRE( Approx( roots[1].imag() ) == 0 );
      REQUIRE( Approx( roots[2].imag() ) == 0 );
      REQUIRE( Approx( roots[3].imag() ) == 0 );
      std::sort( roots, roots + 4, complex_less_pred );
      REQUIRE( Approx( roots[0].real() ) == -2.5 );
      REQUIRE( Approx( roots[1].real() ) == -std::sqrt( 3 ) );
      REQUIRE( Approx( roots[2].real() ) == 1 );
      REQUIRE( Approx( roots[3].real() ) == std::sqrt( 3 ) );

   }
   SECTION( "test solve equation on polynomial" )
   {
      double coeffs1[] = { 3, -2, 8, 4 };
      spline::SimplePolynomial<3, double> poly_lhs( coeffs1 );
      double coeffs2[] = { -8, 7, 5 };
      spline::SimplePolynomial<2, double> poly_rhs( coeffs2 );
      int nsols = poly_lhs.solveEquation( poly_rhs, roots, tolerance );
      REQUIRE( nsols == 3 );
      std::sort( roots, roots + 3, complex_less_pred );
      REQUIRE( roots[0].real() == Approx( -2.2725554037 ) );
      REQUIRE( roots[1].real() == Approx( roots[2].real() ) );
      REQUIRE( roots[1].real() == Approx( 0.761278 ) );
      REQUIRE( roots[1].imag() == Approx( -roots[2].imag() ) );
      REQUIRE( std::abs( roots[1].imag() ) == Approx( 0.79407 ) );
   }

}
