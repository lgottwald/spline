#include "catch/catch.hpp"
#include <spline/KnotParameterTransform.hpp>
#include <spline/PiecewisePolynomial.hpp>
#include <spline/BSplineCurve.hpp>
#include <cmath>
#include <vector>


TEST_CASE( "Variation diminishing coefficients are computed correctly", "[spline]" )
{
   spline::BSplineCurve<3, double> curve( 1, 6, 5 );

   for( int i = 0; i < 5; ++i )
      curve.getKnot( i ) = i + 1;

   curve.approximateVariationDiminishing( []( double x )
   {
      return x;
   } );
   auto &coeffs = curve.getCoefficients();
   REQUIRE( coeffs[0] == Approx( ( ( curve.getKnot( 0 ) + curve.getKnot( 0 ) + curve.getKnot( 0 ) ) / 3.0 ) ) );
   REQUIRE( coeffs[1] == Approx( ( ( curve.getKnot( 0 ) + curve.getKnot( 0 ) + curve.getKnot( 1 ) ) / 3.0 ) ) );
   REQUIRE( coeffs[2] == Approx( ( ( curve.getKnot( 0 ) + curve.getKnot( 1 ) + curve.getKnot( 2 ) ) / 3.0 ) ) );
   REQUIRE( coeffs[3] == Approx( ( ( curve.getKnot( 1 ) + curve.getKnot( 2 ) + curve.getKnot( 3 ) ) / 3.0 ) ) );
   REQUIRE( coeffs[4] == Approx( ( ( curve.getKnot( 2 ) + curve.getKnot( 3 ) + curve.getKnot( 4 ) ) / 3.0 ) ) );
   REQUIRE( coeffs[5] == Approx( ( ( curve.getKnot( 3 ) + curve.getKnot( 4 ) + curve.getKnot( 5 ) ) / 3.0 ) ) );
   REQUIRE( coeffs[6] == Approx( ( ( curve.getKnot( 4 ) + curve.getKnot( 5 ) + curve.getKnot( 6 ) ) / 3.0 ) ) );
}

TEST_CASE( "Parameter log transformation is computed and inversed correctly", "[spline]" )
{
   spline::BSplineCurve<3, double> curve( 0, 5, 6 );

   for( int i = 0; i < 6; ++i )
      curve.getKnot( i ) = i;

   spline::BSplineCurve<3, double> curve2( 0, 5, 6 );

   simd::aligned_vector<double> params;
   spline::parameter_log_transform( curve, params );
   REQUIRE( params.size() == 4 );
   REQUIRE( Approx( params[0] ) == 0.0 );
   REQUIRE( Approx( params[1] ) == 0.0 );
   REQUIRE( Approx( params[2] ) == 0.0 );
   REQUIRE( Approx( params[3] ) == 0.0 );
   spline::inverse_parameter_log_transform( &params[0], 0.0, 5.0, curve2 );
   REQUIRE( curve2.getNumKnots() == curve.getNumKnots() );

   for( std::size_t i = 0; i < curve2.getNumKnots(); ++i )
   {
      REQUIRE( Approx( curve2.getKnot( i ) ) == curve.getKnot( i ) );
   }

}
