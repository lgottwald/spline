#include "catch/catch.hpp"
#include <spline/PiecewisePolynomial.hpp>
#include <vector>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <fstream>

TEST_CASE( "Piecewise linear polynomial test", "[spline]" )
{
   double x[] = {1,   2,   3,   4};
   double y[] = {0.2, 0.0, 0.0, 0.05};
   spline::PiecewisePolynomial<1, double> pp( x, x + 4, y, y + 4 );

   SECTION( "Check if interval search is working" )
   {
      unsigned i = pp.findInterval( 2.8, 0 );
      REQUIRE( pp.getInfimum( i ) == 2.0 );
      REQUIRE( pp.getSupremum( i ) == 3.0 );
      unsigned j = pp.findInterval( 2.8, 4 );
      REQUIRE( i == j );
      REQUIRE( pp.findInterval( 3.0, 3 ) == i + 1 );
   }
   SECTION( "Check if roots are computed correctly" )
   {
      std::vector<double> roots = pp.getRealRoots( -10, 10, 1e-12 );
      REQUIRE( roots.size() == 2 );
      REQUIRE( Approx( roots[0] ) == 2 );
      REQUIRE( Approx( roots[1] ) == 3 );
   }
   SECTION( "Check if evaluated correctly" )
   {
      REQUIRE( pp( -1.0 ) == Approx( 0.2 ) );
      REQUIRE( pp( 1.0 ) == Approx( 0.2 ) );
      REQUIRE( pp( 1.5 ) == Approx( 0.1 ) );
      REQUIRE( pp( 2.0 ) == Approx( 0.0 ) );
      REQUIRE( pp( 3.0 ) == Approx( 0.0 ) );
      REQUIRE( pp( 3.5 ) == Approx( 0.025 ) );
      REQUIRE( pp( 4.0 ) == Approx( 0.05 ) );
   }
}

TEST_CASE( "serialization test", "[spline]" )
{
   double x[] = {1,   2,   3,   4};
   double y[] = {0.2, 0.0, 0.0, 0.05};
   spline::PiecewisePolynomial<1, double> pp( x, x + 4, y, y + 4 );
   {
      std::ofstream ofs( "testspline.dat" );
      boost::archive::binary_oarchive oa( ofs );
      oa << pp;
   }
   spline::PiecewisePolynomial<1, double> pp2;
   {
      std::ifstream ifs( "testspline.dat" );
      boost::archive::binary_iarchive ia( ifs );
      ia >> pp2;
   }

   for( std::size_t i = 0; i < pp2.getBreakpoints().size(); ++i )
      REQUIRE( pp.getBreakpoints()[i] == pp2.getBreakpoints()[i] );

   for( std::size_t i = 0; i < pp2.getCoefficients().size(); ++i )
      REQUIRE( pp.getCoefficients()[i] == pp2.getCoefficients()[i] );
}
