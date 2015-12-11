#ifndef _UNIVARIATE_SPLINE_ESTIMATOR_HPP_
#define _UNIVARIATE_SPLINE_ESTIMATOR_HPP_

#include "UnivariateSplineCurve.hpp"
#include "SimplePolynomial.hpp"
#include "SplineIntervalTraits.hpp"
#include "FloatingPointTraits.hpp"
#include "DegreeTraits.hpp"
#include <map>
#include <vector>
#include <limits>

namespace spline
{
template < typename Impl >
class UnivariateSplineEstimator;

template < typename Impl >
struct DegreeTraits< UnivariateSplineEstimator<Impl> >
{
   constexpr static int value = get_degree<Impl>();
};

template < typename Impl >
struct FloatingPointTraits< UnivariateSplineEstimator<Impl> >
{
   using type = floating_point<Impl>;
};

template < typename Impl >
struct IntervalTypeTraits< UnivariateSplineEstimator<Impl> >
{
   struct bitangent
   {
      floating_point<UnivariateSplineEstimator<Impl>> start;
      SimplePolynomial<1, floating_point<UnivariateSplineEstimator<Impl>> > line;
   };

   using type = typename std::map<floating_point<UnivariateSplineEstimator<Impl>> , bitangent>::const_iterator;
};

/**
 * This class represents a convex underestimator or a concave overestimator
 * and is created from an interval and an UnivariateSplineCurve.
 */
template < typename Impl >
class UnivariateSplineEstimator : public UnivariateSplineCurve< UnivariateSplineEstimator<Impl> >
{
public:
   using real_t = floating_point< UnivariateSplineEstimator<Impl> >;
   using interval_t = interval_type< UnivariateSplineEstimator<Impl> >;
   using bitangent = typename IntervalTypeTraits< UnivariateSplineEstimator<Impl> >::bitangent;

   UnivariateSplineEstimator( const UnivariateSplineCurve<Impl> &curve, real_t l, real_t r, bool overestimate, real_t epsilon )
      : curve( curve ), a( l ), b( l )
   {
      enum curvature_t
      {
         NIL = 0,
         CONVEX = 1,
         CONCAVE = 2,
         LINEAR = CONVEX | CONCAVE
      };

      auto derivative2 = differentiate<2>( curve );

      curvature_t curvature = NIL;
      real_t val = derivative2( a );

      if( val >= 0.0 )
         curvature = curvature_t( curvature | CONVEX );

      if( val <= 0.0 )
         curvature = curvature_t( curvature | CONCAVE );

      bool curvature_ok = overestimate ? curvature & CONCAVE : curvature & CONVEX;

      real_t sweepline = a;

      std::vector<real_t> eventpoint_candidates = derivative2.getRealRoots( l, r, epsilon );

      for( auto iter = eventpoint_candidates.begin();
            iter != eventpoint_candidates.end();
            ++iter )
      {
         curvature_t c = NIL;
         val = derivative2( *iter + epsilon );

         if( val >= 0.0 )
            c = curvature_t( c | CONVEX );;

         if( val <= 0.0 )
            c = curvature_t( c | CONCAVE );;

         bool eventpoint = overestimate ? // check if we got an eventpoint
                           ( curvature & CONCAVE ) ^ ( c & CONCAVE ) : // if overestimating and concavity changed we got one
                           ( curvature & CONVEX ) ^ ( c & CONVEX ); // and if underestimating and convexity changed too

         if( eventpoint )
         {
            if( curvature_ok ) // if the curvature was ok till here we need to update the estimator
            {
               /* if the sweepline moved past the current endpoint of the estimator
                * we need to add a bitangent since it means there was a region where
                * the curvature was not ok but then we entered a region were the curvature
                * was ok again and we are leaving it now. The bitangent will connect
                * the region that we are leaving now with the current estimator.
                * If the sweepline is at b or before it the region we are leaving
                * can safely be added to the current estimator.
                */
               if( sweepline > b )
                  addBitangent( sweepline, *iter, overestimate, epsilon );

               // set estimator region to current point
               b = *iter;
            }

            sweepline = *iter; // sweepline is set to current eventpoint
            curvature_ok = !curvature_ok; //curvature is negated
         }

         // update curvature state
         curvature = c;
      }

      if( sweepline > b )
      {
         //connect region starting at sweepline and ending in the right enpoint
         //to the estimator
         addBitangent( sweepline, r, overestimate, epsilon );
      }
      else if( overestimate ? !( curvature & CONCAVE ) : !( curvature & CONVEX ) )
      {
         //connect right endpoint to estimator
         addBitangent( r, r, overestimate, epsilon );
      }

      b = r;
   }

   void addBitangent(
      const real_t l,
      const real_t r,
      bool overestimate,
      const real_t epsilon
   )
   {
      bitangent bitan;
      bitan.start = a;
      bool start_fixed = b - a <= epsilon;
      real_t end = r;
      bool end_fixed = r - l  <= epsilon;

      if( start_fixed && end_fixed ) // if both fixed the secant will be used
      {

         real_t y1 = ( *this )( bitan.start );
         real_t y2 = curve( end );
         real_t slope = ( y2 - y1 ) / ( end - bitan.start );

         //this case happens when the function is concave and should be underestimated
         //or it is convex and should be overestimated. Then the secant of the endpoints is
         //used as estimator
         bitan.line.setCoeff( 1, slope );
         bitan.line.setCoeff( 0, y2 - slope * r );
         addBitangent( end, bitan );
         return;
      }

      auto curveDeriv = differentiate<1>( curve );
      SimplePolynomial< 0, real_t > rhs;

      real_t y1 = ( *this )( bitan.start );
      real_t y2 = curve( end );

      //compute slope of line connecting the function values
      real_t slope = ( y2 - y1 ) / ( end - bitan.start );

      while( 1 )
      {
         //check for termination, i.e. slope matches curve slope at both parts
         //or one of the parts is fixed
         if( ( start_fixed || std::abs( curveDeriv( bitan.start ) - slope ) <= epsilon ) &&
               ( end_fixed || std::abs( curveDeriv( end ) - slope ) <= epsilon ) )
         {
            bitan.line.setCoeff( 1, slope ); // use current slope
            bitan.line.setCoeff( 0, y2 - end * slope ); // compute constant value using  y2 = end*slope+c <=> c = y2 - end*slope
            //add tangent
            addBitangent( end, bitan );
            return;
         }

         //if not terminating compute new endpoints for tangent
         rhs.setCoeff( 0, slope );

         if( !start_fixed )
         {
            bool changed = false;
            std::vector<real_t> sol = curveDeriv.solveEquation( rhs, a, b, epsilon );
            bool nosol = true;

            //filter out solutions that are cut off by current tangents
            //after that there should be 1 solution or 0.
            for( auto i = sol.begin(); i != sol.end(); ++i )
            {
               interval_t interv = this->findInterval( *i );

               if( interv == bitangents.end() )
               {
                  nosol = false;

                  if( bitan.start != *i )
                  {
                     bitan.start = *i;
                     changed = true;
                  }

                  break;
               }
            }

            if( nosol ) // if there was no solution, i.e. all solutions were cut off by current tagents we choose an interval endpoint
            {
               real_t newstart;

               if( curve( a ) - a * slope > curve( b ) - b * slope )
                  newstart = overestimate ? a : b;
               else
                  newstart = overestimate ? b : a;

               if( newstart != bitan.start )
               {
                  changed = true;
                  bitan.start = newstart;
               }
            }

            if( changed )
            {
               y1 = ( *this )( bitan.start );
               slope = ( y2 - y1 ) / ( end - bitan.start );
            }
            else
            {
               start_fixed = true;
            }
         }

         if( !end_fixed )
         {
            bool changed = false;
            //should always have size 1 or 0 since curve is convex/concave in [l,r]
            std::vector<real_t> sol = curveDeriv.solveEquation( rhs, l, r, epsilon );

            if( sol.empty() )
            {
               real_t newend;

               if( curve( l ) - slope * l > curve( r ) - slope * r )
                  newend = overestimate ? l : r;
               else
                  newend = overestimate ? r : l;

               if( newend != end )
               {
                  end = newend;
                  changed = true;
               }
            }
            else
            {
               if( end != sol[0] )
               {
                  end = sol[0];
                  changed = true;
               }
            }

            if( changed )
            {
               y2 = curve( end );
               slope = ( y2 - y1 ) / ( end - bitan.start );
            }
            else
            {
               end_fixed = true;
            }
         }
      }
   }

   real_t getInfimumImpl( const interval_t interval ) const
   {
      return interval == bitangents.end() ? a : interval->second.start;
   }

   real_t getSupremumImpl( const interval_t interval ) const
   {
      return interval == bitangents.end() ? b : interval->first;
   }

   real_t getInfimumImpl() const
   {
      return curve.getInfimum();
   }

   real_t getSupremumImpl() const
   {
      return curve.getSupremum();
   }

   real_t getStart() const
   {
      return a;
   }

   real_t getEnd() const
   {
      return b;
   }

   std::size_t numIntervalsImpl() const
   {
      return bitangents.size() + 1;
   }

   template < typename T >
   interval_t findIntervalImpl( const T &x, const interval_t hint ) const
   {
      interval_t b = bitangents.upper_bound( x );

      if( b == bitangents.end() )
         return b;

      if( b->second.start <= x )
         return b;

      return bitangents.end();
   }

   template < int D, typename T >
   T evaluate( const T &x, const interval_t interval ) const
   {
      if( interval == bitangents.end() )
         return curve.template evaluate<D>( x, curve.findInterval( x ) );

      return interval->second.line.template derivative<D>( x );
   }

private:

   void addBitangent( real_t end, const bitangent &bitan )
   {
      bitangents.erase( bitangents.upper_bound( bitan.start ), bitangents.end() );
      bitangents[ end ] = bitan;
   }


   //store the curve and the bitangents sorted by their endpoints
   std::map<real_t, bitangent> bitangents;
   const UnivariateSplineCurve<Impl> &curve;
   //interval end points the curve is over/under estimated in
   real_t a;
   real_t b;

};


} // namespace spline

#endif
