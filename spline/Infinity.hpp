#ifndef _SPLINE_INFINITY_HPP_
#define _SPLINE_INFINITY_HPP_

#include <limits>

namespace spline
{

/**
 * If plattform supports it return infinity for
 * given floating point type. Else return
 * max value for given floating point type.
 */
template<class REAL>
constexpr REAL infinity()
{
   return std::numeric_limits<REAL>::has_infinity ? std::numeric_limits<REAL>::infinity() : std::numeric_limits<REAL>::max();
}

} //spline


#endif
