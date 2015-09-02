#ifndef _FLOATING_POINT_HPP_
#define _FLOATING_POINT_HPP_

#include <type_traits>

namespace spline
{

struct floating_point_undefined {};

/**
 * Traits class to define the floating point type
 * used by the given class.
 */
template<typename T>
struct FloatingPointTraits
{
   using type = floating_point_undefined;
};

template<typename T>
using floating_point = typename FloatingPointTraits<T>::type;

template<typename T>
constexpr bool floating_point_type_defined()
{
   return !std::is_same<floating_point<T>, floating_point_undefined>::value;
}

}

#endif
