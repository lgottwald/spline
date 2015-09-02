#ifndef _SPLINE_INTERVAL_TRAITS_HPP_
#define _SPLINE_INTERVAL_TRAITS_HPP_

#include <type_traits>

namespace spline
{

struct interval_type_undefined {};

/**
 * Traits class to define the type
 * used for representation of an
 * interval of a spline type.
 */
template<typename T>
struct IntervalTypeTraits
{
   using type = interval_type_undefined;
};

template<typename T>
using interval_type = typename IntervalTypeTraits<T>::type;

template<typename T>
constexpr bool interval_type_defined()
{
   return !std::is_same<interval_type<T>, interval_type_undefined>::value;
}

}

#endif
