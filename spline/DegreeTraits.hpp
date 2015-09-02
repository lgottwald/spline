#ifndef _DEGREE_TRAITS_HPP_
#define _DEGREE_TRAITS_HPP_

#include <type_traits>
namespace spline
{

struct degree_undefined {};

/**
 * Trait class to define the degree of
 * of spline and polynomial classes.
 */
template<typename T>
struct DegreeTraits
{
   constexpr static degree_undefined value {};
};

template<typename T>
constexpr int get_degree()
{
   return DegreeTraits<T>::value;
}

template<typename T>
constexpr bool has_degree()
{
   return !std::is_same<decltype( DegreeTraits<T>::value ), degree_undefined>::value;
}

}

#endif
