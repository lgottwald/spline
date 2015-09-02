#ifndef _SIMPLE_POLYNOMIAL_HPP_
#define _SIMPLE_POLYNOMIAL_HPP_

#include "Polynomial.hpp"

namespace spline
{

template <int DEGREE, typename REAL>
class SimplePolynomial;

template <int DEGREE, typename REAL>
struct DegreeTraits<SimplePolynomial<DEGREE, REAL>>
{
   constexpr static int value = DEGREE;
};

template <int DEGREE, typename REAL>
struct FloatingPointTraits<SimplePolynomial<DEGREE, REAL>>
{
   using type = REAL;
};

/**
 * Simple polynomial class represented by coefficients.
 */
template <int DEGREE, typename REAL>
class SimplePolynomial : public Polynomial< SimplePolynomial<DEGREE, REAL> >
{
public:
   SimplePolynomial() {}
   SimplePolynomial( const REAL coefs[DEGREE + 1] )
   {
      std::copy( coefs, coefs + DEGREE + 1, coeffs );
   }

   void setCoeff( int order, REAL val )
   {
      coeffs[order] = val;

   }

   REAL getCoeff( int order ) const
   {
      return coeffs[order];
   }

   template<typename T>
   T operator()( const T &x ) const
   {
      return evaluatePolynomial<DEGREE>( coeffs, x );
   }

   template<int D, typename T>
   T derivative( const T &x ) const
   {
      return evaluatePolynomial<DEGREE, D>( coeffs, x );
   }


private:
   REAL coeffs[DEGREE + 1];
};

} // namespace spline

#endif
