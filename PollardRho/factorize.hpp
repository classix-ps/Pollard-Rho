#pragma once

#include <vector>
#include <iostream>
#include <numeric>
#include <chrono>
#include <map>

#include <gmpxx.h> // mpz_class https://gmplib.org/manual/C_002b_002b-Interface-General
#include <mpreal.h> // mpfr::mpreal https://github.com/advanpix/mpreal
//#include <mpir.h> // wrapped by mpz_class for C++

#include <gnuplot-iostream.h>
#include <fstream>

#include <boost/parameter/keyword.hpp>
#include <boost/parameter/name.hpp>
#include <boost/parameter/preprocessor.hpp>

struct Result {
	mpz_class value;
	mpz_class gcdEvaluations;
	mpz_class iterations;
	std::chrono::microseconds elapsed;
};

enum class PollardRho { Floyd, FloydImproved, Brent };

inline mpz_class f(const mpz_class& x, const mpz_class& c) { return x * x + c; }

Result pollardRhoFloyd(const mpz_class& n, const mpz_class& x0, const mpz_class& c);
Result pollardRhoFloydImproved(const mpz_class& n, const mpz_class& x0, const mpz_class& c);
Result pollardRhoBrent(const mpz_class& n, const mpz_class& x0, const mpz_class& c);

Result pollardRhoOne(const mpz_class& n, const mpz_class& s);

namespace Factorize {
	void _removeSmallFactors(std::vector<mpz_class>& factors, mpz_class& n, const mpz_class& b);

	void _getAllFactors(std::vector<mpz_class>& factors, mpz_class n, const mpz_class& b, const mpz_class& x0, const mpz_class& c);

	namespace keywords {
		BOOST_PARAMETER_NAME(n)
		BOOST_PARAMETER_NAME(b)
		BOOST_PARAMETER_NAME(s)
		BOOST_PARAMETER_NAME(x0)
		BOOST_PARAMETER_NAME(c)
		BOOST_PARAMETER_NAME(pRho)
	}

	/** This function returns a vector of all prime factors n, or n itself if n is prime.
	 * @param n (REQUIRED) The number to factor
	 * @param b (OPTIONAL) The bound for small factors to search for
	 * @param s (OPTIONAL) The bound for which to find factors - 1 which are s-powersmooth
	 * @param x0 (OPTIONAL) The initial value for Pollard Rho factoring
	 * @param c (OPTIONAL) The constant for Pollard Rho factoring (Avoid 0, -2, and c = k * N - x0 * (x0 +- 1) where k is an integer)
	 * @param pRho (OPTIONAL) The type of Pollard Rho algorithm to use when applying Pollard Rho
	 */
	BOOST_PARAMETER_FUNCTION(
		(std::vector<mpz_class>),
		findFactors,
		keywords::tag,
		(required
			(n, (mpz_class))
		)
		(optional
			(b, (mpz_class), 1699)

			(s, (mpz_class), 2000)

			(x0, (mpz_class), 2)
			(c, (mpz_class), 1)
			(pRho, (PollardRho), PollardRho::FloydImproved)
		)
	)
	{
		std::vector<mpz_class> factors;

		// Find factors below bound b
		_removeSmallFactors(factors, n, b);

		if (n == 1)
			return factors;

		// Find factor - 1 which is s-powersmooth
		mpz_class factor = pollardRhoOne(n, s).value;

		if (factor != n && factor != 1) {
			_getAllFactors(factors, factor, b, x0, c);
			n /= factor;
		}

		// Find all remaining factors
		_getAllFactors(factors, n, b, x0, c);

		return factors;
	}
}