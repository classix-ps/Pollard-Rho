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

struct Result {
	mpz_class value;
	size_t gcdEvaluations;
	std::chrono::microseconds elapsed;
};

mpz_class f(mpz_class x, mpz_class c);

Result pollardRhoFloyd(mpz_class n, mpz_class x0, mpz_class c);
Result pollardRhoFloydImproved(mpz_class n, mpz_class x0, mpz_class c);

void findSmallFactors(std::vector<mpz_class>& factors, mpz_class& n, mpz_class b);

void findLargeFactors(std::vector<mpz_class>& factors, mpz_class n, mpz_class b, mpz_class x0, mpz_class c);

std::vector<mpz_class> findFactors(mpz_class n, mpz_class b, mpz_class x0, mpz_class c);