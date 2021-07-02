#include "factorize.hpp"

inline mpz_class f(const mpz_class& x, const mpz_class& c) {
	return x * x + c;
}

Result pollardRhoFloyd(const mpz_class& n, const mpz_class& x0, const mpz_class& c) {
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	mpz_class x1 = x0;
	mpz_class x2 = x0;
	mpz_class diff;
	mpz_class d = 1;
	mpz_class gcdEvaluations = 0;
	mpz_class iterations = 0;
	for (; d == 1; gcdEvaluations++, iterations += 3) {
		x1 = f(x1, c) % n;
		x2 = f(f(x2, c), c) % n;
		diff = x1 - x2;
		mpz_gcd(d.get_mpz_t(), diff.get_mpz_t(), n.get_mpz_t());
	}

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);

	return { d, gcdEvaluations, iterations, elapsed };
}

Result pollardRhoFloydImproved(const mpz_class& n, const mpz_class& x0, const mpz_class& c) {
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	mpz_class x1Save, x2Save;
	mpz_class q;
	mpz_sqrt(q.get_mpz_t(), n.get_mpz_t());
	mpz_sqrt(q.get_mpz_t(), q.get_mpz_t());
	mpz_sqrt(q.get_mpz_t(), q.get_mpz_t());
	mpz_sqrt(q.get_mpz_t(), q.get_mpz_t());

	mpz_class x1 = x0;
	mpz_class x2 = x0;
	mpz_class diff;
	mpz_class d = 1;
	mpz_class gcdEvaluations = 0;
	mpz_class iterations = 0;
	for (; d == 1; gcdEvaluations++) {
		diff = 1;
		for (size_t i = 0; i < q; i++, iterations += 3) {
			x1 = f(x1, c) % n;
			x2 = f(f(x2, c), c) % n;
			diff *= (x1 - x2);
			if (i == 0) {
				x1Save = x1;
				x2Save = x2;
			}
		}
		mpz_gcd(d.get_mpz_t(), diff.get_mpz_t(), n.get_mpz_t());
	}

	if (d == n) {
		x1 = x1Save;
		x2 = x2Save;
		d = 1;
		for (; d == 1; gcdEvaluations++, iterations += 3) {
			x1 = f(x1, c) % n;
			x2 = f(f(x2, c), c) % n;
			diff = x1 - x2;
			mpz_gcd(d.get_mpz_t(), diff.get_mpz_t(), n.get_mpz_t());
		}
	}

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);

	return { d, gcdEvaluations, iterations, elapsed };
}

Result pollardRhoBrent(const mpz_class& n, const mpz_class& x0, const mpz_class& c) {
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	mpz_class powerOfTwo = 1;
	mpz_class xSave;

	mpz_class x1 = f(x0, c) % n;
	mpz_class x2 = f(x1, c) % n;
	mpz_class diff = x1 - x2;
	mpz_class d;
	mpz_gcd(d.get_mpz_t(), diff.get_mpz_t(), n.get_mpz_t());

	mpz_class gcdEvaluations = 1;
	mpz_class iteration = 2;
	for (mpz_class x = x2; d == 1; powerOfTwo *= 2) {
		xSave = x;

		for (; iteration < 3 * powerOfTwo; iteration++) {
			x = f(x, c) % n;
		}

		for (; d == 1 && iteration < powerOfTwo * 4; iteration++, gcdEvaluations++) {
			x = f(x, c) % n;
			diff = xSave - x;
			mpz_gcd(d.get_mpz_t(), diff.get_mpz_t(), n.get_mpz_t());
		}
	}

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);

	return { d, gcdEvaluations, iteration, elapsed };
}

void findSmallFactors(std::vector<mpz_class>& factors, mpz_class& n, const mpz_class& b) {
	std::ifstream primorialFile("C:/Users/psusk/source/repos/C++/PollardRho/b034386.txt"); // https://oeis.org/A034386/b034386.txt
	std::string line;
	size_t i;
	for (i = 0; i < b && std::getline(primorialFile, line); i++);
	if (i != b) {
		std::cout << "Bound to search for small prime numbers is too large.";
		return;
	}

	std::getline(primorialFile, line);
	std::string primorial = line.substr(line.find(' ') + 1);

	mpz_class p, g;

	p.set_str(primorial, 10);

	for (mpz_gcd(g.get_mpz_t(), p.get_mpz_t(), n.get_mpz_t()); g > 1; mpz_gcd(g.get_mpz_t(), p.get_mpz_t(), n.get_mpz_t())) {
		std::vector<int> basicDivisors = { 2, 3, 5 };
		std::vector<int> divisors = { 1, 7, 11, 13, 17, 19, 23, 29 };

		for (const int& divisor : basicDivisors) {
			if (g % divisor == 0) {
				factors.push_back(divisor);
			}
		}

		for (mpz_class expression = 7, k = 0; expression < b; k++) {
			for (size_t i = k == 0 ? 1 : 0; i < 8 && expression < b; i++) {
				expression = 30 * k + divisors[i];
				if (g % expression == 0) {
					factors.push_back(expression);
				}
			}
		}

		n /= g;
	}
}

void findLargeFactors(std::vector<mpz_class>& factors, mpz_class n, const mpz_class& b, const mpz_class& x0, const mpz_class& c) {
	// Check if n is prime
	int isPrime = mpz_probab_prime_p(n.get_mpz_t(), 10); // 2 denotes guaranteed prime, 1 denotes probably prime, 0 denotes guaranteed composite
	if (isPrime == 2) {
		factors.push_back(n);
		return;
	}
	else {
		// Check if n is a perfect power
		mpfr::mpreal max_k = floor(log(mpfr::mpreal(n.get_str())) / log(mpfr::mpreal(b.get_str()))).toULLong();
		for (size_t k = max_k.toULLong(); k >= 2; k--) {
			mpz_class result;
			mpz_class remainder;
			mpz_rootrem(result.get_mpz_t(), remainder.get_mpz_t(), n.get_mpz_t(), k);
			if (remainder == 0) {
				std::vector<mpz_class> perfectPowerFactors(k, result);
				factors.insert(factors.end(), perfectPowerFactors.begin(), perfectPowerFactors.end());
				return;
			}
		}

		// Apply Pollard rho
		mpz_class factor;
		if (isPrime == 1) {
			size_t run = 0;
			size_t runs = 10;
			//for (factor = pollardRhoFloyd(n, x0, c); run < runs && factor == n; x0 = (x0 + 1) % n, factor = pollardRhoFloyd(n, x0, c), run++);
			if (factor == n) {
				factors.push_back(factor);
				return;
			}
		}
		else {
			//for (factor = pollardRhoFloyd(n, x0, c); factor == n; x0 = (x0 + 1) % n, factor = pollardRhoFloyd(n, x0, c));
		}
		findLargeFactors(factors, factor, b, x0, c);
		findLargeFactors(factors, n / factor, b, x0, c);
	}
}

std::vector<mpz_class> findFactors(mpz_class n, mpz_class b, mpz_class x0, mpz_class c) {
	std::vector<mpz_class> factors;

	// Find factors below bound b
	findSmallFactors(factors, n, b);

	// Find remaining large factors
	findLargeFactors(factors, n, b, x0, c);

	return factors;
}