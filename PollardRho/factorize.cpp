#include "factorize.hpp"

Result pollardRhoFloyd(const mpz_class& n, const mpz_class& x0, const mpz_class& c) {
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	mpz_class x1 = x0;
	mpz_class x2 = x0;
	mpz_class diff;
	mpz_class d = 1;
	mpz_class gcdEvaluations = 0;
	mpz_class iteration = 0;
	for (; d == 1; gcdEvaluations++, iteration += 3) {
		x1 = f(x1, c) % n;
		x2 = f(f(x2, c), c) % n;
		diff = x1 - x2;
		mpz_gcd(d.get_mpz_t(), diff.get_mpz_t(), n.get_mpz_t());
	}

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);

	return { d, gcdEvaluations, iteration, elapsed };
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
	mpz_class iteration = 0;
	for (; d == 1; gcdEvaluations++) {
		diff = 1;
		for (size_t i = 0; i < q; i++, iteration += 3) {
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
		for (; d == 1; gcdEvaluations++, iteration += 3) {
			x1 = f(x1, c) % n;
			x2 = f(f(x2, c), c) % n;
			diff = x1 - x2;
			mpz_gcd(d.get_mpz_t(), diff.get_mpz_t(), n.get_mpz_t());
		}
	}

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);

	return { d, gcdEvaluations, iteration, elapsed };
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

Result pollardRhoBrentImproved(const mpz_class& n, const mpz_class& x0, const mpz_class& c) {
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	std::pair<std::pair<mpz_class, mpz_class>, std::pair<mpz_class, mpz_class>> save;
	mpz_class q;
	mpz_sqrt(q.get_mpz_t(), n.get_mpz_t());
	mpz_sqrt(q.get_mpz_t(), q.get_mpz_t());
	mpz_sqrt(q.get_mpz_t(), q.get_mpz_t());
	mpz_sqrt(q.get_mpz_t(), q.get_mpz_t());

	mpz_class powerOfTwo = 1;
	mpz_class xSave;

	mpz_class x1 = f(x0, c) % n;
	mpz_class x2 = f(x1, c) % n;
	mpz_class diff = x1 - x2;
	mpz_class d;
	mpz_gcd(d.get_mpz_t(), diff.get_mpz_t(), n.get_mpz_t());

	mpz_class i = 0;

	mpz_class gcdEvaluations = 1;
	mpz_class iteration = 2;
	for (mpz_class x = x2; d == 1; powerOfTwo *= 2) {
		xSave = x;

		for (; iteration < 3 * powerOfTwo; iteration++) {
			x = f(x, c) % n;
		}

		diff = 1;
		for (; d == 1 && iteration < powerOfTwo * 4; i++, iteration++) {
			if (i == q) {
				mpz_gcd(d.get_mpz_t(), diff.get_mpz_t(), n.get_mpz_t());
				gcdEvaluations++;

				save = std::make_pair(std::make_pair(powerOfTwo, xSave), std::make_pair(iteration, x));
				i = 0;
			}

			x = f(x, c) % n;
			diff *= xSave - x;
		}
	}

	mpz_class iterations = iteration;
	if (d == n) {
		powerOfTwo = save.first.first;
		xSave = save.first.second;
		iteration = save.second.first;
		for (mpz_class x = save.second.second; d == 1; powerOfTwo *= 2) {
			for (; iteration < 3 * powerOfTwo; iteration++, iterations++) {
				x = f(x, c) % n;
			}

			for (; d == 1 && iteration < powerOfTwo * 4; iteration++, iterations++, gcdEvaluations++) {
				x = f(x, c) % n;
				diff = xSave - x;
				mpz_gcd(d.get_mpz_t(), diff.get_mpz_t(), n.get_mpz_t());
			}
		}
	}

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);

	return { d, gcdEvaluations, iteration, elapsed };
}

Result pollardPOne(const mpz_class& n, const mpz_class& s) {
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	mpz_class g;
	if (mpz_odd_p(n.get_mpz_t())) {
		g = 2;
	}
	else {
		for (g = 3; n % g != 0; g += 2);
	}

	mpz_class r = g;
	mpz_class iteration = 0;
	for (mpz_class i = 2; i <= s; mpz_nextprime(i.get_mpz_t(), i.get_mpz_t()), iteration++) {
		mpfr::mpreal alpha = floor(log(mpfr::mpreal(s.get_str())) / log(mpfr::mpreal(i.get_str())));
		mpz_class q;
		mpz_pow_ui(q.get_mpz_t(), i.get_mpz_t(), alpha.toULong());

		mpz_powm(r.get_mpz_t(), r.get_mpz_t(), q.get_mpz_t(), n.get_mpz_t());
	}

	r -= 1;
	mpz_class d;
	mpz_gcd(d.get_mpz_t(), r.get_mpz_t(), n.get_mpz_t());

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);

	return { d, 1, iteration, elapsed };
}

void Factorize::_removeSmallFactors(std::vector<mpz_class>& factors, mpz_class& n, const mpz_class& b) {
	mpz_class p = primorial(b); // Alternatively, can access OEIS bFile containing primorials up to 2000
	mpz_class g;

	std::vector<int> basicDivisors = { 2, 3, 5 };
	std::vector<int> divisors = { 1, 7, 11, 13, 17, 19, 23, 29 };
	for (mpz_gcd(g.get_mpz_t(), p.get_mpz_t(), n.get_mpz_t()); g > 1; mpz_gcd(g.get_mpz_t(), p.get_mpz_t(), n.get_mpz_t())) {
		n /= g;

		for (const int& divisor : basicDivisors) {
			if (g % divisor == 0) {
				factors.push_back(divisor);
				g /= divisor;
			}
		}

		for (mpz_class expression = 7, k = 0; expression < b; k++) {
			for (size_t i = k == 0 ? 1 : 0; i < 8 && expression < b; i++) {
				expression = 30 * k + divisors[i];
				if (g % expression == 0) {
					factors.push_back(expression);
					g /= expression;
				}
			}
		}
	}
}

void Factorize::_getAllFactors(std::vector<mpz_class>& factors, mpz_class n, const mpz_class& b, const mpz_class& x0, const mpz_class& c, PollardRho pRho) {
	// Check if n is prime
	int isPrime = mpz_probab_prime_p(n.get_mpz_t(), 10); // 2 denotes guaranteed prime, 1 denotes probably prime, 0 denotes guaranteed composite
	if (isPrime == 2) {
		factors.push_back(n);
		return;
	}
	else {
		// Check if n is a perfect power (alternatively, mpz_perfect_power_p(), but we want to use the fact that there are no factors below b left in N to our advantage)
		mpfr::mpreal max_k = b > 1 ? floor(log(mpfr::mpreal(n.get_str())) / log(mpfr::mpreal(b.get_str()))) : floor(log(mpfr::mpreal(n.get_str())) / log(mpfr::mpreal(2)));
		for (mpz_class k = 2; k <= max_k.toULong(); mpz_nextprime(k.get_mpz_t(), k.get_mpz_t())) {
			mpz_class result, remainder;
			mpz_rootrem(result.get_mpz_t(), remainder.get_mpz_t(), n.get_mpz_t(), k.get_ui());
			if (remainder == 0) {
				std::vector<mpz_class> perfectPowerFactors(k.get_ui(), result);
				factors.insert(factors.end(), perfectPowerFactors.begin(), perfectPowerFactors.end());
				return;
			}
		}

		// Apply factoring algorithm
		Result factor;
		mpz_class cInc(c);
		std::function<Result(const mpz_class&, const mpz_class&, const mpz_class&)> rhoFunction;
		switch (pRho) {
		case PollardRho::Floyd:
			rhoFunction = pollardRhoFloyd;
			break;
		case PollardRho::FloydImproved:
			rhoFunction = pollardRhoFloydImproved;
			break;
		case PollardRho::Brent:
			rhoFunction = pollardRhoBrent;
			break;
		}
		// This is the most dangerous spot to be in: we aren't sure if N is prime, so we are forced to try Pollard Rho (multiple times unless we find a factor, since the algorithm can fail for some (consecutive) c)
		// We thereby must run the algorithm arbitrarily many times and, upon failing every time, assume that N is prime. One must note that the runtime for failure for prime N is O(sqrt(N)) rather than O(sqrt(p)) (failure due to c is still O(sqrt(p))
		if (isPrime == 1) {
			size_t runs = 5;
			for (factor = rhoFunction(n, x0, cInc); runs && factor.value == n; factor = rhoFunction(n, x0, cInc), runs--) {
				// 3 cases of c we want to avoid: 0 or -2 mod N, and staying at x0. (last occurs when c = - x0 * (x0 +- 1) mod N)
				// We want to increment c because incrementing x0 is less likely to avoid consecutive failures
				for (cInc = (cInc + 1) % n; cInc % n == 0 || (cInc + 2) % n == 0 || (cInc + x0 * (x0 + 1)) % n == 0 || (cInc + x0 * (x0 - 1)) % n == 0; cInc = (cInc + 1) % n);
			}
			if (factor.value == n) {
				factors.push_back(factor.value);
				return;
			}
		}
		// Here we know N is composite, so we apply Pollard Rho until we find a factor
		else {
			for (factor = rhoFunction(n, x0, cInc); factor.value == n; factor = rhoFunction(n, x0, cInc)) {
				// 3 cases of c we want to avoid: 0 or -2 mod N, and staying at x0. (last occurs when c = - x0 * (x0 +- 1) mod N),
				// We want to increment c because incrementing x0 is less likely to avoid consecutive failures
				for (cInc = (cInc + 1) % n; cInc % n == 0 || (cInc + 2) % n == 0 || (cInc + x0 * (x0 + 1)) % n == 0 || (cInc + x0 * (x0 - 1)) % n == 0; cInc = (cInc + 1) % n);
			}
		}
		Factorize::_getAllFactors(factors, factor.value, b, x0, c, pRho);
		Factorize::_getAllFactors(factors, n / factor.value, b, x0, c, pRho);
	}
}