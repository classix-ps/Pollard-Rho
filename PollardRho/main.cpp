#include "test.hpp"

void displayBasics() {
	Result factor;

	std::cout << "Primes are allowed, but they'll take longer than normal: O(sqrt(n)) instead of O(sqrt(p)), and leave you uncertain if the algorithm failed." << std::endl;
	factor = pollardRhoFloyd(23, 2, 1);
	std::cout << "Factoring 23: " << factor.value.get_str() << ", gcd evaluations: " << factor.gcdEvaluations << std::endl << std::endl;

	std::cout << "Perfect powers are allowed but unfactorizable, all the factors are the same, so they all collide at the same time, meaning N is always found." << std::endl;
	factor = pollardRhoFloyd(8, 2, 1);
	std::cout << "Factoring 8: " << factor.value.get_str() << ", gcd evaluations: " << factor.gcdEvaluations << std::endl << std::endl;

	std::cout << "This obviously still holds when perfect powers are multiplied with additional different factors, you'll never be able to find the root of the perfect power as a factor,"
		"yet the likelihood of finding the perfect power is governed by its root" << std::endl << std::endl;
	factor = pollardRhoFloyd(8 * 3, 2, 1);
	std::cout << "Factoring 8 * 3, x0 = 2, c = 1: " << factor.value.get_str() << ", gcd evaluations: " << factor.gcdEvaluations << std::endl << std::endl; // here we find 3
	factor = pollardRhoFloyd(8 * 3, 2, 2);
	std::cout << "Factoring 8 * 3, x0 = 2, c = 2: " << factor.value.get_str() << ", gcd evaluations: " << factor.gcdEvaluations << std::endl << std::endl; // here we find 8
	factor = pollardRhoFloyd(8 * 3, 2, 3);
	std::cout << "Factoring 8 * 3, x0 = 2, c = 3: " << factor.value.get_str() << ", gcd evaluations: " << factor.gcdEvaluations << std::endl << std::endl; // here we find 8
}

int main() {
	/*
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	//std::vector<mpz_class> factors = findFactors(509747, 30, 30, 1);
	//std::vector<mpz_class> factors = findFactors(10403, 102, 30, 1);
	//std::vector<mpz_class> factors = findFactors(49, 5, 30, 1);
	//std::vector<mpz_class> factors = findFactors(4307, 5, 2, 1);
	std::vector<mpz_class> factors = findFactors(mpz_class("5915587277") * mpz_class("3267000013"), 5, 2, 1);
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
	std::cout << "Elapsed time: " << elapsed.count() << " [microseconds]" << std::endl;
	std::cout << "Factors are:" << std::endl;
	for (const mpz_class& factor : factors) {
		std::string fStr = factor.get_str();
		std::cout << "\t" << factor.get_str() << std::endl;
	}
	*/

	//findSameMod(59, 73); // 4307
	//findSameMod(4, 6); // 24
	//testAverageP(2, 9, 50, 100);
	//testAverageP(2, 9, 50, 1000);
	//testAverageP(3, 3, 50, 10000);
	//testAverageC(3, 9, 50, 10000);
	//testAverageC(3, 3, 100, 50000);
	//visualizeIterationsC(3, 2197);
	//visualizeIterationsC(3, 2000);
	//visualizeIterationsC(3, 2209);
	//visualizeIterationsP(-5, 3);
	//visualizeIterationsP(1, 3);
	//visualizeIterationsP(-20, 3);
	
	//displayBasics();

	Result factor;
	factor = pollardRhoFloyd(mpz_class("5915587277") * mpz_class("3267000013"), 2, 1);
	std::cout << factor.value.get_str() << ", gcd evaluations: " << factor.gcdEvaluations << ", elapsed: " << factor.elapsed.count() << std::endl;
	factor = pollardRhoFloydImproved(mpz_class("5915587277") * mpz_class("3267000013"), 2, 1);
	std::cout << factor.value.get_str() << ", gcd evaluations: " << factor.gcdEvaluations << ", elapsed: " << factor.elapsed.count() << std::endl;

	//testPollardFactor(59 * 73);
	//testPollardFactor(3 * 59 * 73);
	//testPollardFactor(13 * 59 * 73);

	//std::cout << findCollision(2209, 2, -6) << std::endl; // this is stupid
	//std::cout << findCollision(2209, 2, 2203) << std::endl; // this is lucky (don't do this for N though! only p)
	
	//std::cout << findCollision(2209, 2, -2) << std::endl; // this is stupid
	//std::cout << findCollision(2209, 2, 2207) << std::endl; // this is lucky (don't do this for N though! only p)

	//std::cout << pollardRhoFloyd(2209 * 5, 2, 2207).get_str() << std::endl; // here's being lucky in action: 2209 is always found in 1 iteration
	//std::cout << pollardRhoFloyd(2209 * 5, 2, -2).get_str() << std::endl; // here's begin stupid in action: 2209 * 5 is always found in 1 iteration
}