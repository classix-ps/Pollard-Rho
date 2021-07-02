#include "test.hpp"

struct Compare {
	int val;
	Compare(const int& n) : val(n) {}
};
bool operator==(const std::pair<int, int>& p, const Compare& c) {
	return c.val == p.second;
}
bool operator==(const Compare& c, const std::pair<int, int>& p) {
	return c.val == p.second;
}

void findSameMod(int p, int q) {
	int n = p * q;
	if (n == 0) {
		std::cout << "Can't find overlapping mods when the product of the factors is 0." << std::endl;
		return;
	}

	std::map<int, std::map<int, int>> mods;
	for (int i = 0; i < n; i += boost::math::sign(n)) {
		int pMod = i % p;
		int qMod = i % q;
		if (mods.find(pMod) != mods.end() && mods[pMod].find(qMod) != mods[pMod].end()) {
			std::cout << mods[pMod][qMod] << ", " << i << std::endl;
		}
		else {
			mods[pMod][qMod] = i;
		}
	}

	return;
}

void testPollardFactor(int n) {
	std::map<int, std::pair<int, double>> factors;

	int x0 = 2; // We have shown that changing the inital value has negligible impact on iterations
	for (int c = 1; c < n; c++) {
		// Staying at x0 causes failure in N and occurs when c = k * N - x0 * (x0 +- 1), where k is an integer
		if (c == n - x0 * (x0 + 1) || c == n - x0 * (x0 - 1)) {
			continue;
		}

		Result factor = pollardRhoFloyd(n, x0, c);
		if (factor.value.fits_sint_p()) { // Should alwys be true, since n is an int and its factors must be smaller, but just to be sure
			factors[factor.value.get_si()].first++;
			factors[factor.value.get_si()].second += (factor.gcdEvaluations.get_si() - factors[factor.value.get_si()].second) / factors[factor.value.get_si()].first;
		}
	}

	std::cout << "Occurences of factors found by Pollard Rho:" << std::endl;
	for (const std::pair<int, std::pair<int, double>>& factor : factors) {
		std::cout << "\t" << factor.first << ", relative frequency: " << double(factor.second.first) / (n - 1) << ", mean gcd evaluations: " << factor.second.second << std::endl;
	}
}

void testPollardRhoRuntime(int maxC, int maxN) {
	std::vector<std::pair<unsigned long, double>> averageGCDEvaluations;
	std::vector<std::pair<unsigned long, double>> averageIterations;
	std::vector<std::pair<unsigned long, double>> fourthRoot;

	std::ifstream file("../b001097 (twin primes).txt");
	std::string line1, line2;
	
	while (std::getline(file, line1)) {
		if (std::getline(file, line2)) {
			mpz_class p1(line1.substr(line1.find(' ') + 1));
			mpz_class p2(line2.substr(line2.find(' ') + 1));
			mpz_class n = p1 * p2;
			if (n > maxN) {
				break;
			}

			std::vector<unsigned long> gcdEvaluations;
			std::vector<unsigned long> iterations;
			for (int c = 1; c < maxC; c++) {
				Result res = pollardRhoFloyd(n, 2, 1);
				gcdEvaluations.push_back(res.gcdEvaluations.get_ui());
				iterations.push_back(res.iterations.get_ui());
			}
			averageGCDEvaluations.push_back(std::make_pair(n.get_ui(), std::accumulate(gcdEvaluations.begin(), gcdEvaluations.end(), 0.0) / gcdEvaluations.size()));
			averageIterations.push_back(std::make_pair(n.get_ui(), std::accumulate(iterations.begin(), iterations.end(), 0.0) / iterations.size()));
			fourthRoot.push_back(std::make_pair(n.get_ui(), sqrt(sqrt(n.get_ui()))));
		}
	}

	Gnuplot gp("gnuplot -persist");
	double maxG = std::max_element(averageGCDEvaluations.begin(), averageGCDEvaluations.end(), [](const std::pair<unsigned long, double>& p1, const std::pair<unsigned long, double>& p2) { return p1.second < p2.second; })->second;
	double maxI = std::max_element(averageIterations.begin(), averageIterations.end(), [](const std::pair<unsigned long, double>& p1, const std::pair<unsigned long, double>& p2) { return p1.second < p2.second; })->second;
	double max4 = std::max_element(fourthRoot.begin(), fourthRoot.end(), [](const std::pair<unsigned long, double>& p1, const std::pair<unsigned long, double>& p2) { return p1.second < p2.second; })->second;
	std::vector<double> maxVals = { maxG, maxI, max4 };

	gp << "set xrange [1:" << maxN << "]\n";
	gp << "set yrange [0:" << *std::max_element(maxVals.begin(), maxVals.end()) << "]\n";
	gp << "plot";
	gp << gp.file1d(averageGCDEvaluations) << "with lines title 'gcdEvaluations',";
	gp << gp.file1d(averageIterations) << "with lines title 'iterations',";
	gp << gp.file1d(fourthRoot) << "with lines title 'linear'";
	gp << std::endl;

	file.close();
}

int findCollision(mpz_class n, mpz_class x0, mpz_class c) {
	c += n * (c < 0); // make c positive mod N to avoid having to check if x1 or x2 are negative after modulo

	mpz_class x1 = f(x0, c) % n;
	mpz_class x2 = f(x1, c) % n;
	int iterations = 1;
	for (; x1 != x2; iterations++) {
		x1 = f(x1, c) % n;
		x2 = f(f(x2, c), c) % n;
	}

	return iterations;
}

void testAverageP(int minX0, int maxX0, int minP, int maxP) {
	std::map<int, std::map<int, double>> averageIterations; // x0, p, average across c
	std::vector<std::pair<int, double>> sqrtFunc;
	std::ofstream file("../notableIterations.txt");

	for (int x0 = minX0; x0 <= maxX0; x0++) {
		for (int p = minP; p <= maxP; p++) {
			sqrtFunc.push_back(std::make_pair(p, sqrt(p)));
			std::map<int, int> iterations;
			for (int c = 1; c < p; c++) {
				int it = findCollision(p, x0, c);
				iterations[c] = it;
			}

			double mean = std::accumulate(iterations.begin(), iterations.end(), 0.0, [](double value, const std::pair<int, int>& p) { return value + p.second; }) / iterations.size();
			averageIterations[x0][p] = mean;
			if (mean > 3 * sqrt(p)) {
				file << "Unusually high iterations for p = " << p << ", x0 = " << x0 << ": " << mean << std::endl;
				auto minmax = std::minmax_element(iterations.begin(), iterations.end(), [](const std::pair<int, double>& p1, const std::pair<int, double>& p2) { return p1.second < p2.second; });
				file << "\t Maximum value: " << minmax.second->second << " at c = " << minmax.second->first << ", occured " << std::count(iterations.begin(), iterations.end(), Compare(minmax.second->second)) << " times" << std::endl;
				file << "\t Minimum value: " << minmax.first->second << " at c = " << minmax.first->first << ", occured " << std::count(iterations.begin(), iterations.end(), Compare(minmax.first->second)) << " times" << std::endl;
				std::vector<double> diff(iterations.size());
				std::transform(iterations.begin(), iterations.end(), diff.begin(), [mean](const std::pair<int, int>& p) { return p.second - mean; });
				file << "\t Standard deviation: " << sqrt(std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0) / diff.size()) << std::endl;
			}
		}
	}

	Gnuplot gp("gnuplot -persist");
	std::vector<double> minIterations;
	std::vector<double> maxIterations;
	for (const auto& averageIteration : averageIterations) {
		minIterations.push_back(std::min_element(averageIteration.second.begin(), averageIteration.second.end(), [](const std::pair<int, double>& p1, const std::pair<int, double>& p2) { return p1.second < p2.second; })->second);
		maxIterations.push_back(std::max_element(averageIteration.second.begin(), averageIteration.second.end(), [](const std::pair<int, double>& p1, const std::pair<int, double>& p2) { return p1.second < p2.second; })->second);
	}

	gp << "set xrange [" << minP << ':' << maxP << "]\n";
	gp << "set yrange [" << *std::min_element(minIterations.begin(), minIterations.end()) << ':' << *std::max_element(maxIterations.begin(), maxIterations.end()) << "]\n";
	gp << "plot";
	for (const auto& averageIteration : averageIterations) {
		std::string title = "with lines title 'avg. iterations across c from 1 to p, p from " + std::to_string(minP) + " to " + std::to_string(maxP) + ", x0 = " + std::to_string(averageIteration.first) + "',";
		gp << gp.file1d(averageIteration.second) << title;
	}
	gp << gp.file1d(sqrtFunc) << "with lines title 'sqrt'";
	gp << std::endl;

	file.close();
}

void testAverageC(int minX0, int maxX0, int minP, int maxP) {
	std::map<int, std::map<int, double>> averageIterations; // x0, c, average across p

	for (int x0 = minX0; x0 <= maxX0; x0++) {
		for (int c = -3; c < minP - 3; c++) {
			std::vector<int> iterations;
			for (int p = minP; p <= maxP; p++) {
				int it = findCollision(p, x0, c);
				iterations.push_back(it);
			}
			averageIterations[x0][c] = std::accumulate(iterations.begin(), iterations.end(), 0.0) / iterations.size();
		}
	}

	Gnuplot gp("gnuplot -persist");
	std::vector<double> minIterations;
	std::vector<double> maxIterations;
	for (const auto& averageIteration : averageIterations) {
		minIterations.push_back(std::min_element(averageIteration.second.begin(), averageIteration.second.end(), [](const std::pair<int, double>& p1, const std::pair<int, double>& p2) { return p1.second < p2.second; })->second);
		maxIterations.push_back(std::max_element(averageIteration.second.begin(), averageIteration.second.end(), [](const std::pair<int, double>& p1, const std::pair<int, double>& p2) { return p1.second < p2.second; })->second);
	}

	gp << "set xrange [-3:" << minP - 3 << "]\n";
	gp << "set yrange [" << *std::min_element(minIterations.begin(), minIterations.end()) << ':' << *std::max_element(maxIterations.begin(), maxIterations.end()) << "]\n";
	gp << "plot";
	for (const auto& averageIteration : averageIterations) {
		std::string title = "with lines title 'avg. iterations across p from " + std::to_string(minP) + " to " + std::to_string(maxP) + ", c from -3 to " + std::to_string(minP - 3) + ", x0 = " + std::to_string(averageIteration.first) + "',";
		gp << gp.file1d(averageIteration.second) << title;
	}
	gp << std::endl;
}

void visualizeIterationsP(int x0, int c) {
	// Dont't avoid 0, -2 (might want to test despite high iterations), staying at x0 (through negative constants) would cause failure in N and occurs when c = -x0 * (x0 +- 1). c = p - x0 * (x0 +- 1) (positive constants) is ok
	if (c == -x0 * (x0 + 1) || c == -x0 * (x0 - 1)) {
		std::cout << "Trivial case where iterators always stay at starting value x0 regardless of the value of p. Can't visualize iterations for given constant. Avoid c = -x0 * (x0 +- 1)." << std::endl;
		return;
	}

	std::map<int, int> iterations;
	std::vector<std::pair<int, double>> sqrtFunc;
	int maxP = 10000;
	for (int p = 1; p <= maxP; p++) {
		int it = findCollision(p, x0, c);
		iterations[p] = it;
		sqrtFunc.push_back(std::make_pair(p, sqrt(p)));
	}

	Gnuplot gp("gnuplot -persist");
	gp << "set xrange [1:" << maxP << "]\n";
	gp << "set yrange [" << std::min_element(iterations.begin(), iterations.end(), [](const std::pair<int, int>& p1, const std::pair<int, int>& p2) { return p1.second < p2.second; })->second << ":" <<
		std::max_element(iterations.begin(), iterations.end(), [](const std::pair<int, int>& p1, const std::pair<int, int>& p2) { return p1.second < p2.second; })->second << "]\n";
	gp << "plot" << gp.file1d(iterations) << "with lines title 'iterations across p, c = " + std::to_string(c) + ", x0 = " + std::to_string(x0) +
		", avg. = " + std::to_string(std::accumulate(iterations.begin(), iterations.end(), 0.0, [](double value, const std::pair<int, int>& p) { return value + p.second / sqrt(p.first); }) / iterations.size()) + "',"; // avg. relative difference between iterations and sqrt(p)
	gp << gp.file1d(sqrtFunc) << "with lines title 'sqrt'" << std::endl;
}

void visualizeIterationsC(int x0, int p) {
	std::map<int, int> iterations;
	for (int c = 1; c < p; c++) {
		int it = findCollision(p, x0, c);
		iterations[c] = it;
	}

	Gnuplot gp("gnuplot -persist");
	gp << "set xrange [1:" << p - 1 << "]\n";
	gp << "set yrange [" << std::min_element(iterations.begin(), iterations.end(), [](const std::pair<int, int>& p1, const std::pair<int, int>& p2) { return p1.second < p2.second; })->second << ":" <<
		std::max_element(iterations.begin(), iterations.end(), [](const std::pair<int, int>& p1, const std::pair<int, int>& p2) { return p1.second < p2.second; })->second << "]\n";
	gp << "plot" << gp.file1d(iterations) << "with lines title 'iterations across c, p = " + std::to_string(p) + ", x0 = " + std::to_string(x0) + "'" << std::endl;
}