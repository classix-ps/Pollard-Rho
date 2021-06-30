#pragma once

#include "factorize.hpp"

#include <boost/math/special_functions/sign.hpp>

void testPollardFactor(int n);

void findSameMod(int p, int q);

int findCollision(mpz_class n, mpz_class x0, mpz_class c);

void testAverageP(int minX0, int maxX0, int minP, int maxP);

void testAverageC(int minX0, int maxX0, int minP, int maxP);

void visualizeIterationsP(int x0, int c);

void visualizeIterationsC(int x0, int p);