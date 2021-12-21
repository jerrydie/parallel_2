#include <iostream>
#include <functional>
#include <chrono>
#include <cmath>
#include "format.hpp"

namespace hse::parallel::lab1
{
	std::uint64_t rdtsc() {
	unsigned int lo, hi;
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	return ((std::uint64_t) hi << 32) | lo;
	}

	double simpsonIntegral(double a, double b, std::size_t n, const std::function<double (double)> &f) {
	    const double width = (b-a)/n;

	    double simpson_integral = 0;
	    for(std::size_t step = 0; step < n; step++) {
		const double x1 = a + step*width;
		const double x2 = a + (step+1)*width;

		simpson_integral += (x2-x1)/6.0*(f(x1) + 4.0*f(0.5*(x1+x2)) + f(x2));
	    }

	    return simpson_integral;
	}
}

int main()
{

	double a;
	double b;
    
	std::cout << "Type lower integration limit a: ";
	std::cin >> a;
	std::cout << "Type upper integration limit b: ";
	std::cin >> b;
	
	std::size_t max_n = pow(10,6);
	auto f = [](double x)->double{ return x*x*x*x +5.*x - 2./(x*x) +3.; };
	
	std::cout << "### SIMPSON INTEGRATION from a: " << a <<" to b: " << b <<'\n' << hse::parallel::lab1::make_header();
	for(std::size_t n = pow(10,2); n <= max_n; n *= 10)
	{
		auto start = std::chrono::system_clock::now();
		std::uint64_t tact_start = hse::parallel::lab1::rdtsc();
	    
		double res = hse::parallel::lab1::simpsonIntegral(a, b, n, f);

		std::uint64_t tact_end = hse::parallel::lab1::rdtsc();
		auto end = std::chrono::system_clock::now();

		std::chrono::duration<double> duration = end - start;
		auto tacts = tact_end - tact_start;
		hse::parallel::lab1::Line_state state {n, duration.count(), n/(duration.count()), tacts, res};
		
		std::cout << state;
	}
	return 0;
}
