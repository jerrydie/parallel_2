#include <iostream>
#include <functional>
#include <immintrin.h>
#include <chrono>
#include <cmath>
#include <omp.h>
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
	    auto g = [](double x1, double x2, const std::function<double (double)> &f_)->double{ return (x2-x1)/6.0*(f_(x1) + 4.0*f_(0.5*(x1+x2)) + f_(x2)); };
	    double simpson_integral_tmp[4] = {0.,0.,0.,0.};
	    __m256d val1 = _mm256_load_pd(&simpson_integral_tmp[0]);
	    __m256d val2;
	    double x1 = a;
	    #pragma omp parallel
	    for( std::size_t step = 0; step < n ; step+=4) {
	        double x1 = a + step*width;
	        double x2 = x1 + width;
	        double x3 = x2 + width;
	        double x4 = x3 + width;
	        double x5 = x4 + width;
	        val2 = _mm256_set_pd (g(x1, x2, f),g(x2, x3, f),g(x3, x4, f),g(x4, x5, f));
	        val1 = _mm256_add_pd (val1, val2);
	    }
	    _mm256_store_pd (simpson_integral_tmp, val1);
	    double simpson_integral = simpson_integral_tmp[0] + simpson_integral_tmp[1] + simpson_integral_tmp[2] + simpson_integral_tmp[3];
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
