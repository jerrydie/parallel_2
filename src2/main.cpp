#include <iostream>
#include <functional>
#include <immintrin.h>
#include <chrono>
#include <cmath>
#include <omp.h>
#include <mpi.h>
#include "format.hpp"
#define FIRST_THREAD 0

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
	    double simpson_integral = 0;
	    double x1 = a;
	    for( std::size_t step = 0; step < n ; step+=4) {
	        double x1 = a + step*width;
	        double x2 = x1 + width;
	        double x3 = x2 + width;
	        double x4 = x3 + width;
	        double x5 = x4 + width;
	        val2 = _mm256_set_pd (g(x1, x2, f),g(x2, x3, f),g(x3, x4, f),g(x4, x5, f));
	        _mm256_store_pd (simpson_integral_tmp, val2);
	        #pragma omp parallel for reduction(+:simpson_integral)
	        for (auto elem: simpson_integral_tmp)
	        	simpson_integral += elem;
	    }
	    return simpson_integral;
	}
	
}

int main(int argc, char **argv)
{

	double a;
	double b;
    	MPI_Init(&argc, &argv);
    	int rank, size, n;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int LAST_THREAD = FIRST_THREAD + size -1;
	if(rank == FIRST_THREAD) {
		std::cout << "Type lower integration limit a: ";
		std::cin >> a;
		std::cout << "Type upper integration limit b: ";
		std::cin >> b;
		std::cout << "### SIMPSON INTEGRATION from a: " << a <<" to b: " << b <<'\n' << hse::parallel::lab1::make_header();
		for (int to_thread = 1; to_thread < size; to_thread++) {
			MPI_Send(&a, 1, MPI_DOUBLE, to_thread, 1, MPI_COMM_WORLD);
			MPI_Send(&b, 1, MPI_DOUBLE, to_thread, 2, MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Status status;
		MPI_Recv(&a, 1, MPI_DOUBLE, FIRST_THREAD, 1, MPI_COMM_WORLD, &status);
		MPI_Recv(&b, 1, MPI_DOUBLE, FIRST_THREAD, 2, MPI_COMM_WORLD, &status);
	}
	std::size_t max_n = pow(10,6)/size;
	double thread_a = a + (b-a)*rank/size;
	double thread_b = a + (b-a)*(rank+1)/size;
	auto f = [](double x)->double{ return x*x*x*x +5.*x - 2./(x*x) +3.; };
	
	for(std::size_t n = pow(10,2)/size; n <= max_n; n *= 10)
	{
		if(rank == FIRST_THREAD) {
			auto start = MPI_Wtime();
			std::uint64_t tact_start = hse::parallel::lab1::rdtsc();
			MPI_Send(&start, 1, MPI_DOUBLE, LAST_THREAD, 3, MPI_COMM_WORLD);
			MPI_Send(&tact_start, 1, MPI_UINT64_T, LAST_THREAD, 4, MPI_COMM_WORLD);
		}
		
		double res = hse::parallel::lab1::simpsonIntegral(thread_a, thread_b, n, f);
		if(rank != LAST_THREAD)
			MPI_Send(&res, 1, MPI_DOUBLE, LAST_THREAD, 0, MPI_COMM_WORLD);
		else {
			MPI_Status status;
		  	double start;
			std::uint64_t tact_start;
		  	MPI_Recv(&start, 1, MPI_DOUBLE, FIRST_THREAD, 3, MPI_COMM_WORLD, &status);
		  	MPI_Recv(&tact_start, 1, MPI_UINT64_T, FIRST_THREAD, 4, MPI_COMM_WORLD, &status);
		  	for(int sender = 0; sender < rank; sender++) {
		  		double sender_result;
		  		MPI_Recv(&sender_result, 1, MPI_DOUBLE, sender, 0, MPI_COMM_WORLD, &status);
		  		res += sender_result;
		  	}
		  	// Конец замера
		  	std::uint64_t tact_end = hse::parallel::lab1::rdtsc();
			auto end = MPI_Wtime();
			
		  	double duration = end - start;
			auto tacts = tact_end - tact_start;
			// Если таймеры разных процесов синхронизированы
			hse::parallel::lab1::Line_state state;
			if(MPI_WTIME_IS_GLOBAL)
				state = {n*size, duration, (n*size)/(duration), tacts, res};
			else
				state = {n*size, 0, 0, 0, res};
			std::cout << state;
		}
	}
	MPI_Finalize();
	return 0;
}
