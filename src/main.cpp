#include <iostream>
#include <random>
#include <limits>
#include <mpi.h>
#include <chrono>

#define FIRST_THREAD 0
#define MAX_NUM 50

namespace hse::parallel_2
{
	std::uint64_t rdtsc() {
	    unsigned int lo, hi;
	    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	    return ((std::uint64_t) hi << 32) | lo;
	}

	static std::mt19937 prng(std::random_device{}());

	std::int32_t generate_int()
	{
	    std::uniform_int_distribution<std::int32_t> distrib(0, MAX_NUM);
	    return distrib(prng);
	}
	
	
	int scalar_mul(int argc, char **argv, std::size_t v_size)
	{
		std::int32_t v1 [v_size];
		std::int32_t v2 [v_size];
		// Генерация векторов
		for(std::size_t i = 0; i < v_size; i++)
		{
			v1[i] = generate_int();
			v2[i] = generate_int();
		}
		int rank, size, n, ibeg, iend;
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		int LAST_THREAD = FIRST_THREAD + size -1;
		int result = 0;	
		n = (v_size - 1) / size + 1;
		ibeg = rank * n;
		iend = (rank + 1) * n;
		// Начало замера
		if(rank == FIRST_THREAD) {
			auto start = MPI_Wtime();
			std::uint64_t tact_start = rdtsc();
			MPI_Send(&start, 1, MPI_DOUBLE, LAST_THREAD, 1, MPI_COMM_WORLD);
			MPI_Send(&tact_start, 1, MPI_UINT64_T, LAST_THREAD, 2, MPI_COMM_WORLD);
		}
		
		for(int i = ibeg; i < ((iend > v_size) ? v_size : iend); i++)
		{
		  result += v1[i]*v2[i];
		}
		if(rank != LAST_THREAD) {
		  	MPI_Send(&result, 1, MPI_INT, LAST_THREAD, 0, MPI_COMM_WORLD);
		}
		else {
		  	MPI_Status status;
		  	double start;
			std::uint64_t tact_start;
		  	MPI_Recv(&start, 1, MPI_DOUBLE, FIRST_THREAD, 1, MPI_COMM_WORLD, &status);
		  	MPI_Recv(&tact_start, 1, MPI_UINT64_T, FIRST_THREAD, 2, MPI_COMM_WORLD, &status);
		  	for(int sender = 0; sender < rank; sender++) {
		  		int sender_result;
		  		MPI_Recv(&sender_result, 1, MPI_INT, sender, 0, MPI_COMM_WORLD, &status);
		  		result += sender_result;
		  	}
		  	// Конец замера
		  	std::uint64_t tact_end = rdtsc();
			auto end = MPI_Wtime();
			// Вывод результата
		  	std::cout << "(v1,v2) = " << result << '\n';
		  	double duration = end - start;
			auto tacts = tact_end - tact_start;
			// Вывод замеров
			std::cout << "Vector size: " << v_size;
			// Если таймеры разных процесов синхронизированы
			if(MPI_WTIME_IS_GLOBAL) {
				std::cout << "\nTotal time: " << duration;
				std::cout << "\nGFLOPS: " << v_size*(((iend > v_size) ? v_size : iend) - ibeg)/(duration);
				std::cout << "\nProcessor clocks: " << tacts;
			}
			std::cout << "\nResult: " << result << "\n\n";
		}
	
		return 0;		
	}
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  // Потворяем замеры по 5 раз на каждый размер вектора
  std::vector<std::size_t> sizes {5, 10, 50, 100, 1000};
  for (auto size : sizes) {
  	for( int i = 0; i<5 ; i++)
  		hse::parallel_2::scalar_mul(argc, argv, size);
  
  }
  MPI_Finalize();
  return 0;
}
