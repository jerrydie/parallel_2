#include <iostream>
#include <random>
#include <limits>
#include <mpi.h>

#define FIRST_THREAD 0
#define MAX_NUM 50
#define VECTOR_SIZE 7

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
		for(std::size_t i = 0; i < v_size; i++)
		{
			v1[i] = generate_int();
			v2[i] = generate_int();
		}
		int rank, size, n, ibeg, iend;
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		int LAST_THREAD = FIRST_THREAD + size -1;
		if(rank == FIRST_THREAD) {
		  std::cout << "Scalar multiplication of:\nv1 = [";
		  for(int i = 0; i<v_size-1; i++)
		  	std::cout << v1[i] << ", ";
		  std::cout << v1[v_size-1] << "]\nv2 = [";
		  for(int i = 0; i<v_size-1; i++)
		  	std::cout << v2[i] << ", ";
		  std::cout << v2[v_size-1] << "]\n\n";
		}
		int result = 0;	
		n = (v_size - 1) / size + 1;
		ibeg = rank * n;
		iend = (rank + 1) * n;
		std::cout << "Process: " << rank << '\n';
		for(int i = ibeg; i < ((iend > v_size) ? v_size : iend); i++)
		{
		  result += v1[i]*v2[i];
		}
		if(rank != LAST_THREAD) {
		  	MPI_Send(&result, 1, MPI_INT, LAST_THREAD, 0, MPI_COMM_WORLD);
		}
		if (rank == LAST_THREAD) {
		  	MPI_Status status;
		  	for(int sender = 0; sender < rank; sender++) {
		  		int sender_result;
		  		MPI_Recv(&sender_result, 1, MPI_INT, sender, 0, MPI_COMM_WORLD, &status);
		  		result += sender_result;
		  	}
		  	std::cout << "(v1,v2) = " << result << '\n';
		}
		MPI_Finalize();
		return 0;		
	}
}

int main(int argc, char **argv)
{
  hse::parallel_2::scalar_mul(argc, argv, VECTOR_SIZE);
  return 0;
}
