#include <iostream>
#include <vector>
#include <mpi.h>

#define FIRST_THREAD 0

int main(int argc, char **argv)
{
  int v1 [7] {1,2,3,4,5,6,7};
  int v2 [7] {3,2,4,1,6,3,2};
  if((sizeof(v1)/sizeof(*v1)) != (sizeof(v2)/sizeof(*v2))) {
  	std::cerr << "Scalar multiplication impossible\n";
  	return 1;
  }
  int rank, size, n, ibeg, iend;
  int v_size = sizeof(v1)/sizeof(*v1);
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
