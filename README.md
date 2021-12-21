Данный отчет и все связанные с ним материалы можно просмотреть по ссылке: https://github.com/jerrydie/parallel_2

# Основы MPI

## Вариант задания и используемые технологии

Задача 2.1: Номер в таблице 24: вариант 1 - задача скалярного произведения двух векторов.

Задача 2.2: Оптимизировать задачу интегрирования на отрезке с помощью MPI (в дополнение к предыдущим способам оптимизации).


Данный отчет включает: 

+ листинг программы на С++ (`src/main.cpp, src2/main.cpp, src2/format.hpp, src2/format.cpp`);
+ оценки производительности.


## Характеристики системы
```
cmd: lscpu output
```

|Характеристика                   |Значение
| ------------------------------- | -----------------------------------------
|Архитектура:                     |x86_64
|CPU op-mode(s):                  |32-bit, 64-bit
|Порядок байт:                    |Little Endian
|Address sizes:                   |39 bits physical, 48 bits virtual
|CPU(s):                          |2
|On-line CPU(s) list:             |0,1
|Потоков на ядро:                 |1
|Ядер на сокет:                   |2
|Сокетов:                         |1
|NUMA node(s):                    |1
|ID прроизводителя:               |GenuineIntel
|Семейство ЦПУ:                   |6
|Модель:                          |158
|Имя модели:                      |Intel(R) Core(TM) i7-7700HQ CPU @ 2.80GHz
|Степпинг:                        |9
|CPU МГц:                         |2808.002

Компилятор: g++ (Ubuntu 9.3.0-17ubuntu1~20.04) 9.3.0

Для данной задачи на виртуальной машине было увеличено количество ядер (по сравнение с первым заданием), все представленные далее программы написаны с учетом возможности использования любого количества процессов, однако тестировалось только на максимально возможном для виртуальной машине - 2.
## 2.1

Листинг программы можно найти в `src/main.cpp`, в `src\simple_main.cpp` показан пример работы алгоритма для одной пары сгенерированных векторов заданного размера (`VECTOR_SIZE 7`), а в `src_with_array` продемонстрирована корректность работы алгоритма с помощью вывода массива попарно перемноженных координат заданной пары векторов.
В конце данного отчета приведен только `src/main.cpp`.

Продемонстрируем работу `src/main.cpp`:

```cmd
polina@polina-VirtualBox:~/Prog/parallel_2/src$ mpic++ main.cpp -o main
polina@polina-VirtualBox:~/Prog/parallel_2/src$ mpiexec -n 2 ./main
```
Вывод программы с замерами (в случае, если таймеры разных процессов синхронизированы):

```
(v1,v2) = 3022
Vector size: 5
Total time: 0
GFLOPS: inf
Processor clocks: 445086
Result: 3022

(v1,v2) = 2465
Vector size: 5
Total time: 0.000399099
GFLOPS: 25056.4
Processor clocks: 1569777
Result: 2465

(v1,v2) = 2378
Vector size: 5
Total time: 0.000801189
GFLOPS: 12481.4
Processor clocks: 2698791
Result: 2378

(v1,v2) = 3320
Vector size: 5
Total time: 0.00109786
GFLOPS: 9108.67
Processor clocks: 3532932
Result: 3320

(v1,v2) = 3274
Vector size: 5
Total time: 0.00139919
GFLOPS: 7147
Processor clocks: 4378249
Result: 3274

(v1,v2) = 7354
Vector size: 10
Total time: 0.00167329
GFLOPS: 29881.3
Processor clocks: 5147732
Result: 7354

(v1,v2) = 1931
Vector size: 10
Total time: 0.0019817
GFLOPS: 25230.9
Processor clocks: 6013796
Result: 1931

(v1,v2) = 7051
Vector size: 10
Total time: 0.00225453
GFLOPS: 22177.6
Processor clocks: 6779983
Result: 7051

(v1,v2) = 4888
Vector size: 10
Total time: 0.00256473
GFLOPS: 19495.2
Processor clocks: 7650842
Result: 4888

(v1,v2) = 7754
Vector size: 10
Total time: 0.00308852
GFLOPS: 16189
Processor clocks: 9121708
Result: 7754

(v1,v2) = 36253
Vector size: 50
Total time: 0.00340254
GFLOPS: 367372
Processor clocks: 10003242
Result: 36253

(v1,v2) = 25624
Vector size: 50
Total time: 0.00369011
GFLOPS: 338743
Processor clocks: 10811009
Result: 25624

(v1,v2) = 33056
Vector size: 50
Total time: 0.00396238
GFLOPS: 315467
Processor clocks: 11575508
Result: 33056

(v1,v2) = 31746
Vector size: 50
Total time: 0.00427369
GFLOPS: 292487
Processor clocks: 12449475
Result: 31746

(v1,v2) = 25160
Vector size: 50
Total time: 0.00477209
GFLOPS: 261940
Processor clocks: 13849108
Result: 25160

(v1,v2) = 51752
Vector size: 100
Total time: 0.00505738
GFLOPS: 988654
Processor clocks: 14650179
Result: 51752

(v1,v2) = 65769
Vector size: 100
Total time: 0.00536211
GFLOPS: 932469
Processor clocks: 15506555
Result: 65769

(v1,v2) = 64224
Vector size: 100
Total time: 0.00565061
GFLOPS: 884861
Processor clocks: 16315910
Result: 64224

(v1,v2) = 57718
Vector size: 100
Total time: 0.00592986
GFLOPS: 843190
Processor clocks: 17100093
Result: 57718

(v1,v2) = 69063
Vector size: 100
Total time: 0.00624223
GFLOPS: 800995
Processor clocks: 17977422
Result: 69063

(v1,v2) = 642421
Vector size: 1000
Total time: 0.00646382
GFLOPS: 7.73537e+07
Processor clocks: 18599476
Result: 642421

(v1,v2) = 615256
Vector size: 1000
Total time: 0.00587313
GFLOPS: 8.51335e+07
Processor clocks: 16939881
Result: 615256

(v1,v2) = 646661
Vector size: 1000
Total time: 0.00560715
GFLOPS: 8.91718e+07
Processor clocks: 16193765
Result: 646661

(v1,v2) = 634466
Vector size: 1000
Total time: 0.00532477
GFLOPS: 9.39007e+07
Processor clocks: 15401176
Result: 634466

(v1,v2) = 624510
Vector size: 1000
Total time: 0.00522021
GFLOPS: 9.57815e+07
Processor clocks: 15107460
Result: 624510
```

## 2.2
Напомним результаты оптимизации (без MPI) задачи интегрирования на отрезке:

```cmd
polina@polina-VirtualBox:~/Prog/parallel_2/src2$ g++ main_pragma.cpp format.cpp -g -O3 -march=native -o test.o -fopenmp
polina@polina-VirtualBox:~/Prog/parallel_2/src2$ ./test.o
Type lower integration limit a: 1
Type upper integration limit b: 3
```
### SIMPSON INTEGRATION from a: 1 to b: 3
|   SPLIT SEGMENTS   |        TIME        |       GFLOPS       |  PROCESSOR CLOCKS  |  OPERATION RESULT  |
|--------------------|--------------------|--------------------|--------------------|--------------------|
|        100         |      0.000235      |   425694.947001    |       658081       |     109.427734     |
|        1000        |      0.000116      |   8610001.377600   |       325808       |     137.834219     |
|       10000        |      0.000586      |  17064962.900771   |      1645018       |     105.283878     |
|       100000       |      0.005879      |  17008850.725564   |      16506053      |     100.871242     |
|      1000000       |      0.073858      |  13539546.951865   |     207387988      |     119.008871     |


Теперь доработаем программу, добавив использование MPI (Листинг программы содержится в `src2/main.cpp`, а также в приложении к данному отчету):

```cmd
polina@polina-VirtualBox:~/Prog/parallel_2/src2$ mpic++ main.cpp format.cpp -g -O3 -march=native -o test.o -fopenmp
polina@polina-VirtualBox:~/Prog/parallel_2/src2$ mpiexec -n 2 ./test.o
Type lower integration limit a: 1
Type upper integration limit b: 3
```
Вывод программы с замерами (в случае, если таймеры разных процессов синхронизированы):

### SIMPSON INTEGRATION from a: 1 to b: 3
|   SPLIT SEGMENTS   |        TIME        |       GFLOPS       |  PROCESSOR CLOCKS  |  OPERATION RESULT  |
|--------------------|--------------------|--------------------|--------------------|--------------------|
|        100         |      0.000000      |        inf         |       152042       |     78.279964      |
|        1000        |      0.000050      |  19843238.416510   |       298014       |     73.066667      |
|       10000        |      0.000356      |  28091860.383454   |      1155535       |     73.066667      |
|       100000       |      0.002242      |  44604246.591901   |      6451373       |     73.066667      |
|      1000000       |      0.023706      |  42183755.137138   |      66718799      |     73.066667      |

По сравнению с оптимизацией без MPI на больших объемах вычислений видно заметное повышение производительности. 

# WARNING!!!
Выявлен баг - без MPI интегрирование вычислено неверно! Хотя функция та же самая.

## Листинги

```C++
//src/main.cpp:

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
```
```C++
//src2/main.cpp:

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
```

```C++
//src2/format.hpp:

#ifndef __FORMAT_HPP__
#define __FORMAT_HPP__

#include <iostream>

namespace hse::parallel::lab1
{
	std::string make_break();
	
	std::string make_cell(const std::string& content);

	std::string make_header();
	
	struct Line_state
	{
		std::size_t n;
		double secs;
		double ops_per_sec;
		std::uint64_t tact_duration;
		double result;
		friend std::ostream& operator <<(std::ostream& os, Line_state const& line);
	};
}
#endif
```

```C++
//src2/format.cpp:

#include "format.hpp"
#define CELL_WIDTH 20
#define COL_NUM 5
#include <iostream>

namespace hse::parallel::lab1
{
	std::string make_break()
	{
		std::string header = "";
		for (int i = 0; i < COL_NUM; i++)
		{
			header += '|' + std::string(CELL_WIDTH, '-');
		}
		return header + "|\n";
	}
	
	std::string make_cell(const std::string& content)
	{
		std::size_t len = content.length();
		std::size_t left_len = (CELL_WIDTH-len)/2;
		std::size_t right_len = (CELL_WIDTH-len) - left_len;
		return std::string(left_len, ' ') + content + std::string(right_len, ' ');
	}

	std::string make_header()
	{
		return  '|' + make_cell("SPLIT SEGMENTS")+
			'|' + make_cell("TIME")+ 
			'|' + make_cell("GFLOPS")+
			'|' + make_cell("PROCESSOR CLOCKS")+
			'|' + make_cell("OPERATION RESULT")+ "|\n" + make_break();
	}
	
	
	std::ostream& operator <<(std::ostream& os, Line_state const& line)
		    {
			return os <<  '|' << make_cell(std::to_string(line.n))<< '|'
				  <<  make_cell(std::to_string(line.secs))<< '|'
				  <<  make_cell(std::to_string(line.ops_per_sec))<< '|'
				  <<  make_cell(std::to_string(line.tact_duration))<< '|'
				  <<  make_cell(std::to_string(line.result)) << "|\n";
		    }

}
```
