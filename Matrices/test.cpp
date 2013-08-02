#include <iostream>
#include "DenseMatrix2D.h"
#include <ctime>


double timeElapsed(std::clock_t start, std::clock_t end)
{
	return (double) (end - start) / CLOCKS_PER_SEC * 1000.0;
}

int main(int argc, char const *argv[])
{
	DenseMatrix2D<int> mat1 (1000,1000, true);

	// for (int i =0 ; i < 1000 ; i++)
	// {
	// 	for (int j =0; j < 100; j++)
	// 	{
	// 		mat1(i,j) = rand()%10;
	// 	}
	// }

	// // std::cout<< mat1 << std::endl;
	// // DenseMatrix2D<int> mat2 (mat1);
	// // std::cout<< mat2 << std::endl;

	// std::clock_t time_start = std::clock();
	// for (int i = 0 ; i < 5000; i++)
	// {
	// 	DenseMatrix2D<int> mat2 (mat1);
	// }
	// std::clock_t time_end = std::clock();

	// std::cout<< timeElapsed(time_start,time_end) << std::endl;

	
	return 0;
}