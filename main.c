#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

extern int make_csv(double** u, int segmentSize, int t_count, int x_count, int commsize, int my_rank);
extern double** solve_scheme(int segmentSize, int t_count, int commsize, int my_rank, double courant, double* timeBoundary, double* distBoundary);

extern double t_max;
extern double x_max;
extern double t_step;
extern double x_step;

int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);
	int my_rank, commsize;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &commsize);

	// Start counting time inside process zero
	// double startTime;
	// if (my_rank == 0)
	// 	startTime = MPI_Wtime();

	int t_count = (int) t_max/t_step;
	int x_count = (int) x_max/x_step;

	double a = 1.0;
	double* timeBoundary = (double*) calloc(t_count, sizeof(double));
	double* distBoundary = (double*) calloc(x_count, sizeof(double));
	int i;
	
	for (i = 0; i < t_count; i++)
		timeBoundary[i] = sin(2. * M_PI * t_step * i);
	for (i = 0; i < x_count; i++)
		distBoundary[i] = sin(2. * M_PI * x_step * i);
	
	double courant = a * t_step / x_step;

	int segmentSize = (int) x_max/(x_step * commsize);
	int startOfSegment = my_rank * segmentSize;

	double** u = solve_scheme(segmentSize, t_count, commsize, my_rank, courant, timeBoundary, distBoundary);

	// Count time
	// double endTime;
	// if (my_rank == 0){
	// 	endTime = MPI_Wtime();
	// 	printf("Time: %lg.\n", endTime - startTime);
	// }

	make_csv(u, segmentSize, t_count, x_count, commsize, my_rank);

	// Free memory
	free(timeBoundary);
	free(distBoundary);
	for (i = 0; i < segmentSize; i++)
		free(u[i]);
	free(u);

	MPI_Finalize();
	return 0;
}