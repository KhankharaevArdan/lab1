#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

double t_max = 1;
double x_max = 1;
double t_step = 1e-2;
double x_step = 1e-2;

double** solve_scheme(int segmentSize, int t_count, int commsize, int my_rank, double courant, double* timeBoundary, double* distBoundary) {
    double** u = (double**) calloc(segmentSize, sizeof(double*));
	for (int i = 0; i < segmentSize; i++)
		u[i] = (double*) calloc(t_count, sizeof(double));

	for (int i = 0; i < segmentSize; i++)
		u[i][0] = distBoundary[i + my_rank * segmentSize];

	for (int i = 1; i < segmentSize; i++)
		u[i][1] = courant * u[i-1][0] + (1. - courant) * u[i][0];

	if (my_rank != commsize - 1)
		MPI_Send(&u[segmentSize-1][0], 1, MPI_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD);

	if (my_rank != 0) {
		double tmp;
		MPI_Recv(&tmp, 1, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		u[0][1] = courant * tmp + (1. - courant) * u[0][0];
	}

	if (my_rank == 0)
		u[0][1] = timeBoundary[1];
	

	int j; double tmp_behind, tmp_front;
	for (int i = 2; i < t_count; i++) {
		for (j = 1; j < segmentSize - 1; j++) {
			u[j][i] = u[j][i - 2] + courant * (u[j-1][i - 1] - u[j+1][i - 1]);
		}

		if (my_rank != commsize - 1)
			MPI_Send(&u[segmentSize-1][i - 1], 1, MPI_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD);
		if (my_rank != 0) {
			MPI_Recv(&tmp_behind, 1, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Send(&(u[0][i - 1]), 1, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD);
		}
		if (my_rank != commsize - 1) {
			MPI_Recv(&tmp_front, 1, MPI_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			u[segmentSize-1][i] = u[segmentSize-1][i - 2] + courant * (u[segmentSize-2][i - 1] - tmp_front);
		}
		if (my_rank != 0)
			u[0][i] = u[0][i - 2] + courant * (tmp_behind - u[1][i - 1]);
		if (!my_rank)
			u[0][i] = timeBoundary[i];
		if (my_rank == commsize - 1)
			u[segmentSize-1][i] = courant * u[segmentSize-2][i - 1] + (1. - courant) * u[segmentSize-1][i - 1];
	}
    return u;
}

int make_csv(double** u, int segmentSize, int t_count, int x_count, int commsize, int my_rank) {
    if (my_rank != 0)
		for (int i = 0; i < segmentSize; i++)
			MPI_Send(u[i], t_count, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

	if (my_rank == 0) {
		double** solution = (double**) calloc(x_count, sizeof(double*));
		for (int i = 0; i < x_count; i++)
			solution[i] = (double*) calloc(t_count, sizeof(double));

		for (int j = 0; j < segmentSize; j++)
			memcpy(solution[j], u[j], t_count * sizeof(double));
		if (commsize > 1) {
			for (int j = segmentSize; j < x_count; j++) {
				MPI_Recv(solution[j], t_count, MPI_DOUBLE, (int)(j / segmentSize), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
				
		}      

		FILE *fp = fopen("data.csv", "w");
        if (fp == NULL) {
            printf("Error opening file!\n");
            return 1;
        }

		for (int j = 0; j < x_count; j++) {
			for (int i = 0; i < t_count; i++)
				fprintf (fp, "%lf ", solution[j][i]);
			fprintf(fp, "\n");
			free(solution[j]);
		}
		
		free (solution);
	}
    return 0;
}