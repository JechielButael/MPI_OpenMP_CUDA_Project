/*
 ============================================================================
 Name        : finalProject.c
 Author      : Yehiel Butael
 Version     : final
 Copyright   : Your copyright notice
 Description : Hello OpenMP World in C
 ============================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "functions.h"
#include "cuda.h"
#include <mpi.h>
#include <omp.h>

#define ROOT 0
#define WEIGHT_SIZE 4

#define FILE_NAME "/home/linuxu/Downloads/finalProject_315016774/finalProject/src/input.txt"


int main(int argc, char *argv[]) {

	//variables for reading file
	float *weights = NULL;
	char *seq1 = NULL;
	char **arrOfSeq2 = NULL;
	int numOfSeq2 = 0;
	int seq1Length = 0;
	int* subSeq2Lengths = NULL;

	//general variables for processors
	double startTime, finishTime;
	Mutant mutant;
	int numOfOffsets = 0;
	Mutant* allMutants;
	Mutant rankers[4];


	//more general variables
	MPI_Status status;
	int my_rank, num_procs,inputSize;
	MPI_Datatype mutantMPIType;
	MPI_Datatype type[4] = { MPI_INT, MPI_INT , MPI_INT , MPI_FLOAT};
	int blocklen[4] = { 1, 1, 1 ,1};
	MPI_Aint disp[4];


	//MPI INIT
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	// NEED ONLY 2 PROCESS
	if(num_procs!=2)
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);

	// Create MPI user data type for partical
	disp[0] = (char *) &mutant.n -	 (char *) &mutant;
	disp[1] = (char *) &mutant.k -	 (char *) &mutant;
	disp[2] = (char *) &mutant.offset -	 (char *) &mutant;
	disp[3] = (char *) &mutant.score - (char *) &mutant;
	MPI_Type_create_struct(WEIGHT_SIZE , blocklen, disp, type, &mutantMPIType);
	MPI_Type_commit(&mutantMPIType);



	if(my_rank == ROOT) //PROCESS 0
	{
		//starting timer
		startTime = MPI_Wtime();

		//reading from file
		arrOfSeq2 = readFromFile(FILE_NAME, &weights, &seq1 , &numOfSeq2 , &seq1Length , &subSeq2Lengths);
		Mutant finalScores[numOfSeq2];
		inputSize = numOfSeq2/2;


		//SEND THE FIRST HALF TO PROCESS 1
		MPI_Send(weights, WEIGHT_SIZE ,MPI_FLOAT,1,0,MPI_COMM_WORLD);
		MPI_Send(&seq1Length,1 ,MPI_INT,1,0,MPI_COMM_WORLD);
		MPI_Send(seq1,seq1Length,MPI_CHAR,1,0,MPI_COMM_WORLD);
		MPI_Send(&inputSize, 1 ,MPI_INT,1,0,MPI_COMM_WORLD);
		for(int i = 0 ; i < inputSize ; i++)
		{
			MPI_Send(&subSeq2Lengths[i], 1 ,MPI_INT,1,0,MPI_COMM_WORLD);
			MPI_Send(arrOfSeq2[i],subSeq2Lengths[i],MPI_CHAR,1,0,MPI_COMM_WORLD);
		}


		Mutant currentMutant , bestMutant;
		int numOfMutants;


		//calculating half of the array of seq2
		for(int i = inputSize ; i < numOfSeq2 ; i++)
		{
			//initializing rankers
			for (int k = 0; k < 4; k++) {
				rankers[k].score = -5000;
			}
			numOfOffsets = seq1Length - (subSeq2Lengths[i]-2) +1;
			allMutants = NULL;
			currentMutant.score = 0;
			bestMutant.score = 0;

			//filling n,k for each mutant
			preprareMutants(&allMutants , subSeq2Lengths[i]);

			omp_set_num_threads(4);


			//running threads to calculate best mutant for each offset
			#pragma omp parallel for private(currentMutant)
			for(int j = 0 ; j < numOfOffsets ; j++)
			{
				currentMutant.offset = j;
				//cuda - calculating best mutant for a given offset
				cudaKernal(allMutants , seq1, arrOfSeq2[i] ,subSeq2Lengths[i] ,seq1Length,  weights , &currentMutant);
				if(currentMutant.score > rankers[omp_get_thread_num()].score)
				{
					//updating scores according to thread_num to prevent collisions
					rankers[omp_get_thread_num()] = currentMutant;
				}
			}


		//finding best mutant out of all the 4 threads
		findMaxMutant(rankers , &bestMutant);
		finalScores[i] = bestMutant;
		}
		for(int i = 0 ; i < inputSize ; i++)
		{
			//receiving  the other half of the array of seq2
			MPI_Recv(&finalScores[i] ,1,mutantMPIType,1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		}


		//printing best mutants 
		printBestMutants(finalScores , numOfSeq2);

		//shutting down timer
		finishTime = MPI_Wtime();
		printf("Execution Time: %f\n", finishTime - startTime);



	}else{ //PROCESS 1
		int temp;

		//RECEIVE THE FIRST HALF FROM PROCESS 0
		weights = (float*)malloc(WEIGHT_SIZE*sizeof(float));
		MPI_Recv(weights,WEIGHT_SIZE,MPI_FLOAT,ROOT,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		MPI_Recv(&seq1Length,1 ,MPI_INT,ROOT,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		seq1 = (char*)malloc(seq1Length*sizeof(char));
		MPI_Recv(seq1, seq1Length ,MPI_CHAR,ROOT,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		MPI_Recv(&inputSize, 1 ,MPI_INT ,ROOT,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		arrOfSeq2 = (char**)malloc(inputSize*sizeof(char*));
		subSeq2Lengths = (int*)malloc(inputSize*sizeof(int));

		for(int i = 0 ; i < inputSize ; i++)
		{
			MPI_Recv(&temp, 1 ,MPI_INT ,ROOT,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			subSeq2Lengths[i] = temp;
			arrOfSeq2[i] = (char*)malloc(temp*sizeof(char));
			MPI_Recv(arrOfSeq2[i], temp ,MPI_CHAR ,ROOT,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		}

		
		Mutant currentMutant , bestMutant;
		Mutant finalScores[inputSize];
		int numOfMutants;


		//calculating half of the array of seq2
		for(int i = 0 ; i < inputSize ; i++)
		{
			//initializing rankers
			for (int k = 0; k < 4; k++) {
				rankers[k].score = -5000;
			}
			numOfOffsets = seq1Length - (subSeq2Lengths[i]-2) +1;
			
			allMutants = NULL;
			currentMutant.score = 0;
			bestMutant.score = 0;

			//filling n,k for each mutant
			preprareMutants(&allMutants , subSeq2Lengths[i]);

			omp_set_num_threads(4);


			//running threads to calculate best mutant for each offset
			#pragma omp parallel for private(currentMutant)
			for(int j = 0 ; j < numOfOffsets ; j++)
			{
				currentMutant.offset = j;

				//cuda - calculating best mutant for a given offset
				cudaKernal(allMutants , seq1, arrOfSeq2[i] ,subSeq2Lengths[i] ,seq1Length,  weights ,&currentMutant);
				if(currentMutant.score > rankers[omp_get_thread_num()].score)
				{
					//updating scores according to thread_num to prevent collisions
					rankers[omp_get_thread_num()] = currentMutant;
				}
			}

			//finding best mutant out of all the 4 threads
			findMaxMutant(rankers , &bestMutant);
			finalScores[i] = bestMutant;
		}


		//sending the other half of the array of seq2
		for(int i = 0 ; i < inputSize ; i++){
			MPI_Send(&finalScores[i],1,mutantMPIType ,ROOT,my_rank,MPI_COMM_WORLD);
		}

	}


	//free all allocations
	free(weights);
	free(allMutants);
	free(seq1);
	for(int i = 0 ; i < numOfSeq2 ; i++)
		free(arrOfSeq2[i]);
	free(subSeq2Lengths);
	free(arrOfSeq2);
	MPI_Finalize();
	return 0;
}
