/*
 * functions.c
 *
 *  Created on: 7 Mar 2022
 *      Author: linuxu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "cuda.h"
#include <mpi.h>
#include <omp.h>


char** readFromFile(const char *fileName, float **weights, char **seq1 , int* numOfSeq2
		, int* seq1Length , int** subSeq2Lengths)
{

	FILE *fp;
	int size = 0;
	char ch;
	float *weightTemp;
	int* tempSubSeq2Lengths;
	char *seq1Temp = NULL;
	char **arrOfSeq2;

	// Open file for reading
	if ((fp = fopen(fileName, "r")) == 0) {
		printf("cannot open file %s for reading\n", fileName);
		exit(0);
	}
	weightTemp = (float*) malloc(4 * sizeof(float));
	if (!(weightTemp)) {
		printf("Problem to allocate memory\n");
		exit(0);
	}

	// reading weights
	for (int i = 0; i < 4; i++) {
		fscanf(fp, "%f ", &weightTemp[i]);
	}

	*weights = weightTemp;

	//reading seq1
	while (EOF != (ch = fgetc(fp)) && ch != '\n') {
		size++;
		seq1Temp = (char*) realloc(seq1Temp, size * sizeof(char));
		if (!seq1Temp) {
			printf("Problem to allocate memory\n");
			exit(0);
		}
		seq1Temp[size - 1] = ch;
	}
	*seq1Length = size -1;
	*seq1 = seq1Temp;

	//reading arrary of seq2
	fscanf(fp, "%d ", numOfSeq2);
	arrOfSeq2 = (char**) malloc(*(numOfSeq2) * sizeof(char*));
	if (!arrOfSeq2) {
		printf("Problem to allocate memory\n");
		exit(0);
	}
	tempSubSeq2Lengths = (int*) malloc(*(numOfSeq2) * sizeof(int));
	if (!tempSubSeq2Lengths) {
		printf("Problem to allocate memory\n");
		exit(0);
	}
	for (int i = 0; i < *(numOfSeq2) ; i++)
	{
		arrOfSeq2[i] = NULL;
		size = 0;

		while (EOF != (ch = fgetc(fp)) && ch != '\n') {
			size++;
			arrOfSeq2[i] = (char*) realloc(arrOfSeq2[i], size * sizeof(char));
			if (!arrOfSeq2[i]) {
				printf("Problem to allocate memory\n");
				exit(0);
			}
			arrOfSeq2[i][size - 1] = ch;
		}
		tempSubSeq2Lengths[i] = size -1;
	}
	*subSeq2Lengths = tempSubSeq2Lengths;
	fclose(fp);
	return arrOfSeq2;
}


void printBestMutants(Mutant* finalScores , int numOfScores)
{
	for(int i = 0 ; i < numOfScores ; i++)
	{
		printf("Seq 2.%d:	MS(%d,%d)	Offset:%d	Score:%f \n" , i, finalScores[i].n , finalScores[i].k
				, finalScores[i].offset , finalScores[i].score);
	}
}


void preprareMutants(Mutant** allMutants , int sizeOfSeq2)
{
	int numOfMutants = (sizeOfSeq2*(sizeOfSeq2-1))/2;
	//allocating all mutants array according to the size
	*allMutants = (Mutant*)malloc(numOfMutants*sizeof(Mutant));						
	if(!(*allMutants)){
		fprintf(stderr, "Mutants Array Allocation Error\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	fillAllMutants(*allMutants , sizeOfSeq2);

}

void fillAllMutants(Mutant* allMutants , int size)
{
	int index = 0;
	for(int n = 0 ; n < (size-1) ; n++){
		for(int k = n+1 ; k < size ; k++){
			allMutants[index].n = n;
			allMutants[index].k = k;
			allMutants[index].score = -5000;
			index++;
		}
	}
}


void findMaxMutant(Mutant* rankers , Mutant* bestMutant)
{
	Mutant max = rankers[0];
	for (int i = 1; i < 4; i++) {
		if(rankers[i].score>max.score){
			max = rankers[i];
		}
	}
	bestMutant->score = max.score;
	bestMutant->k = max.k;
	bestMutant->n = max.n;
	bestMutant->offset = max.offset;

}