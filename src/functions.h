/*
 * functions.h
 *
 *  Created on: 7 Mar 2022
 *      Author: linuxu
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_


typedef struct{
	int n;
	int k;
	int offset;
	float score;
}Mutant;

char** readFromFile(const char *fileName, float **weights, char **seq1 , int* numOfSeq2 ,
		int* seq1Length ,int** subSeq2Lengths);
void printBestMutants(Mutant* finalScores , int numOfScores);
void fillAllMutants(Mutant* allMutants , int size);
void preprareMutants(Mutant** allMutants , int sizeOfSeq2);
void findMaxMutant(Mutant* rankers , Mutant* bestMutant);



#endif /* FUNCTIONS_H_ */
