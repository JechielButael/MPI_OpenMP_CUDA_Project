/*
 * cuda.h
 *
 *  Created on: 15 Mar 2022
 *      Author: linuxu
 */

#ifndef CUDA_H_
#define CUDA_H_

#include "functions.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


#define MAX_THREADS 1024


void cudaKernal(Mutant* allMutants , char* seq1, char* seq2 , int sizeOfSeq2 ,int seq1Length, float*  weights
                     ,Mutant* currentMutant);

#endif /* CUDA_H_ */
