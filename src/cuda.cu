#include "cuda.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <cuda_runtime_api.h> 
#include <cooperative_groups.h>


//CHECK FUNCTION
__host__ void checkStatus(cudaError_t cudaStatus, const char err[])
{
    if(cudaStatus != cudaSuccess)
    {
        printf("%s",err);
        exit(1);
    }
}

__device__ int checkCoservativeGroup(char seq1 , char seq2)
{
    const char* CONSERVATIVE_GROUP[] = {"NDEQ" , "NEQK" , "STA", "FYW",
        "MILV", "MILF" , "QHRK" , "NHQK" , "HY"};
	//strlen(CONSERVATIVE_GROUP) doesn't work, went for 9
    int j ,k;
	for(int i = 0 ; i < 9 ; i++){
        j=0;
        k=0;
        while(CONSERVATIVE_GROUP[i][j] != '\0')
        {
            if(CONSERVATIVE_GROUP[i][j] == seq1)
            {
                 while(CONSERVATIVE_GROUP[i][k] != '\0')
                {
                    if(CONSERVATIVE_GROUP[i][k] == seq2)
                    {
                        return 1;
                    }
                    k++;
                }
            }
            j++;
        }
	}
	return 0;
}




__device__ int checkSemiCoservativeGroup(char seq1 , char seq2)
{
    const char* SEMI_CONSERVATIVE_GROUPS[] = {"SAG" , "ATV", "CSA" , "SGND", "STPA",
        "STNK" , "NEQHRK" , "NDEQHK" , "SNDEQK" , "HFY" , "FVLIM"};
    int j,k;

	for(int i = 0 ; i < 11 ; i++){
        j=0;
        k=0;
        while(SEMI_CONSERVATIVE_GROUPS[i][j] != '\0')
        {
            if(SEMI_CONSERVATIVE_GROUPS[i][j] == seq1)
            {
                while(SEMI_CONSERVATIVE_GROUPS[i][k] != '\0')
                {
                    if(SEMI_CONSERVATIVE_GROUPS[i][k] == seq2)
                    {
                        return 1;
                    }
                    k++;
                }
            }
            j++;
        }
	}
	return 0;
}


__device__ void calcScore(float* weights, Mutant* currentMutant, char* seq1 , 
                char* seq2, int  sizeOfSeq2, int seq1Index)
{
	int spaces = 0;
    int points = 0;
    int colons = 0;
    int stars = 0;


	for(int seq2Index = 0 ; seq2Index < sizeOfSeq2 ; seq2Index++)
	{
		if(seq2Index == currentMutant->n || seq2Index == currentMutant->k){
			continue;
		}
		else if(seq1[seq1Index] == seq2[seq2Index])
		{
			stars++;
			seq1Index++;

		}
		else if(checkCoservativeGroup(seq1[seq1Index],seq2[seq2Index]))
		{
			colons++;
			seq1Index++;

		}
		else if(checkSemiCoservativeGroup(seq1[seq1Index],seq2[seq2Index]))
		{
			points++;
			seq1Index++;

		}
		else{
			spaces++;
			seq1Index++;
		}
	}

	currentMutant->score = (float)(weights[0]*stars) - (float)(weights[1]*colons) - (float)(weights[2]*points) - (float)(weights[3]*spaces);
}



__device__ void getMax(Mutant* mutants,Mutant* best , int numOfMutants){
	Mutant max = mutants[0];
	for (int i = 1; i < numOfMutants; i++) {
		if(mutants[i].score>max.score){
			max = mutants[i];
		}
	}
	best->score = max.score;
	best->k = max.k;
	best->n = max.n;
}



__global__ void cudaCalculateBestScore (float* weights ,char* seq1, char* seq2, Mutant* cudaMutants , int sizeOfSeq2 
                                ,int seq1Length , Mutant* bestMutant , int offset)
{

    //each i represents a mutant
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int sizeOfMutants = (sizeOfSeq2*(sizeOfSeq2-1))/2;
    
    //each mutant's score will be calcualted if he's lower than max num of mutants
    if (i < sizeOfMutants)
        calcScore(weights, &cudaMutants[i] ,seq1 , seq2, sizeOfSeq2, offset);

    //__syncthreads();

    //getting random mutant to calculate the bestMutant out of all the mutants:
    if(i==0)
        getMax(cudaMutants , bestMutant , sizeOfMutants);

    //__syncthreads();
}  



void cudaKernal(Mutant* allMutants , char* seq1, char* seq2 , int sizeOfSeq2 ,int seq1Length, float*  weights
                     ,Mutant* currentMutant)

{
    //variables for cuda
    cudaError_t cudaStatus;
    char* cudaSeq1;
    char* cudaSeq2;
    Mutant* cudaMutants = NULL;
    float* cudaWeights;
    Mutant* bestMutant = NULL;

    //num of mutants
    int numOfMutants = (sizeOfSeq2*(sizeOfSeq2-1))/2;

    //sizes of each cuda variable
    size_t sizeSeq1 = (seq1Length) * sizeof(char);
    size_t sizeSeq2 = (sizeOfSeq2) * sizeof(char);
    size_t sizeWeights = 4 * sizeof(float);
    int sizeMutants = numOfMutants * sizeof(Mutant);

    //num of blocks and num of threads per block
    int threadsPerBlock = MAX_THREADS;
    int blocksPerGrid = (numOfMutants + threadsPerBlock - 1) / threadsPerBlock;

    //copying and allocating seq2 to cuda:
    cudaStatus = cudaMalloc((void**)&cudaSeq2,sizeSeq2);
    checkStatus(cudaStatus,"Cuda Malloc Failed!");
    cudaStatus = cudaMemcpy(cudaSeq2,seq2,sizeSeq2,cudaMemcpyHostToDevice);
    checkStatus(cudaStatus, "Cuda MEMCPY failed!");
    
    //copying  and allocating seq1 to cuda:
	cudaStatus = cudaMalloc((void**)&cudaSeq1,sizeSeq1);
    checkStatus(cudaStatus,"Cuda Malloc Failed!");
    cudaStatus = cudaMemcpy(cudaSeq1,seq1,sizeSeq1,cudaMemcpyHostToDevice);
    checkStatus(cudaStatus, "Cuda MEMCPY failed!");

    //copying  and allocating array of all mutants (that each of them contains n,k) to cuda:
	cudaStatus = cudaMalloc((void**)&cudaMutants,sizeMutants);
    checkStatus(cudaStatus,"Cuda Malloc Failed!");
    cudaStatus = cudaMemcpy(cudaMutants,allMutants,sizeMutants,cudaMemcpyHostToDevice);
    checkStatus(cudaStatus, "Cuda MEMCPY failed!");

    //copying weights to cuda:
	cudaStatus = cudaMalloc((void**)&cudaWeights,sizeWeights);
    checkStatus(cudaStatus,"Cuda Malloc Failed!");
    cudaStatus = cudaMemcpy(cudaWeights,weights,sizeWeights,cudaMemcpyHostToDevice);
    checkStatus(cudaStatus, "Cuda MEMCPY failed!");

    //copying and allocating bestmutant to cuda:
	cudaStatus = cudaMalloc((void**)&bestMutant,sizeof(Mutant));
    checkStatus(cudaStatus,"Cuda Malloc Failed!");
    cudaStatus = cudaMemcpy(bestMutant,currentMutant,sizeof(Mutant),cudaMemcpyHostToDevice);
    checkStatus(cudaStatus, "Cuda MEMCPY failed!");

    //calling global func to find best mutant
    cudaCalculateBestScore<<<blocksPerGrid,threadsPerBlock>>>(cudaWeights ,cudaSeq1,cudaSeq2,cudaMutants ,sizeOfSeq2 
                                    ,seq1Length , bestMutant , currentMutant->offset);
    cudaStatus = cudaDeviceSynchronize();
    checkStatus(cudaStatus, "Cuda Failed!");

    //copying best mutant from cuda to current mutant of given offset in openmp
    cudaStatus = cudaMemcpy(currentMutant,bestMutant,sizeof(Mutant),cudaMemcpyDeviceToHost);
    checkStatus(cudaStatus, "Cuda MEMCPY failed!");


    //Free all allocations from cuda
    cudaStatus = cudaFree(cudaSeq2);
    checkStatus(cudaStatus,"Cuda Free Failed!");
    cudaStatus = cudaFree(cudaSeq1);
    checkStatus(cudaStatus,"Cuda Free Failed!");
    cudaStatus = cudaFree(cudaWeights);
    checkStatus(cudaStatus,"Cuda Free Failed!");
    cudaStatus = cudaFree(cudaMutants);
    checkStatus(cudaStatus,"Cuda Free Failed!");
    cudaStatus = cudaFree(bestMutant);
    checkStatus(cudaStatus,"Cuda Free Failed!");
}


