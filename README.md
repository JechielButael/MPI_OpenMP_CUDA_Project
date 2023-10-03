# MPI_OpenMP_CUDA_Project
Sequence Alignment – a way to estimate a similarity of two strings of letters - is an important field in bioinformatics1. Sequence is a string of capital letters, for example:
PSHLQYHERTHTGEKPYECHQCGQAFKKCSLLQRHKRTHTGEKPYE

Each letter in the sequence represents DNA, RNA, or protein. Identification of region of similarity of set of Sequences is extremely time consuming. This project deals with a simplified version of Sequence Alignment of two sequences. The purpose of the project is to parallelize the basic algorithm to produce an efficient computation within MPI, OpenMP and CUDA environment.

Alignment Score Definition with pair-wise comparison
1. Similarity of two sequences Seq1 and Seq2 of equal length is defined as follows:
• Two Sequences are places one under another:
PSHLQYHERTHTGEKPYECHQCGQAFKKCSLLQRHKRTHTGEKPYE
HSHLQCHKRTHTGEKPYECNQCGKAFSQHGLLQRHKRTHTGEKPYM
 **** *:***********:***:**.: .***************
 • Each letter from Seq1 is compared with the correspondent letter from Seq2. If these letters are identical the pair is marked with Star sign (*).
 • Otherwise, the additional check is provided. The letters are checked if they both present at least in one of 9 groups called Conservative Groups:
NDEQ  NEQK  STA
MILV  QHRK  NHQK
FYW   HY    MILF  

In case that the pair is found in one of Conservative Group it is marked with Colon sign (:). For example, the pair (E, K) is marked with sign : because they both were found in group NEQK
• If no Conservative Group is found, the pair is checked against 11 Semi-Conservative Groups
SAG     ATV     CSA
SGND    STPA    STNK
NEQHRK  NDEQHK  SNDEQK
HFY     FVLIM

If the pair do presents in one of Semi-Conservative Groups, it is marked with Point sign (.). For example, the pair (K,S) is marked with sign . because they both were found in group STNK
• If the letters in the pair are not equal, do not present both not in Conservative nor in Semi-Conservative groups – the pair is marked with Space sign (‘ ‘).
At the end of the check process the whole Sequence of Signs is obtained. This Sequence is used to estimate the similarity of two sequences – Seq1 and Seq2. For this project following formula is used to estimate the Alignment Score: S = W1*NumberOfStars - W2*NumberOfColons - W3*NumberOfPoints - W4*NumberOfSpaces where Wi are the given weight coefficients (the values of Wi are positive).

2. Similarity of two sequences Seq1 and Seq2 in case that Seq2 is shorter than Seq1, is defined as follows:
• The Sequence Seq2 is places under the Sequence Seq1 with offset n from the start of the Sequence Seq1. The Sequence Seq2 do not allowed to pass behind the end of Seq1. • The letters from Seq1 that do not have a corresponding letter from Seq2 are ignored. • The Alignment Score is calculated according to the pair-wise procedure described above.

For example, Sequence Seq2 is placed at different offsets under Seq1:
PSHLQYHERTHTGEKPYECHQCGQAFKKCSLLQRHKRTHTGEKPYE 
     PYECNQCGKAFSQHGLLQRHKRTHTGEKPYM
      :* .: *:   :     :  :. :  : : 

PSHLQYHERTHTGEKPYECHQCGQAFKKCSLLQRHKRTHTGEKPYE 
               PYECNQCGKAFSQHGLLQRHKRTHTGEKPYM
               ****:***:**.: .***************



Mutant Sequence Definition
For a given Sequence S we define a Mutant Sequence MS(n, k) which is received by removing n-th and k-th letter in S. The indices n and k are positive, letter location is counted from left to right, n < k. For example, for a Sequence S = PSHLQY a set of possible Mutant Sequences (not a full set is displayed):

MS(1, 2) = HLQY
MS(1, 3) = SLQY
MS(1, 4) = SHQY
MS(2, 3) = PLQY
MS(2, 4) = PHQY
MS(3, 6) = PSLQ

Problem Definition
Let Seq1 and Seq2 be a given Sequences of letters.
For all Mutant Sequences of Seq2 find an offset which produce a maximum Alignment Score against Seq1. Among all results, choose the Mutant Sequence MS(n, k) with the best Alignment Score. For example, for Seq1 = PSHLQY and Seq2 = SHQP, the Mutant Sequence MS(3, 4) = SH produce the best Alignment Score at offset = 1

PSHLQY
 SH


Input data and Output Result of the project
You will be provided with the input file input.txt with Weight Coefficients, Sequence Seq1 and few Sequences Seq2. The data is taken from2
The output must be written to file output.txt. For each Sequence Seq2 from input file provide an offset of the Mutant Sequence MS(n, k) with the best Alignment Score and pair (n, k) of this Mutant Sequence
Input File format
First line - weight coefficients W1, W2, W3, W4
Next line – Seq1 (not more than 5000 chars in line)
Next line – the number NS2 of Sequences Seq2 to check against Seq1
Next NS2 lines - Seq2 in each line (not more than 3000 chars in each line)
Output File format
This file contains NS2 lines with offset, n and k found for each Sequence Seq2 from the input file, in order Sequences Seq2 appear in the input file. offset, n and k defines the location of the Mutant Sequence MS(n, k) with the best Alignment Score.

Requirements
• Implement the Simplified Sequence Alignment algorithm explained in the class (see above).
• The input file input.txt initially is known for one machine only. The results must be written to the file output.txt on the same machine
• The computation time of the parallel program must be faster than sequential solution.
• For a given Seq1 and Seq2, the computation stops when the best Alignment Score is achieved.



