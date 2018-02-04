/**
 *  \file parallel-mergesort.cc
 *
 *  \brief Implement your parallel mergesort in this file.
 */
//#include <omp.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "sort.hh"

/*  returns start if T is empty or if x <= T[start] or
    returns mid if x > T[mid -1]
*/
int binarySearch(keytype x, keytype* T, int start, int end)
{
  int low = start;
  int high = start;
  if(start < end + 1)
    high = end + 1;

  int mid = 0;
  while(low < high){
    mid = (low + high) / 2;
    if(x <= T[mid])
      high = mid;
    else
      low = mid + 1;
  }
  return high;
}//binarySearch

/*  merges subarrays of T to form A;
    calls binarySearch
*/
// void seq_merge(keytype* T, int lower1, int upper1, int lower2, int upper2,
//             keytype* A, int lowerOutput)
// {
//   int n1 = upper1 - lower1 + 1;
//   int n2 = upper2 - lower2 + 1;
//
//
// }

void seqMergenew(keytype* T,int p1,int r1,int p2,int r2,keytype* A,int p3)
{
	int n1 = r1-p1+1;
	int n2 = r2-p2+1;
	if(n1<n2)
	{
		int temp=0;
		temp = p1;
		p1 = p2;
		p2 = temp;

		temp=r1;
		r1 = r2;
		r2 = temp;

		temp = n1;
		n1 = n2;
		n2 = temp;

	}
	if(n1 ==0) return;
	else
	{
		int position_a = p1;
		int position_b = p2;
		int position_arr = p3;
		while (position_a<=r1 && position_b<=r2) {
			if(T[position_a] < T[position_b])
			{
				A[position_arr] = T[position_a];
				position_a ++ ;
				position_arr++;
			}
			else
			{
				A[position_arr] = T[position_b];
				position_b ++;
				position_arr++;

			}

		}

		while(position_a<=r1)
		{
			A[position_arr] = T[position_a];
			position_a ++ ;
			position_arr++;
		}

		while(position_b<=r2)
		{
			A[position_arr] = T[position_b];
			position_b ++;
			position_arr++;
		}


	}

}



void merge( keytype* T, int lower1, int upper1, int lower2, int upper2,
            keytype* A, int lowerOutput)
{
  int n1 = upper1 - lower1 + 1;
  int n2 = upper2 - lower2 + 1;
  int lower_bound = 100;
  if(n1<lower_bound && n2<lower_bound)
  {

        seqMergenew(T,lower1,upper1,lower2,upper2,A,lowerOutput);
        return;
  }


  if(n1 < n2){
    int myp = lower2;
    lower2 = lower1;
    lower1 = myp;
    int myr = upper2;
    upper2 = upper1;
    upper1 = myr;
    int myn = n2;
    n2 = n1;
    n1 = myn;
  }
  if (n1 == 0){
    return;
    printf("return\n");
  }
  else {
    int mid1 = (lower1 + upper1) / 2;
    int split2 = binarySearch(T[mid1], T, lower2, upper2);
    int indx_Divide = lowerOutput + (mid1 - lower1) + (split2 - lower2);
    A[indx_Divide] = T[mid1];
	//	#pragma omp parallel num_threads(3)
	//	{
		//#pragma omp single nowait
		//{
		//#pragma omp task
		{
    merge(T, lower1, mid1 - 1, lower2, split2 - 1 , A, lowerOutput);
		merge(T, mid1 + 1, upper1, split2, upper2, A, indx_Divide + 1);

//    seqMergenew(T, lower1, mid1 - 1, lower2, split2 - 1 , A, lowerOutput);
//    seqMergenew(T, mid1 + 1, upper1, split2, upper2, A, indx_Divide + 1);


		}
		//#pragma omp taskwait
	//	}
	//	}
	}
}//merge

/*  recursively orders A; calls merge
*/
void parallelSort(keytype* A, int start, int end, keytype* B, int startOutput)
{
  int n = end - start + 1;
  int mid = 0;
  int notQ = 0;
  int lower_bound = 10;

  if(n == 1)
    B[startOutput] = A[start];
  else {
    keytype* T = newKeys(n);
    mid = (start + end) / 2;
    notQ = mid - start + 1;


    if(n<lower_bound)
    {
      parallelSort(A, start, mid, T, 1);
      parallelSort(A, mid + 1, end, T, notQ + 1);
  		seqMergenew(T,1,notQ,notQ+1,n,B,startOutput);
    }
	//	#pragma omp parallel
	//	{
	//	#pragma omp task
	//	{
  else
  {
    parallelSort(A, start, mid, T, 1);
    parallelSort(A, mid + 1, end, T, notQ + 1);
		//}
	//	#pragma omp taskwait
	//	}
//		#pragma omp parallel
//		{
	//	#pragma omp single nowait
	//	{
		merge(T, 1, notQ, notQ + 1, n, B, startOutput);

		//}
	//	}
  }

	}
}//parallelSort

/* eof */
