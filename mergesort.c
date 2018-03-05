#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include<omp.h>

#include <x86intrin.h>

#define NBEXPERIMENTS   7

static long long unsigned int experiments [NBEXPERIMENTS] ;

static   unsigned int N ;

typedef  int  *array_int ;

static array_int X ;

void init_array (array_int T)
{
  register int i ;

  for (i = 0 ; i < N ; i++)
    {
      T [i] = N - i ;
    }
}

void print_array (array_int T, int N)
{
  register int i ;

  for (i = 0 ; i < N ; i++)
    {
      printf ("%d ", T[i]) ;
    }
  printf ("\n") ;
}

int is_sorted (array_int T)
{
  register int i ;
  
  for (i = 1 ; i < N ; i++)
    {
        /* test designed specifically for our usecase */
        if (T[i-1] +1  != T [i] )
            return 0 ;
    }
  return 1 ;
}

long long unsigned int average (long long unsigned int *exps)
{
  unsigned int i ;
  long long unsigned int s = 0 ;

  for (i = 2; i < (NBEXPERIMENTS-2); i++)
    {
      s = s + exps [i] ;
    }

  return s / (NBEXPERIMENTS-2) ;
}

void merge (int *T, int size) 
{

    int L[size/2], R[size/2];
    register int i = 0;
    register int j = 0;
    register int k = 0;

    for (i = 0; i < size/2 ; i++) {
        L[i] = T[i];
        R[i] = T[(size/2) + i];
    }

    i = 0;
    j = 0;

    while (i < size/2  && j < size/2) {
      if (L[i] <= R[j]) {
            T[k] = L[i];
            i++;
        } else {
            T[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < size/2) {
        T[k] = L[i];
        i++;
        k++;
    }
 
    while (j < size/2) {
        T[k] = R[j];
        j++;
        k++;
    }
} 
 

void merge_sort (int *T, int size)
{
    /* sequential version of the merge sort algorithm */

    if (size >= 2) {
        merge_sort(T , size/2);
        merge_sort(T + size/2, size/2);
        merge(T, size);
    }
}


void parallel_merge_sort_tasks (int *T, const int size)
{
    /* parallel version of the merge sort algorithm */

    if (size > 1){
        #pragma omp parallel
        {
	    #pragma omp single
	    {
            	#pragma omp task
            	parallel_merge_sort_tasks(T, size/2);
            	#pragma omp task
            	parallel_merge_sort_tasks(T + size/2, size/2); 
	    }
        }
	merge(T, size);
    } 
}

void parallel_merge_sort (int *T, const int size, const int number_of_threads)
{
    /* parallel version of the merge sort algorithm */

    if (number_of_threads == 1){
        merge_sort(T, size);
    } else {
        #pragma omp parallel
        {
	    #pragma omp single
	    {
            	#pragma omp task
            	parallel_merge_sort(T, size/2, number_of_threads/2);
            	#pragma omp task
            	parallel_merge_sort(T + size/2, size/2, number_of_threads/2); 
	    }
        }
        merge(T, size);
    } 
}


int main (int argc, char **argv)
{
  unsigned long long int start, end, residu ;
  unsigned long long int av ;
  unsigned int exp ;

    if (argc != 2)
    {
      fprintf (stderr, "mergesort N \n") ;
      exit (-1) ;
    }

  N = 1 << (atoi(argv[1])) ;
  X = (int *) malloc (N * sizeof(int)) ;

  printf("--> Sorting an array of size %u\n",N);
  
  start = _rdtsc () ;
  end   = _rdtsc () ;
  residu = end - start ; 

  // print_array (X) ;

  printf("sequential sorting ...\n");


    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
      init_array (X) ;
      start = _rdtsc () ;

               merge_sort (X, N) ;
     
      end = _rdtsc () ;
      experiments [exp] = end - start ;
      
      if (! is_sorted (X))
	{
            fprintf(stderr, "ERROR: the array is not properly sorted\n") ;
            exit (-1) ;
	}      
    }
  
  av = average (experiments) ;  
  printf ("\n merge sort serial\t\t %Ld cycles\n\n", av-residu) ;

  printf("parallel sorting...\n");
  
  for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
      init_array (X) ;
      
      start = _rdtsc () ;

           parallel_merge_sort_tasks (X, N) ;
     
      end = _rdtsc () ;
      experiments [exp] = end - start ;

      if (! is_sorted (X))
	{
            print_array(X, N);
            fprintf(stderr, "ERROR: the array is not properly sorted\n") ;
            exit (-1) ;
	}      
    }
  
  av = average (experiments) ;
  printf ("\n merge sort parallel with tasks\t %Ld cycles\n\n", av-residu) ;
	
  printf("parallel sorting...\n");
  
  for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
      init_array (X) ;
      
      start = _rdtsc () ;

           parallel_merge_sort (X, N, omp_get_max_threads()) ;
     
      end = _rdtsc () ;
      experiments [exp] = end - start ;

      if (! is_sorted (X))
	{
            print_array(X, N);
            fprintf(stderr, "ERROR: the array is not properly sorted\n") ;
            exit (-1) ;
	}      
    }
  
  av = average (experiments) ;
  printf ("\n merge sort parallel with tasks and seq\t %Ld cycles\n\n", av-residu) ;

}
