#include <iostream>
#include <ctime>
#include <omp.h>
#include <cmath>
#include <emmintrin.h>
#include <immintrin.h>

#define size 1000

double a[size][size];
double b[size][size];
double transpose[size][size];
double output[size][size];

void sequential();

void default_parallel();

void optimized_parallel();

int main() {
/**
 * initializing values for matrices
 */
    // values are being initialized in parallel
#pragma omp parallel for collapse(2)
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            a[i][j] = (double) rand() / rand();
            b[i][j] = (double) rand() / rand();
        }
    }

    sequential();
    default_parallel();
    optimized_parallel();

    return 0;
}

void sequential() {

/**
 * Sequential Code
 */
    double begin_seq = omp_get_wtime(); //start time count for sequential calculation

    for (int k = 0; k < size; ++k) {
        for (int i = 0; i < size; ++i) {
            output[k][i] = 0;
            for (int j = 0; j < size; ++j) {
                output[k][i] += a[k][j] * b[j][i];      // perform sequential matrix calculation
            }
        }
    }

    double end_seq = omp_get_wtime();   //get terminal time count for sequential calculation
    std::cout << "Sequential Elapsed Time : " << (end_seq - begin_seq) << " secs" << std::endl;
}

void default_parallel() {
    /**
     * Default Parallel Code
     */

    double begin = omp_get_wtime(); //start time count for parallel calculation
//openmp parallel block for loops
#pragma omp parallel for
    for (int k = 0; k < size; ++k) {
        for (int i = 0; i < size; i++) {
            output[k][i] = 0;
            for (int j = 0; j < size; j++) {
                output[k][i] += a[k][j] * b[j][i];     //calculate matrix multiplication
            }
        }
    }

    double end = omp_get_wtime();       //terminal time count for parallel calculation
    std::cout << "Default Parallel Elapsed Time : " << (end - begin) << " secs" << std::endl;
}

void optimized_parallel() {

    /**
     * Optimized Parallel Code
     */

    double begin_transpose = omp_get_wtime();           //start time count

    /**
     * Initially get the transpose of the second matrix. As in A * B, get the transpose of B
     */
#pragma omp parallel for num_threads(8) collapse(2)
    for (int l = 0; l < size; ++l) {
        for (int i = 0; i < size; ++i) {
            transpose[l][i] = b[i][l];
        }
    }
/**
 * do matrix calculation via AVX Instruction Parallelsm
 */
    //get the openmp parallel looop block
#pragma omp parallel for
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            __m256d c = _mm256_setzero_pd();        //initialize a __m256 double vector.
            for (int k = 0; k < size; k += 4) {
                /*
                 * here load 4 double values from each matrix and matrix transpose being loaded into SIMD registers.
                 * get the dot product between corresponding matrices.
                 * assign output vector into 'c'
                 *
                 * */

                c = _mm256_add_pd(c, _mm256_mul_pd(_mm256_load_pd(&a[i][k]), _mm256_load_pd(&transpose[j][k])));
            }
            c = _mm256_hadd_pd(c, c);
            // get the upper part(first two elements) from the c vector
            __m128d acc1 = _mm256_extractf128_pd(c, 0);

            // get the lower part(first two elements) from the c vector
            __m128d acc2 = _mm256_extractf128_pd(c, 1);

            //add lower part and upper pert of the c vector and assign the magnitude in output matrix.
            acc1 = _mm_add_sd(acc1, acc2);
            _mm_store_sd(&output[i][j], acc1);

        }
    }
    double end_transpose = omp_get_wtime();     //get terminal time.
    std::cout << "Optimized Parallel Elapsed Time : " << (end_transpose - begin_transpose) << " secs" << std::endl;
}
