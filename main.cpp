#include <iostream>
#include <ctime>
#include <omp.h>


#define n 2000
double a[n][n];
double b[n][n];
double output[n][n];
double output_parallel[n][n];

int main() {


    double temp = 0;

//#pragma omp parallel for collapse(2)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            b[i][j] = rand() / rand();
            a[i][j] = rand() / rand();
        }
    }


// matrix multiplication

    clock_t begin = omp_get_wtime();

#pragma omp parallel for reduction(+:temp) collapse(2)
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                temp += a[k][j] * b[j][i];
            }
            output_parallel[k][i] = temp;
            temp = 0;
        }
    }

    clock_t end = omp_get_wtime();

//    std::cout << "Parallel Output" << std::endl;
//
//    for (int l = 0; l < n; ++l) {
//        for (int i = 0; i < n; ++i) {
//            std::cout << output_parallel[l][i] << "\t";
//        }
//        std::cout << std::endl;
//    }

    clock_t begin_seq = omp_get_wtime();

//    for (int k = 0; k < n; ++k) {
//        for (int i = 0; i < n; ++i) {
//            for (int j = 0; j < n; ++j) {
//                temp += a[k][j] * b[j][i];
//            }
//            output[k][i] = temp;
//            temp = 0;
//        }
//    }

    clock_t end_seq = omp_get_wtime();

    std::cout << "Sequential Elapsed Time : " << (end_seq - begin_seq) << " secs" << std::endl;

//    std::cout << "Sequential Output" << std::endl;
//
//    for (int l = 0; l < n; ++l) {
//        for (int i = 0; i < n; ++i) {
//            std::cout << output[l][i] << "\t";
//        }
//        std::cout << std::endl;
//    }

    std::cout << "Parallel Elapsed Time : " << (end - begin) << " secs" << std::endl;
    return 0;
}