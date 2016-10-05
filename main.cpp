#include <iostream>
#include <ctime>
#include <omp.h>
#include <cmath>


#define n 1000
double a[n][n];
double b[n][n];
double transpose[n][n];
double output[n][n];
double output_parallel[n][n];

//parallel
double c_matrix[n][n];
double t_matrix[n][n];


int main() {


//    double temp = 0;
//
#pragma omp parallel for collapse(2)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            b[i][j] = (double) rand() / rand();
            a[i][j] = (double) rand() / rand();
        }
    }

//
//
//// matrix multiplication
//
//    /*
//     * Parallel Implementation
//     */
//
//    double begin = omp_get_wtime();
//
//#pragma omp parallel for reduction(+:temp) collapse(2)
//    for (int k = 0; k < n; ++k) {
//        for (int i = 0; i < n; ++i) {
//            for (int j = 0; j < n; ++j) {
//                temp += a[k][j] * b[j][i];
//            }
//            output_parallel[k][i] = temp;
//            temp = 0;
//        }
//    }
//
//    double end = omp_get_wtime();
//
//    std::cout << "A Output" << std::endl;
//
//    for (int l = 0; l < n; ++l) {
//        for (int i = 0; i < n; ++i) {
//            std::cout << a[l][i] << "\t";
//        }
//        std::cout << std::endl;
//    }
//    std::cout << "B Output" << std::endl;
//
//    for (int l = 0; l < n; ++l) {
//        for (int i = 0; i < n; ++i) {
//            std::cout << b[l][i] << "\t";
//        }
//        std::cout << std::endl;
//    }
//    std::cout << omp_get_num_procs();
//    /**
//     * Sequaential Implementation
//     */
//
//    double begin_seq = omp_get_wtime();
//
//    for (int k = 0; k < n; ++k) {
//        for (int i = 0; i < n; ++i) {
//            for (int j = 0; j < n; ++j) {
//                temp += a[k][j] * b[j][i];
//            }
//            output[k][i] = temp;
//            temp = 0;
//        }
//    }
//
//    double end_seq = omp_get_wtime();
//
//    std::cout << "Sequential Elapsed Time : " << (end_seq - begin_seq) << " secs" << std::endl;
//
////    std::cout << "Sequential Output" << std::endl;
////
////    for (int l = 0; l < n; ++l) {
////        for (int i = 0; i < n; ++i) {
////            std::cout << output[l][i] << "\t";
////        }
////        std::cout << std::endl;
////    }
//
//    std::cout << "Parallel Elapsed Time : " << (end - begin) << " secs" << std::endl;
//
////    transpose B matrix
//#pragma omp parallel for collapse(2)
//    for (int l = 0; l < n; ++l) {
//        for (int i = 0; i < n; ++i) {
//            transpose[l][i] = b[i][l];
//        }
//    }
//
//    double begin_transpose = omp_get_wtime();
//
//#pragma omp parallel for reduction(+:temp) collapse(2)
//    for (int k = 0; k < n; ++k) {
//        for (int i = 0; i < n; ++i) {
//            for (int j = 0; j < n; ++j) {
//                temp += a[k][j] * transpose[i][j];
//            }
//            output_parallel[k][i] = temp;
//            temp = 0;
//        }
//    }
//
//    double end_transpose = omp_get_wtime();
//
//    std::cout << "Parallel(Transpose) Elapsed Time : " << (end_transpose - begin_transpose) << " secs" << std::endl;

    double begin_opt = omp_get_wtime();

#pragma omp parallel num_threads(8)
    {
        int thread_id = omp_get_thread_num();

        int temp_1 = thread_id / 4;
        int temp_3 = thread_id % 2;
        int temp_2 = (thread_id % 4) / 2;
        double tmp = 0;

        int b_j = (n / 2) * temp_3;
        int b_i = (n / 2) * temp_1;
        int a_j = (n / 2) * b_i;
        int a_i = (n / 2) * temp_2;


        if (thread_id < 4) {
            for (int i = a_i; i < a_i + n / 2; ++i) {
                for (int j = b_j; j < b_j + n / 2; ++j) {
                    for (int k = 0; k < n / 2; ++k) {
                        tmp = +a[i][a_j + k] * b[b_i + k][j];
                    }
                    c_matrix[i][j] = tmp;
                    tmp = 0;
                }
            }
        } else {
            for (int i = a_i; i < a_i + n / 2; ++i) {
                for (int j = b_j; j < b_j + n / 2; ++j) {
                    for (int k = 0; k < n / 2; ++k) {
                        tmp = +a[i][a_j + k] * b[b_i + k][j];
                    }
                    t_matrix[i][j] = tmp;
                    tmp = 0;
                }
            }
        }
    }

#pragma omp parallel for collapse(2)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            c_matrix[i][j] = c_matrix[i][j] + t_matrix[i][j];
        }
    }

    double end_opt = omp_get_wtime();

    std::cout << "Optimized (Transpose) Elapsed Time : " << (end_opt - begin_opt) << " secs" << std::endl;

//    std::cout << "Optimized Output" << std::endl;
//
//    for (int l = 0; l < n; ++l) {
//        for (int i = 0; i < n; ++i) {
//            std::cout << c_matrix[l][i] << "\t";
//        }
//        std::cout << std::endl;
//    }

    return 0;
}