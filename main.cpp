//#include <iostream>
//#include <ctime>
//#include <omp.h>
//#include <cmath>
//
//
//#define n 1000
//
//
//double a[n][n];
//double b[n][n];
//double transpose[n][n];
//double output[n][n];
//double output_parallel[n][n];
//
//
////parallel
//double count_1 = 0;
//double count_2 = 0;
//
//
//int main() {
//
//    double temp = 0;
////
//    for (int l = 0; l < 35; ++l) {
//#pragma omp parallel for collapse(2)
//        for (int i = 0; i < n; ++i) {
//            for (int j = 0; j < n; ++j) {
//                a[i][j] = (double) rand() / rand();
//                b[i][j] = (double) rand() / rand();
//            }
//        }
//
//        double begin = omp_get_wtime();
//
//#pragma omp parallel for
//        for (int k = 0; k < n; ++k) {
//            for (int i = 0; i < n; i++) {
//                output_parallel[k][i] = 0;
//                for (int j = n; j--;) {
//                    output_parallel[k][i] += a[k][n - j] * b[n - j][i];
//                }
//            }
//        }
//
//        double end = omp_get_wtime();
//        count_1 += (end - begin);
//
//
//
//        double begin_seq = omp_get_wtime();
//
//        for (int k = 0; k < n; ++k) {
//            for (int i = 0; i < n; ++i) {
//                output[k][i] = 0;
//                for (int j = 0; j < n; ++j) {
//                    output[k][i] += a[k][j] * b[j][i];
//                }
//            }
//        }
////
//        double end_seq = omp_get_wtime();
//        count_2 += (end_seq - begin_seq);
//    }
//
//            std::cout << "Sequential Elapsed Time : " << (double)count_2/35  << " secs" << std::endl;
//            std::cout << "Parallel Elapsed Time : " << (double)count_1/35 << " secs" << std::endl;
//    exit(0);
////
////
////// matrix multiplication
////
////    /*
////     * Parallel Implementation
////     */
////
////    double begin = omp_get_wtime();
//////
////#pragma omp parallel for
////    for (int k = 0; k < n; ++k) {
////        for (int i = 0; i < n; i++) {
////            output_parallel[k][i] = 0;
////            for (int j = n; j--;) {
////                output_parallel[k][i] += a[k][n - j] * b[n - j][i];
////            }
////
////        }
////    }
//////
////    double end = omp_get_wtime();
////    std::cout << "Parallel Elapsed Time : " << (end - begin) << " secs" << std::endl;
//
//
////    std::cout << "A Output" << std::endl;
////
////    for (int l = 0; l < n; ++l) {
////        for (int i = 0; i < n; ++i) {
////            std::cout << a[l][i] << "\t";
////        }
////        std::cout << std::endl;
////    }
////    std::cout << "B Output" << std::endl;
//
////    for (int l = 0; l < n; ++l) {
////        for (int i = 0; i < n; ++i) {
////            std::cout << b[l][i] << "\t";
////        }
////        std::cout << std::endl;
////    }
////    std::cout << sizeof(a)/(1024*1024)<<std::endl;
////    /**
////     * Sequaential Implementation
////     */
////
//    double begin_seq = omp_get_wtime();
//
//    for (int k = 0; k < n; ++k) {
//        for (int i = 0; i < n; ++i) {
//            output[k][i] = 0;
//            for (int j = 0; j < n; ++j) {
//                output[k][i] += a[k][j] * b[j][i];
//            }
//        }
//    }
////
//    double end_seq = omp_get_wtime();
//
//    std::cout << "Sequential Elapsed Time : " << (end_seq - begin_seq) << " secs" << std::endl;
//
//
////
//////    std::cout << "Sequential Output" << std::endl;
//////
//////    for (int l = 0; l < n; ++l) {
//////        for (int i = 0; i < n; ++i) {
//////            std::cout << output[l][i] << "\t";
//////        }
//////        std::cout << std::endl;
//////    }
////
//
////
//////    transpose B matrix
////    int _temp = 0, temp_1 = 0, _temp_1 = 0, temp_2 = 0, _temp_2 = 0, temp_3 = 0, _temp_3 = 0, temp_4 = 0, _temp_4 = 0;
////    double begin_transpose = omp_get_wtime();
////#pragma omp parallel for num_threads(8) collapse(2)
////    for (int l = 0; l < n; ++l) {
////        for (int i = 0; i < n; ++i) {
////            transpose[l][i] = b[i][l];
////        }
////    }
//////
//////    __m128d *p1, *p2, v1, v2;
//////    double prod[2];
//////    double *p, *ptmp;
//////    p = &output[0][0];
////
////#pragma omp parallel for reduction(+:temp) collapse(2)
////    for (int k = 0; k < n; ++k) {
//////        ptmp = p;
////        for (int i = 0; i < n; i++) {
//////            p1 = (__m128d *) &a[k][0];
//////            p2 = (__m128d *) &transpose[i][0];
//////            v1 = _mm_sub_pd(v1, v1);
////            for (int j = 0; j < n; j ++) {
////                temp += a[k][j] * transpose[i][j];
//////                v2 = _mm_mul_pd(*p1, *p2);
//////                v1 = _mm_add_pd(v1, v2);
//////                p1++;
//////                p2++;
////            }
////            output[k][i] = temp;
////            temp = 0;
//////            _mm_store_pd(prod, v1);
//////            *ptmp += prod[0];
//////            *ptmp += prod[1];
//////            ptmp++;
////        }
//////        p += n;
////    }
////
////    double end_transpose = omp_get_wtime();
////
////    std::cout << "Parallel(Transpose) Elapsed Time : " << (end_transpose - begin_transpose) << " secs" << std::endl;
//
//    //parallel opt
//
//
////    double begin_opt = omp_get_wtime();
////    double tmp_1 = 0;
////#pragma omp parallel for num_threads(8) collapse(2)
////    for (int i = 0; i < n; ++i) {
////        for (int j = i + 1; j < n; ++j) {
////            tmp_1 = b[i][j];
////            b[i][j] = b[j][i];
////            b[j][i] = tmp_1;
////        }
////    }
////
////#pragma omp parallel num_threads(8)
////    {
////        int thread_id = omp_get_thread_num();
////
////        int temp_1 = thread_id / 4;
////        int temp_3 = thread_id % 2;
////        int temp_2 = (thread_id % 4) / 2;
////        double tmp = 0;
////
////        int a_i = (n / 2) * temp_1;
////        int b_j = a_i;
////        int b_i = (n / 2) * temp_2;
////        int a_j = temp_3;
////
//////        printf("Thread:%d %d %d %d %d\n",thread_id, a_i, a_j, b_i, b_j);
////
////        if (thread_id < 4) {
////
////            for (int i = a_i; i < a_i + n / 2; ++i) {
////                for (int j = b_i; j < b_i + n / 2; ++j) {
////                    for (int k = 0; k < n / 2; ++k) {
////                        tmp += a[i][a_j + k] * b[b_i][b_j + k];
////                    }
////                    c_matrix[i][j] = tmp;
////                    tmp = 0;
////                }
////            }
////        } else {
////            for (int i = a_i; i < a_i + n / 2; ++i) {
////                for (int j = b_i; j < b_i + n / 2; ++j) {
////                    for (int k = 0; k < n / 2; ++k) {
////                        tmp += a[i][a_j + k] * b[b_i][b_j + k];
////                    }
////                    t_matrix[i][j] = tmp;
////                    tmp = 0;
////                }
////            }
////        }
////    }
////
////#pragma omp parallel for collapse(2)
////    for (int i = 0; i < n; ++i) {
////        for (int j = 0; j < n; ++j) {
////            c_matrix[i][j] = c_matrix[i][j] + t_matrix[i][j];
////        }
////    }
////
////    double end_opt = omp_get_wtime();
////
////    std::cout << "Optimized Elapsed Time : " << (end_opt - begin_opt) << " secs" << std::endl;
//
////    std::cout << "Optimized Output" << std::endl;
////
////    for (int l = 0; l < n; ++l) {
////        for (int i = 0; i < n; ++i) {
////            std::cout << output[l][i] << "\t";
////        }
////        std::cout << std::endl;
////    }
//
//    //cache optimization
//
////    int Tile = 50;
////
////    double begin_cache = omp_get_wtime();
////
//////#pragma omp parallel for
////    for (int i = 0; i < n; i += Tile) {
////        for (int j = 0; j < n; j += Tile) {
////            for (int k = 0; k < n; k += Tile) {
////                for (int y = i; y < i + Tile; y++) {
////                    for (int x = j; x < j + Tile; x++) {
////                        for (int z = k; z < k + Tile; z++) {
////                            temp += a[y][z] * b[z][x];
////                        }
////                        output_parallel[y][x] = temp;
////                        temp = 0;
////                    }
////                }
////            }
////        }
////    }
////    double end_cache = omp_get_wtime();
////    std::cout << "Optimized cache Elapsed Time : " << (end_cache - begin_cache) << " secs" << std::endl;
//
//    return 0;
//}
