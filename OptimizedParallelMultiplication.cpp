//#include <iostream>
//#include <omp.h>
//#include <emmintrin.h>
//#include <immintrin.h>
//
//
//#define size 1000
//
//double a[size][size];
//double b[size][size];
//double transpose[size][size];
//double output[size][size];
//
//int main() {
//
//#pragma omp parallel for
//    for (int i = 0; i < size; ++i) {
//        for (int j = 0; j < size; ++j) {
//            a[i][j] = (double) rand() / rand();
//            b[i][j] = (double) rand() / rand();
//        }
//    }
//
//    double begin_transpose = omp_get_wtime();
//
//#pragma omp parallel for num_threads(8) collapse(2)
//    for (int l = 0; l < size; ++l) {
//        for (int i = 0; i < size; ++i) {
//            transpose[l][i] = b[i][l];
//        }
//    }
//
//#pragma omp parallel for
//    for (int i = 0; i < size; i++) {
//        for (int j = 0; j < size; j++) {
//            __m256d c = _mm256_setzero_pd();
//            for (int k = 0; k < size; k += 4) {
//                c = _mm256_add_pd(c, _mm256_mul_pd(_mm256_load_pd(&a[i][k]), _mm256_load_pd(&transpose[j][k])));
//            }
//            c = _mm256_hadd_pd(c, c);
//
//            __m128d acc1 = _mm256_extractf128_pd(c, 0);
//            __m128d acc2 = _mm256_extractf128_pd(c, 1);
//
//            acc1 = _mm_add_sd(acc1, acc2);
//            _mm_store_sd(&output[i][j], acc1);
//
//        }
//    }
//    double end_transpose = omp_get_wtime();
//
//    std::cout << "Optimized Parallel Elapsed Time : " << (end_transpose - begin_transpose) << " secs" << std::endl;
//
//    return 0;
//}