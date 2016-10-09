//#include <iostream>
//#include <ctime>
//#include <omp.h>
//#include <cmath>
//
//#define n 1000
//
//double a[n][n];
//double b[n][n];
//double transpose[n][n];
//double output[n][n];
//
//int main() {
//
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
//                output[k][i] = 0;
//                for (int j = n; j--;) {
//                    output[k][i] += a[k][n - j] * b[n - j][i];
//                }
//            }
//        }
//
//        double end = omp_get_wtime();
//
//    std::cout << "Parallel Elapsed Time : " << (end - begin) << " secs" << std::endl;
//
//    return 0;
//}
