//#include <iostream>
//#include <ctime>
//#include <omp.h>
//
//
//
//#define n 1000
//
//
//double a[n][n];
//double b[n][n];
//
//double output[n][n];
//
//int main() {
//
//#pragma omp parallel for
//    for (int i = 0; i < n; ++i) {
//        for (int j = 0; j < n; ++j) {
//            a[i][j] = (double) rand() / rand();
//            b[i][j] = (double) rand() / rand();
//        }
//    }
//
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
//
//    double end_seq = omp_get_wtime();
//    std::cout << "Sequential Elapsed Time : " << (end_seq - begin_seq) << " secs" << std::endl;
//
//    return 0;
//}
