//#include <string>
//#include <vector>
//#include <list>
//#include <cmath>
//#include <algorithm>
//#include <map>
//#include <omp.h>
//#include <set>
//#include <iostream>
//#include <pmmintrin.h>
//
////#include <stdlib.h>
////#include <cstdio>
////#include <cstdlib>
////#include <sys/time.h>
//
//using namespace std;
//
////#define VIT(i, v) for (i = 0; i < v.size(); i++)
////#define IT(it, ds) for (it = ds.begin(); it != ds.end(); it++)
////#define FUP(i, n) for (i = 0; i < n; i++)
//
//typedef vector <int> IVec;
//typedef vector <double> DVec;
//typedef vector <string> SVec;
//
//void usage(string s)
//{
//    fprintf(stderr, "usage: mm r1 c1/r2 c2 blocksize seed(-1=time(0)) print(y|n)\n");
//    if (s != "") fprintf(stderr, "%s\n", s.c_str());
//    exit(1);
//}
//
//class MM {
//public:
//    vector <double> M1;
//    vector <double> M2;
//    vector <double> P;
//    int r1, c1, c2;
//    int Print;
//    int Blocksize;
//    void MultiplyBlock(int row1, int col1, int col2);
//    void Multiply();
//    void PrintAll();
//};
//
//void MM::MultiplyBlock(int row1, int col1, int col2)
//{
//    double *m1top, *m1;
//    double *m2top, *m2base, *m2;
//    double *p, *ptmp, *ptmp2;
//    double prod[2];
//        int tmp, tmptop;
//    __m128d *p1, *p2, *p3, v1, v2, v3, v4;
//
//    m1 = &(M1[row1 * c1 + col1]);
//    m1top = m1 + (c1 * Blocksize);
//    if (m1top > &(M1[r1 * c1])) m1top = &(M1[r1 * c1]);
//
//    m2base = &(M2[col2 * c1 + col1]);
//    m2top = m2base + (c1 * Blocksize);
//    if (m2top > &(M2[c2 * c1])) m2top = &(M2[c2 * c1]);
//
//    p = &(P[row1 * c2 + col2]);
//
//    tmptop = col1 + Blocksize;
//    if (tmptop > c1) tmptop = c1;
//
//    for ( ; m1 < m1top; m1 += (c1*2)) {
//        ptmp = p;
//        ptmp2 = p+c2;
//        for (m2 = m2base; m2 < m2top; m2 += c1) {
//            p1 = (__m128d *) m1;
//            p2 = (__m128d *) m2;
//            p3 = (__m128d *) (m1+c1);
//            v1 = _mm_sub_pd(v1, v1);   /* Zero v1 */
//            v3 = _mm_sub_pd(v3, v3);   /* Zero v3 */
//
//            for (tmp = col1; tmp < tmptop; tmp += 2) {
//                v2 = _mm_mul_pd(*p1, *p2);
//                v4 = _mm_mul_pd(*p3, *p2);
//                v1 = _mm_add_pd(v1, v2);
//                v3 = _mm_add_pd(v3, v4);
//                p1++;
//                p2++;
//                p3++;
//            }
//            _mm_store_pd(prod, v1);
//            *ptmp += prod[0];
//            *ptmp += prod[1];
//            _mm_store_pd(prod, v3);
//            *ptmp2 += prod[0];
//            *ptmp2 += prod[1];
//            ptmp++;
//            ptmp2++;
//        }
//        p += (c2*2);
//    }
//}
//
//void MM::Multiply()
//{
//    int row1, col1, col2;
//
//    for (row1 = 0; row1 < r1; row1 += Blocksize) {
//        for (col2 = 0; col2 < c2; col2 += Blocksize) {
//            for (col1 = 0; col1 < c1; col1 += Blocksize) {
//                MultiplyBlock(row1, col1, col2);
//            }
//        }
//    }
//}
//
//void MM::PrintAll()
//{
//    int i, j;
//
//    printf("M1: %d x %d\n\n", r1, c1);
//    for (i = 0; i < r1; i++) {
//        for (j = 0; j < c1; j++) printf(" %6.4lf", M1[i*c1+j]);
//        printf("\n");
//    }
//    printf("\n");
//
//    printf("M2: %d x %d\n\n", c1, c2);
//    for (i = 0; i < c1; i++) {
//        for (j = 0; j < c2; j++) printf(" %6.4lf", M2[j*c1+i]);
//        printf("\n");
//    }
//    printf("\n");
//
//    printf("P: %d x %d\n\n", r1, c2);
//    for (i = 0; i < r1; i++) {
//        for (j = 0; j < c2; j++) printf(" %6.4lf", P[i*c2+j]);
//        printf("\n");
//    }
//}
//
//int main(int argc, char **argv)
//{
//    MM *M;
//    int r1=1000, c1=1000, c2=1000, i, j, bs=128;
//    string s;
//    long seed;
//    double t0, t1;
//    struct timeval tv;
//
////    if (argc != 7) usage("");
//    M = new MM;
//
////    if (sscanf(argv[1], "%d", &r1) == 0 || r1 <= 0) usage("Bad r1");
////    if (sscanf(argv[2], "%d", &c1) == 0 || c1 <= 0) usage("Bad c1/r2");
////    if (sscanf(argv[3], "%d", &c2) == 0 || c2 <= 0) usage("Bad c2");
////    if (sscanf(argv[4], "%d", &bs) == 0 || bs <= 0) usage("Bad bs");
////
////    if (r1 % 2 == 1) { fprintf(stderr, "This only works on even matrix dimensions.\n"); exit(0); }
////    if (c1 % 2 == 1) { fprintf(stderr, "This only works on even matrix dimensions.\n"); exit(0); }
////    if (c2 % 2 == 1) { fprintf(stderr, "This only works on even matrix dimensions.\n"); exit(0); }
////    if (bs % 2 == 1) { fprintf(stderr, "This only works on even block sizes.\n"); exit(0); }
//
//    M->r1 = r1;
//    M->c1 = c1;
//    M->c2 = c2;
//    M->Blocksize = bs;
//
//    seed = 0;
////    sscanf(argv[5], "%ld", &seed);
////    if (seed == -1) seed = time(0);
//    srand48(seed);
//
////    s = argv[6];
////    s = "y";
////    if (s == "y") {
////        M->Print = 1;
////    } else if (s == "n") {
////        M->Print = 0;
////    } else usage("Bad print");
//
//    M->M1.resize(r1*c1);
//    M->M2.resize(c1*c2);
//    M->P.resize(r1*c2, 0);
//
//    if (random) {
//        for(i = 0; i < r1*c1; i++) M->M1[i] = drand48()*2.0;
//        for(i = 0; i < c1; i++) {
//            for(j = 0; j < c2; j++) M->M2[j*c1+i] = drand48()*2.0;
//        }
//    }
//
//    double begin_seq = omp_get_wtime();
//
//    M->Multiply();
//
//    double end_seq = omp_get_wtime();
////
//    std::cout << "Sequential Elapsed Time : " << (end_seq - begin_seq) << " secs" << std::endl;
////    if (M->Print) M->PrintAll();
//
////    printf("Time: %.4lf\n", t1-t0);
//    exit(0);
//}
