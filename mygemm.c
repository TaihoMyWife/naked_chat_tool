#include "mygemm.h"

/**
 * 
 * Implement all functions here in this file.
 * Do NOT change input parameters and return type.
 * 
 **/

//Register Reuse part 1
void dgemm0(const double* A, const double* B, double* C, const int n)
{
    int i, j, k;
    for (i = 0; i < n; i++) { 
        for (j = 0; j < n; j++) { 
            for (k = 0; k < n; k++) { 
                C[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }
}

void dgemm1(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            register double sum = C[i * n + j]; 
            for (k = 0; k < n; k++) {
                sum += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = sum;
        }
    }
}
//Register Reuse part 1 End

//Register Reuse part 2
void dgemm2(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (i = 0; i < n; i += 2) {
        for (j = 0; j < n; j += 2) {
            register double c00 = C[i * n + j];
            register double c01 = C[i * n + j + 1];
            register double c10 = C[(i + 1) * n + j];
            register double c11 = C[(i + 1) * n + j  + 1];
            for (k = 0; k < n; k += 2) {
                register double a00 = A[i * n + k];
                register double a10 = A[(i + 1) * n + k]; 
                register double b00 = B[k * n + j];
                register double b01 = B[k * n + j + 1];

                c00 += a00 * b00; c01 += a00 * b01;
                c10 += a10 * b00; c11 += a10 * b01;

                a00 = A[i * n + k + 1];       //a01
                a10 = A[(i + 1) * n + k + 1]; //a11
                b00 = B[(k + 1) * n + j];     //b10
                b01 = B[(k + 1) * n + j + 1]; //b11

                c00 += a00 * b00; c01 += a00 * b01;
                c10 += a10 * b00; c11 += a10 * b01;
            }

            C[i * n + j] = c00;
            C[i * n + j + 1] = c01;
            C[(i + 1) * n + j] = c10;
            C[(i + 1) * n + j + 1] = c11;
        }
    }
}
//Register Reuse part 2 End

//Register Reuse part 3
void dgemm3(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (i = 0; i < n; i += 3) {
        for (j = 0; j < n; j += 3) {
            register double c00 = C[i * n + j];
            register double c01 = C[i * n + j + 1];
            register double c02 = C[i * n + j + 2];
            register double c10 = C[i * n + j + n];
            register double c11 = C[i * n + j + n + 1];
            register double c12 = C[i * n + j + n + 2];
            register double c20 = C[i * n + j + 2 * n];
            register double c21 = C[i * n + j + 2 * n + 1];
            register double c22 = C[i * n + j + 2 * n + 2];
            for (k = 0; k < n; k += 3) {
                register double a00 = A[i * n + k];
                register double a10 = A[i * n + k + n]; 
                register double a20 = A[i * n + k + 2 * n];
                register double b00 = B[k * n + j];
                register double b01 = B[k * n + j + 1];
                register double b02 = B[k * n + j + 2];
                c00 += a00 * b00; c01 += a00 * b01; c02 += a00 * b02;
                c10 += a10 * b00; c11 += a10 * b01; c12 += a10 * b02;
                c20 += a20 * b00; c21 += a20 * b01; c22 += a20 * b02;

                a00 = A[i * n + k + 1];                 //a01
                a10 = A[(i + 1) * n + k + 1];           //a11
                a20 = A[(i + 2) * n + k + 1];           //a21
                b00 = B[(k + 1) * n + j];               //b10
                b01 = B[(k + 1) * n + j + 1];           //b11
                b02 = B[(k + 1) * n + j + 2];           //b12

                c00 += a00 * b00; c01 += a00 * b01; c02 += a00 * b02;
                c10 += a10 * b00; c11 += a10 * b01; c12 += a10 * b02;
                c20 += a20 * b00; c21 += a20 * b01; c22 += a20 * b02;

                a00 = A[i * n + k + 2];                 //a02
                a10 = A[(i + 1) * n + k + 2];           //a12
                a20 = A[(i + 2) * n + k + 2];           //a22
                b00 = B[(k + 2) * n + j];               //b20
                b01 = B[(k + 2) * n + j + 1];           //b21
                b02 = B[(k + 2) * n + j + 2];           //b22

                c00 += a00 * b00; c01 += a00 * b01; c02 += a00 * b02;
                c10 += a10 * b00; c11 += a10 * b01; c12 += a10 * b02;
                c20 += a20 * b00; c21 += a20 * b01; c22 += a20 * b02;
            }

            C[i * n + j] = c00;
            C[i * n + j + 1] = c01;
            C[i * n + j + 2] = c02;
            C[(i + 1) * n + j] = c10;
            C[(i + 1) * n + j + 1] = c11;
            C[(i + 1) * n + j + 2] = c12;
            C[(i + 2) * n + j] = c20;
            C[(i + 2) * n + j + 1] = c21;
            C[(i + 2) * n + j + 2] = c22;
        }
    }

}
//Register Reuse part 3 End

//Cache Reuse part 3
void ijk(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (i = 0; i < n; i++) { 
        for (j = 0; j < n; j++) { 
            register double sum = C[i * n + j]; 
            for (k = 0; k < n; k++) { 
                sum += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = sum;
        }
    }
}

void bijk(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i, j, k, i1, j1, k1;
    for (i = 0; i < n; i += b) { 
        for (j = 0; j < n; j += b) { 
            for (k = 0; k < n; k += b) { 
                /* B x B */
                for (i1 = i; i1 < i + b; i1++) { 
                    for (j1 = j; j1 < j + b; j1++) { 
                        register double sum = C[i1 * n + j1];
                        for (k1 = k; k1 < k + b; k1++) {
                            sum += A[i1 * n + k1] * B[k1 * n + j1];
                        }
                        C[i1 * n + j1] = sum;
                    }
                }
            }
        }
    }
}

void jik(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (j = 0; j < n; j++) { 
        for (i = 0; i < n; i++) { 
            register double sum = C[i * n + j]; 
            for (k = 0; k < n; k++) {
                sum += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = sum;
        }
    }
}

void bjik(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i, j, k, i1, j1, k1;
    for (j = 0; j < n; j += b) { 
        for (i = 0; i < n; i += b) { 
            for (k = 0; k < n; k += b) { 
                /* B x B  */
                for (i1 = i; i1 < i + b; i1++) { 
                    for (j1 = j; j1 < j + b; j1++) { 
                        register double sum = C[i1 * n + j1];
                        for (k1 = k; k1 < k + b; k1++) {
                            sum += A[i1 * n + k1] * B[k1 * n + j1];
                        }
                        C[i1 * n + j1] = sum;
                    }
                }
            }
        }
    }
}

void kij(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (k = 0; k < n; k++) { 
        for (i = 0; i < n; i++) { 
            register double T = A[i * n + k];
            for (j = 0; j < n; j++) {
                C[i * n + j] += T * B[k * n + j];
            }
        }
    }
}

void bkij(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i, j, k, i1, j1, k1;
    for (k = 0; k < n; k += b) { 
        for (i = 0; i < n; i += b) { 
            for (j = 0; j < n; j += b) { 
                /* B x B */
                for (k1 = k; k1 < k + b; k1++) { 
                    for (i1 = i; i1 < i + b; i1++) { 
                        register double T = A[i1 * n + k1];
                        for (j1 = j; j1 < j + b ; j1++) {
                            C[i1 * n + j1] += T * B[k1 * n + j1];
                        }
                    }
                }
            }
        }
    }
}


void ikj(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (i=0; i<n; i++) { 
        for (k=0; k<n; k++) { 
            register double T = A[i * n + k]; 
            for (j=0; j<n; j++) {
                C[i * n + j] += T * B[k * n + j];
            }
        }
    }
}

void bikj(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i, j, k, i1, j1, k1;
    for (i = 0; i < n; i += b) { 
        for (k = 0; k < n; k += b) { 
            for (j = 0; j < n; j += b) { 
                /* B x B */
                for (i1 = i; i1 < i + b; i1++) { 
                    for (k1 = k; k1 < k + b; k1++) {
                        register double T = A[i1 * n + k1];
                        for (j1 = j; j1 < j + b ; j1++) {
                            C[i1 * n + j1] += T * B[k1 * n + j1];
                        } 
                    }
                }
            }
        }
    }
}

void jki(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (j = 0; j < n; j++) { 
        for (k = 0; k < n; k++) { 
            register double T = B[k * n + j];
            for (i = 0; i < n; i++) {
                C[i * n + j] += A[i * n + k] * T;
            }
        }
    }
}

void bjki(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i, j, k, i1, j1, k1;
    for (j = 0; j < n; j += b) { 
        for (k = 0; k < n; k += b) { 
            for (i = 0; i < n; i += b) { 
                /* B x B  */
                for (j1 = j; j1 < j + b; j1++) { 
                    for (k1 = k; k1 < k + b; k1++) {
                        register double T = B[k1 * n + j1];
                        for (i1 = i; i1 < i + b; i1++) {
                             C[i1 * n + j1] += A[i1 * n + k1] * T;
                        }
                    }
                }
            }
        }
    }
}

void kji(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (k=0; k<n; k++) { 
        for (j = 0; j < n; j++) { 
            register double T = B[k * n + j];
            for (i = 0; i < n; i++) {
                C[i * n + j] += A[i * n + k] * T;
            }
        }
    }

}

void bkji(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i, j, k, i1, j1, k1;
    for (k = 0; k < n; k += b) { 
        for (j = 0; j < n; j += b) { 
            for (i = 0; i < n; i += b) { 
                /* B x B */
                for (k1 = k; k1 < k + b; k1++) { 
                    for (j1 = j; j1 < j + b; j1++) {
                        register double T = B[k1 * n + j1];
                        for (i1 = i; i1 < i + b; i1++) {
                             C[i1 * n + j1] += A[i1 * n + k1] * T;
                        } 
                    }
                }
            }
        }
    }
}
//Cache Reuse part 3 End 

//Cache Reuse part 4
void optimal(const double* A, const double* B, double *C, const int n, const int b)
{
    int i, j, k, i1, j1, k1;
    for (k = 0; k < n; k += b) { 
        for (i = 0; i < n; i += b) { 
            for (j=0; j<n; j+=b) { 
                /* B x B */
                for (i1 = i; i1 < i + b; i1 += 2) { 
                    for (j1=j; j1<j+b; j1+=2) { 
                        /* 2 X 2  register optimize */
                        register double c00 = C[i1 * n + j1];
                        register double c01 = C[i1 * n + j1 + 1];
                        register double c10 = C[i1 * n + j1 + n];
                        register double c11 = C[i1 * n + j1 + n + 1];
                        for (k1 = k; k1 < k + b; k1 += 2) {
                            register double a00 = A[i1 * n + k1];
                            register double a10 = A[i1 * n + k1 + n]; 
                            register double b00 = B[k1 * n + j1];
                            register double b01 = B[k1 * n + j1 + 1];

                            c00 += a00 * b00; c01 += a00 * b01;
                            c10 += a10 * b00; c11 += a10 * b01;

                            a00 = A[i1 * n + k1 + 1];     
                            a10 = A[i1 * n + k1 + n + 1]; 
                            b00 = B[k1 * n + j1 + n];     
                            b01 = B[k1 * n + j1 + 1 + n]; 

                            c00 += a00 * b00; c01 += a00 * b01;
                            c10 += a10 * b00; c11 += a10 * b01;
                        }
                        C[i1 * n + j1] = c00;
                        C[i1 * n + j1 + 1] = c01;
                        C[i1 * n + j1 + n] = c10;
                        C[i1 * n + j1 + n + 1] = c11;
                    }
                }
            }
        }
    }
    
}
void optimal1(const double* A, const double* B, double *C, const int n, const int b){
    int i, j, k, i1, j1, k1;    
    for (k = 0; k < n; k += b) { 
        for (i = 0; i < n; i += b) {            
            for (j=0; j<n; j+=b) {      /* B x B */                
                for (i1 = i; i1 < i + b; i1 += 2) {
                    for (j1 = j; j1 < j + b; j1 +=2) { 
                        register double c00 = C[i1 * n + j1];
                        register double c01 = C[i1 * n + j1 + 1];
                        register double c10 = C[i1 * n + j1 + n];
                        register double c11 = C[i1 * n + j1 + n + 1];
                        for (k1 = k; k1 < k + b; k1 += 2) { 
                            register double a00 = A[i1 * n + k1];
                            register double a10 = A[i1 * n + k1 + n];
                            register double b00 = B[k1 * n + j1]; 
                            register double b01 = B[k1 * n + j1 + 1];
                            register double a01 = A[i1 * n + k1 + 1];     //a01
                            register double a11 = A[i1 * n + k1 + n + 1]; //a11
                            register double b10 = B[k1 * n + j1 + n];     //b10
                            register double b11 = B[k1 * n + j1 + 1 + n]; //b11
                            double m1=(a00+a11)*(b00+b11);
                            double m2=(a10+a11)*b00;
                            double m3=(b01-b11)*a00;
                            double m4=(b10-b00)*a11;
                            double m5=(a00+a01)*b11;
                            double m6=(a10-a00)*(b00+b01);
                            double m7=(a01-a11)*(b10+b11);
                            c00 += m1+m4-m5+m7;
                            c01 += m3+m5;
                            c10 += m2+m4;
                            c11 += m1-m2+m3+m6;
                        }   
                        C[i1 * n + j1] = c00;
                        C[i1 * n + j1 + 1] = c01;
                        C[i1 * n + j1 + n] = c10;                        
                        C[i1 * n + j1 + n + 1] = c11;                    
                    }                
                 }            
             }        
       }    
    }
}
