#include "../include/for_you_to_do.h"

int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */
    return 3;
  
}

/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/


int mydgetrf(double *A, int *ipiv, int n) 
{
    /* add your code here */
    int i, max_index;
    double max;

    for(i = 0; i < n - 1; i++){
        max_index = i;
        max = fabs(A[i*n + i]);
        int t;
        // find the max
        for(t = i + 1; t < n; t++){
            if(fabs(A[t * n + i]) > max){
                max_index = t;
                max = fabs(A[t * n + i]);
            }
        }

        if(max == 0){
            printf("coefficient matrix is singular");
            return -1;
        }
        else{
            if(max_index != i){
                //swap info
                int temps= ipiv[i];
                ipiv[i] = ipiv[max_index];
                ipiv[max_index] = temps;
                //swap row 
                double trow[n];
                memcpy(trow, A + i * n, n * sizeof(double));
                memcpy(A + i * n, A + max_index * n, n * sizeof(double));
                memcpy(A + max_index * n, trow, n * sizeof(double));
            }

        }
        //update A
        int j;
        for(j = i + 1; j < n;j++){
            A[j * n + i] = A[j * n + i] / A[i * n + i];
            int k;
            for(k =  i + 1; k < n; k++){
                A[j * n + k] = A[j * n + k] - A[j * n + i] * A[i * n + k];
            }
        }
    }

    return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    int i,j;
    double sum;
    /* add your code here */
    if(UPLO == 'L'){
        double y[n];
        y[0] = B[ipiv[0]];
        for(i = 1; i < n; i++ ){
           sum = 0;
           for(j = 0; j < i; j++ ){
               sum += y[j] * A[i*n + j];
           } 
           y[i] = B[ipiv[i]] - sum;
        }

        for(i = 0; i < n; i++){
            B[i] = y[i];
        }
    }
    else if(UPLO == 'U'){
        double x[n];
        int i,j;
        double sum;
        x[n-1] = B[n-1] / A[(n-1) * n + n - 1];
        for(i = n - 2; i >= 0; i-- ){
            sum = 0;
            for(j = i + 1; j < n; j++){
                sum += x[j] * A[i*n + j];
            }
            x[i] = (B[i] - sum) / A[i*n + i];
        }
        for(i = 0; i < n; i++){
            B[i] = x[i];
        }
    }
    return;
}

/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
    int i2,j2,k2,i1,j1,k1;
    register int m;
    int block_size = 3;
    for(i1 = i; i1 < n; i1 += b){
        for(j1 = j;j1 < n; j1 += b){
            for(k1 = k; k1 < k + b; k1 += b){
                for (i2 = i1; i2 < (i1 + b); i2 += block_size){
                    for (j2 = j1; j2 < (j1 + b); j2 += block_size) {

                        //9 register for C
                        register double c00 = C[i2 * n + j2];
                        register double c01 = C[i2 * n + (j2 + 1)];
                        register double c02 = C[i2 * n + (j2 + 2)];

                        register double c10 = C[(i2 + 1) * n + j2];
                        register double c11 = C[(i2 + 1) * n + (j2 + 1)];
                        register double c12 = C[(i2 + 1) * n + (j2 + 2)];

                        register double c20 = C[(i2 + 2) * n + j2];
                        register double c21 = C[(i2 + 2) * n + (j2 + 1)];
                        register double c22 = C[(i2 + 2) * n + (j2 + 2)];

                        //6 reg for A,B
                        register double a00;
                        register double a10;
                        register double a20;
                        register double b00;
                        register double b01;
                        register double b02;

                        for (k2 = k1; k2 < (k1 + b); k2 += block_size){
                            for(m = 0; m < block_size; m++){
                                a00 = A[i2 * n + k2 + m];
                                a10 = A[(i2 + 1)*n + k2 + m];
                                a20 = A[(i2 + 2)*n + k2 + m];
                                b00 = B[(k2 + m) * n + (j2)];
                                b01 = B[(k2 + m) * n + (j2 + 1)];
                                b02 = B[(k2 + m) * n + (j2 + 2)];
                                c00 -= a00 * b00;
                                c01 -= a00 * b01;
                                c02 -= a00 * b02;
                                c10 -= a10 * b00;
                                c11 -= a10 * b01;
                                c12 -= a10 * b02;
                                c20 -= a20 * b00;
                                c21 -= a20 * b01;
                                c22 -= a20 * b02;
                            }
                            

                        }
                        //Write back 
                        C[i2 * n + j2] = c00;
                        C[i2 * n + (j2 + 1)] = c01;
                        C[i2 * n + (j2 + 2)] = c02;
                        C[(i2 + 1) * n + j2] = c10;
                        C[(i2 + 1) * n + (j2 + 1)] = c11;
                        C[(i2 + 1) * n + (j2 + 2)] = c12;
                        C[(i2 + 2) * n + j2] = c20;
                        C[(i2 + 2) * n + (j2 + 1)] = c21;
                        C[(i2 + 2) * n + (j2 + 2)] = c22;

                    }
                }

            }
        }
    }
    return;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{

    /* add your code here */
    int i,j,k,ic,t, max_index;
    double max;
    for(ic = 0; ic < n - 1;ic +=b){
        for(i = ic; i < ic + b ; i++){
            max_index = i;
            max = fabs(A[i * n + i]);
            for(t = i + 1; t < n; t++){
                if(fabs(A[t * n + i]) > max){
                    max_index = t;
                    max = fabs(A[t * n + i]);
                }
            }
            if(max == 0){
                printf("coefficient matrix is single");
                return -1;
            }
            else{
                //swap 
                if(max_index != i){
                    //swap pivoting info
                    int temps= ipiv[i];
                    ipiv[i] = ipiv[max_index];
                    ipiv[max_index] = temps;
                    //printf("switch max:%d, i:%d\n",max_index,i);
                    double trow[n];
                    memcpy(trow, A + i * n, n*sizeof(double));
                    memcpy(A + i * n, A + max_index * n, n*sizeof(double));
                    memcpy(A + max_index * n, trow, n*sizeof(double));
                }

            }
            //update A(ic:end , ic:end) and A(end+1:n , ic:end)
            for(j = i + 1; j <n;j++){
                A[j*n + i] = A[j*n + i] / A[i*n + i];
                //block
                for(k = i + 1; k < ic + b; k++){
                    //printf("A[J,k] -= A[j,i]* A[i,k] i:%d, j:%d, k:%d\n", i, j, k);
                    A[j*n + k] -= A[j*n + i] * A[i*n + k];
                }
            }
        }

        //update A(ic:end, end+1:n)
        register double total;
        //end = ic + b
        for(i = ic; i < ic + b; i++){
            for(j= ic + b;j < n;j++){
                total = 0;
                for(k = ic; k < i; k++){
                    //printf("total += A[i,k]* A[k,j] i:%d, j:%d, k:%d\n", i, j, k);
                    total += A[i*n + k] * A[k*n + j];
                }
                //printf("A[i,j] -= total i:%d, j:%d\n", i, j);
                A[i*n + j] -= total;
            }
        }
        mydgemm(A, A, A,n, ic + b, ic + b, ic, b);
    }

    

    return 0;
}

