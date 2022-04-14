/* Jiwon Kim, 201924445 */
/* 
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */ 
#include <stdio.h>
#include "cachelab.h"

int is_transpose(int M, int N, int A[N][M], int B[M][N]);

char transpose_32_desc[] = "32x32 matrix";
void transpose_32(int M, int N, int A[N][M], int B[M][N]){
	int i,j;
	int row, col;
	int temp = 0;
	int diagonal = 0;
	int size = 8;
	
	for(col=0;col<M;col+=size){
		for(row=0;row<N;row+=size){
			for(i=row;i<(row+size);i++){
				for(j=col;j<(col+size);j++){
					if(i!=j) B[j][i] = A[i][j];
					else{
						diagonal = i;
						temp = A[i][j];
					}
				}
				if(row==col) B[diagonal][diagonal]=temp;
			}
		}
	}
}

char transpose_64_desc[] = "64x64 matrix";
void transpose_64(int M, int N, int A[N][M], int B[M][N]){
	int col,row,i,j;
	int a0,a1,a2,a3,a4,a5,a6,a7;
	
	for(row=0;row<N;row+=8){
		for(col=0;col<M;col+=8){
			for(i=0;i<4;i++){
				a0 = A[row+i][col+0];
				a1 = A[row+i][col+1];
				a2 = A[row+i][col+2];
				a3 = A[row+i][col+3];
				a4 = A[row+i][col+4];
				a5 = A[row+i][col+5];
				a6 = A[row+i][col+6];
				a7 = A[row+i][col+7];
				
				B[col+0][row+i] = a0;
				B[col+1][row+i] = a1;
				B[col+2][row+i] = a2;
				B[col+3][row+i] = a3;
				
				B[col+0][row+i+4] = a4;
				B[col+1][row+i+4] = a5;
				B[col+2][row+i+4] = a6;
				B[col+3][row+i+4] = a7;
			}
			
			for(i=4;i<8;i++){
				a0 = B[col+i-4][row+0+4];
				a1 = B[col+i-4][row+1+4];
				a2 = B[col+i-4][row+2+4];
				a3 = B[col+i-4][row+3+4];
				
				for(j=0;j<4;j++) B[col+i-4][row+j+4] = A[row+j+4][col+i-4];
				
				B[col+i][row+0] = a0;
				B[col+i][row+1] = a1;
				B[col+i][row+2] = a2;
				B[col+i][row+3] = a3;
				
				for(j=4;j<8;j++) B[col+i][row+j] = A[row+j][col+i];
			}
		}
	}	
}

char transpose_default_desc[] = "Any matrix";
void transpose_default(int M, int N, int A[N][M], int B[M][N]){
	int i,j;
	int row, col;
	int temp = 0;
	int diagonal = 0;
	int size = 16;
	
	for(col=0;col<M;col+=size){
		for(row=0;row<N;row+=size){
			for(i=row;(i<N)&&i<(row+size);i++){
				for(j=col;(j<M)&&j<(col+size);j++){
					if(i!=j) B[j][i] = A[i][j];
					else{
						diagonal = i;
						temp = A[i][j];
					}
				}
				if(row==col) B[diagonal][diagonal]=temp;
			}
		}
	}
}

/* 
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded. 
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N])
{
	if(N==32) transpose_32(M,N,A,B);
	else if(N==64) transpose_64(M,N,A,B);
	else transpose_default(M,N,A,B);
}



/* 
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started. 
 */ 

char transpose_test_1_desc[] = "Transpose test";
void transpose_test_1(int M, int N, int A[N][M], int B[M][N])
{
	int block_size, block_row, block_col;
	int row,col;
	int temp = 0, diagonal = 0;
	
	if(N==32){
		block_size = 8;
		for(block_col = 0; block_col < M; block_col += block_size){
			for(block_row = 0; block_row<N; block_row += block_size){
				for(row = block_row; row<block_row + block_size; row++){
					for(col = block_col; col<block_col + block_size; col++){
						if(row!=col) B[col][row] = A[row][col];
						else{
							temp = A[row][col];
							diagonal = row;
						}
					}
					if(block_row == block_col) B[diagonal][diagonal] = temp;
				}
			}
		}
	}
	else if (N==64){
		block_size = 4;
		int a0, a1, a2, a3, a4;
		
		for(col=0; col<M; col+=block_size){
			for(row = 0; row<N; row+=block_size){
				a0 = A[row][col];
				a1 = A[row+1][col];
				a2 = A[row+2][col];
				a3 = A[row+2][col+1];
				a4 = A[row+2][col+2];
				
				B[col+3][row] = A[row][col+3];
				B[col+3][row+1] = A[row+1][col+3];
				B[col+3][row+2] = A[row+2][col+3];
				
				B[col+2][row] = A[row][col+2];
				B[col+2][row+1] = A[row+1][col+2];
				B[col+2][row+2] = a4;
				a4 = A[row+1][col+1];
				
				B[col+1][row] = A[row][col+1];
				B[col+1][row+1] = a4;
				B[col+1][row+2] = a3;
				
				B[col][row] = a0;
				B[col][row+1] = a1;
				B[col][row+2] = a2;
				
				B[col][row+3] = A[row+3][col];
				B[col+1][row+3] = A[row+3][col+1];
				B[col+2][row+3] = A[row+3][col+2];
				a0 = A[row+3][col+3];
				
				B[col+3][row+3] = a0;		
			}
		}
	}
	
	else{
		block_size = 16;
		
		for(block_col = 0; block_col<M; block_col+=block_size){
			for(block_row = 0; block_row<N; block_row+=block_size){
				for(row = block_row; (row<N)&&(row<block_row + block_size); row++){
					for(col = block_col; (col<M)&&(col<block_col + block_size); col++){
						if(row!=col) B[col][row]=A[row][col];
						else{
							temp = A[row][col];
							diagonal = row;
						}	
					}
					if(block_row == block_col) B[diagonal][diagonal] = temp;
				}
			}
		}
	}	
}

char transpose_test_2_desc[] = "Transpose test";
void transpose_test_2(int M, int N, int A[N][M], int B[M][N])
{
	int col, row;
	int i=0, j=0;
	int diagonal=0;
	int temp;
	int size;
	if(M==32){
		size = 8;
		for(row = 0; row<N; row+=size){
			for(col = 0; col<M; col+=size){
				for(i = row; i<row+size; i++){
					for(j=col; j<col+size; j++){
						if(i==j){
							diagonal = i;
							temp = A[i][j];
						}
						else{
							B[j][i]=A[i][j];
						}
					}
					if(row==col) B[diagonal][diagonal]=temp;
				}
			}
		}
	}
	
	else if (M==64){
		size=4;
		for(row = 0; row<N; row+=size){
			for(col = 0; col<M; col+=size){
				for(i = row; i<row+size; i++){
					for(j=col; j<col+size; j++){
						if(i==j){
							diagonal=i;
							temp = A[i][j];
						}
						else{
							B[j][i]=A[i][j];
						}
					}
					
					if(row==col){
						B[diagonal][diagonal]=temp;
					}
				}
			}
		}
	}
	
	else{
		size=16;
		for(row = 0; row<N; row+=size){
			for(col = 0; col<M; col+=size){
				for(i = row; (i<N)&&(i<row+size); i++){
					for(j=col; (j<M)&&(j<col+size); j++){
						if(i==j){
							diagonal=i;
							temp = A[i][j];
						}
						else{
							B[j][i]=A[i][j];
						}
					}
					
					if(row==col){
						B[diagonal][diagonal]=temp;
					}
				}
			}
		}
	}
}
 

/* 
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, tmp;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            tmp = A[i][j];
            B[j][i] = tmp;
        }
    }    

}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc); 

    /* Register any additional transpose functions */
    registerTransFunction(trans, trans_desc); 

}

/* 
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; ++j) {
            if (A[i][j] != B[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}

