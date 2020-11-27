#include "matrix.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// Include SSE intrinsics
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#include <x86intrin.h>
#endif

/* Below are some intel intrinsics that might be useful
 * void _mm256_storeu_pd (double * mem_addr, __m256d a)
 * __m256d _mm256_set1_pd (double a)
 * __m256d _mm256_set_pd (double e3, double e2, double e1, double e0)
 * __m256d _mm256_loadu_pd (double const * mem_addr)
 * __m256d _mm256_add_pd (__m256d a, __m256d b)
 * __m256d _mm256_sub_pd (__m256d a, __m256d b)
 * __m256d _mm256_fmadd_pd (__m256d a, __m256d b, __m256d c)
 * __m256d _mm256_mul_pd (__m256d a, __m256d b)
 * __m256d _mm256_cmp_pd (__m256d a, __m256d b, const int imm8)
 * __m256d _mm256_and_pd (__m256d a, __m256d b)
 * __m256d _mm256_max_pd (__m256d a, __m256d b)
*/

/*
 * Generates a random double between `low` and `high`.
 */
double rand_double(double low, double high) {
    double range = (high - low);
    double div = RAND_MAX / range;
    return low + (rand() / div);
}

/*
 * Generates a random matrix with `seed`.
 */
void rand_matrix(matrix *result, unsigned int seed, double low, double high) {
    srand(seed);
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            set(result, i, j, rand_double(low, high));
        }
    }
}

/*
 * Allocate space for a matrix struct pointed to by the double pointer mat with
 * `rows` rows and `cols` columns. You should also allocate memory for the data array
 * and initialize all entries to be zeros. Remember to set all fieds of the matrix struct.
 * `parent` should be set to NULL to indicate that this matrix is not a slice.
 * You should return -1 if either `rows` or `cols` or both have invalid values, or if any
 * call to allocate memory in this function fails. If you don't set python error messages here upon
 * failure, then remember to set it in numc.c.
 * Return 0 upon success and non-zero upon failure.
 */
int allocate_matrix(matrix **mat, int rows, int cols) {
    /* TODO: YOUR CODE HERE */
    
    // add Runtime Error

    if (rows <= 0 || cols <= 0) {
        //throw ValueError
        return -1;
    }

    *mat = (matrix *) malloc(sizeof(matrix));
    
    if (mat == NULL) {
        //throw TypeError?
        return -1;
    }

    (*mat)->rows = rows;
    (*mat)->cols = cols;
    //(*mat)->is_1d = 0; //maybe fix
    (*mat)->parent = NULL;
    (*mat)->ref_cnt = 1; //maybe fix

    if (rows == 1 || cols == 1) {
        (*mat)->is_1d = 1;
    } else {
        (*mat)->is_1d = 0;
    }

    int i;
    double **outer = (double **)calloc(rows, sizeof(double*));
    for (i = 0; i < rows; i++) {
    	double *inner = (double *)calloc(cols, sizeof(double));
	outer[i] = inner;
    }
    
    //double *newCols = (double *) calloc(cols, sizeof(double));
    //double **newData = (double **) calloc(rows, sizeof(newCols));

    if (outer == NULL) {
        return -1;
    }

    (*mat)->data = outer;

    return 0;

}

/*
 * Allocate space for a matrix struct pointed to by `mat` with `rows` rows and `cols` columns.
 * This is equivalent to setting the new matrix to be
 * from[row_offset:row_offset + rows, col_offset:col_offset + cols]
 * If you don't set python error messages here upon failure, then remember to set it in numc.c.
 * Return 0 upon success and non-zero upon failure.
 */
int allocate_matrix_ref(matrix **mat, matrix *from, int row_offset, int col_offset,
                        int rows, int cols) {

    if (rows <= 0 || cols <= 0 || rows > from->rows || cols > from->cols || row_offset < 0 || col_offset < 0) {
        return -1;
    }
    
    *mat = (matrix *) malloc(sizeof(matrix));
    
    if (mat == NULL) {
        return -1;
    }
    
    (*mat)->rows = rows - row_offset;
    (*mat)->cols = cols - col_offset;
    if (rows- row_offset == 1 || cols- col_offset == 1) {
	    (*mat)->is_1d = 1;
    } else {
	    (*mat)->is_1d = 0;
    }
    //(*mat)->is_1d = 0;
    (*mat)->parent = from;
    (*mat)->ref_cnt = 1;
    from->ref_cnt += 1; //increment that of parent


    //double *newCols = (double *) malloc(sizeof(double) * cols);
    
   // if (newCols == NULL) {
     //   return -1;
   // }
    
   // double **newData = (double **) malloc(sizeof(newCols) * rows);
    
	
    int w;
    double **outer = (double**)malloc(sizeof(double*)*rows);
    for (w = 0; w < rows - row_offset; w++) {
	    double* inner = (double*)malloc(sizeof(double)*cols);
	    outer[w] = inner;
    }

    if (outer == NULL) {
        return -1;
    }

    double **parentData = from->data;
    int i, j;

    for (i = 0; i < rows - row_offset; i++) {
        for (j = 0; j < cols - col_offset; j++) {
	    //consider edge case when offset +rows > from->rows
            outer[i][j] = parentData[row_offset + i][col_offset + j];
        }
    }

    (*mat)->data = outer;
    return 0;
    /* TODO: YOUR CODE HERE */
}

/*
 * This function will be called automatically by Python when a numc matrix loses all of its
 * reference pointers.
 * You need to make sure that you only free `mat->data` if no other existing matrices are also
 * referring this data array.
 * See the spec for more information.
 */
void deallocate_matrix(matrix *mat) {
    /* TODO: YOUR CODE HERE */
	if (mat == NULL) {
		return;
	}

    if (mat->ref_cnt == 0) {
        int i;
        double **data = mat->data;
        for (i = 0; i < mat->rows; i++) {
            free(data[i]);
        }
        free(data);
    }
}

/*
 * Return the double value of the matrix at the given row and column.
 * You may assume `row` and `col` are valid.
 */
double get(matrix *mat, int row, int col) {
    /* TODO: YOUR CODE HERE */
    return mat->data[row][col];
}

/*
 * Set the value at the given row and column to val. You may assume `row` and
 * `col` are valid
 */
void set(matrix *mat, int row, int col, double val) {
    /* TODO: YOUR CODE HERE */
    mat->data[row][col] = val;
    //column[row] = val;
}

/*
 * Set all entries in mat to val
 */
void fill_matrix(matrix *mat, double val) {
    /* TODO: YOUR CODE HERE */
    int col = mat->cols;
    int row = mat->rows;
    double **data = mat->data;
    int i, j;
    for (i=0; i < row; i++) {
    	for (j=0; j < col; j++) {
    		data[i][j] = val;
    	}
    }
}

/*
 * Store the result of adding mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int add_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    int rowMat1 = mat1->rows;
    int colMat1 = mat1->cols;
    int rowMat2 = mat2->rows;
    int colMat2 = mat2->cols;
    int colMatRes = result->cols;
    int rowMatRes = result->rows;
    int i, j;
    double **data1 = mat1->data;
    double **data2 = mat2->data;
    double **data3 = result->data;

    if (!(rowMat1 == rowMat2 && colMat1 == colMat2 
    	&& rowMat1 == rowMatRes && colMat1 == colMatRes)) {
    	PyErr_SetString(PyExc_ValueError, "Not valid dimensions");
        return -1;
    }

    for (i = 0; i < rowMatRes; i++) {

    	for (j = 0; j < colMatRes; j++) {
    		data3[i][j] = data1[i][j] + data2[i][j];
    	}
    }

    return 0;

}

/*
 * Store the result of subtracting mat2 from mat1 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int sub_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    int rowMat1 = mat1->rows;
    int colMat1 = mat1->cols;
    int rowMat2 = mat2->rows;
    int colMat2 = mat2->cols;
    int colMatRes = result->cols;
    int rowMatRes = result->rows;
    int i, j;
    double **data1 = mat1->data;
    double **data2 = mat2->data;
    double **data3 = result->data;

    if (!(rowMat1 == rowMat2 && colMat1 == colMat2 
    	&& rowMat1 == rowMatRes && colMat1 == colMatRes)) {
        PyErr_SetString(PyExc_ValueError, "Not valid dimensions");
    	return -1;
    }

    for (i = 0; i < rowMatRes; i++) {
    	for (j = 0; j < colMatRes; j++) {
    		data3[i][j] = data1[i][j] - data2[i][j];
    	}
    }
    
    return 0;
}

/*
 * Store the result of multiplying mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that matrix multiplication is not the same as multiplying individual elements.
 */
int mul_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    int matrix1Rows = mat1->rows;
    int matrix1Cols = mat1->cols;

    int matrix2Rows = mat2->rows;
    int matrix2Cols = mat2->cols;

    int matrixResRows = result->rows;
    int matrixResCols = result->cols;

    double **mat1Data = mat1->data;
    double **mat2Data = mat2->data;
    double **resData = result->data;
    int i, j, w;
    double sum;


    if (matrix1Cols != matrix2Rows || matrixResCols != matrix2Cols 
        || matrixResRows != matrix1Rows) {
        PyErr_SetString(PyExc_ValueError, "Not valid dimensions");
        return -1;
    }

    for (i = 0; i < matrix1Rows; i++) {

        for (j = 0; j < matrix2Rows; j++) {
            sum = 0;
            for (w = 0; w < matrix1Cols; w++) {
                sum += mat1Data[w][i] * mat2Data[j][w];
            }
            resData[j][i] = sum;
        }
    }
    return 0;

}

/*
 * Store the result of raising mat to the (pow)th power to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that pow is defined with matrix multiplication, not element-wise multiplication.
 */
//0th power is I
int pow_matrix(matrix *result, matrix *mat, int pow) {
    /* TODO: YOUR CODE HERE */
    if (pow < 0 || result->rows != mat->rows 
        || result->cols != mat->cols || mat->rows != mat->cols) {
        PyErr_SetString(PyExc_ValueError, "Not valid dimensions");
        return -1;
    }

    int first = pow;
    while (pow - 1 > 0) {
        if (pow == first) {
            mul_matrix(result, mat, mat);
        } else {
	    matrix *copy = NULL;
	    allocate_matrix_ref(&copy, result, 0, 0, result->rows, result->cols);
            mul_matrix(result, copy, mat);
	    deallocate_matrix(copy);
        }
        pow--; 
    }
    return 0;
}

/*
 * Store the result of element-wise negating mat's entries to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int neg_matrix(matrix *result, matrix *mat) {
    /* TODO: YOUR CODE HERE */
    int originalRow = mat->rows;
    int originalCol = mat->cols;
    int newRow = result->rows;
    int newCol = result->cols;
    double **dataOriginal = mat->data;
    double **dataNew = result->data;
    int i, j;

    if (!(originalRow == newRow && originalCol == newCol)) {
        PyErr_SetString(PyExc_ValueError, "Not valid dimensions");
    	return -1;
    }

    for (i = 0; i < newRow; i++) {
    	for (j = 0; j < newCol; j++) {
    		dataNew[i][j] = -dataOriginal[i][j];
    	}
    }

    return 0;
}

/*
 * Store the result of taking the absolute value element-wise to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int abs_matrix(matrix *result, matrix *mat) {
    /* TODO: YOUR CODE HERE */
    int originalRow = mat->rows;
    int originalCol = mat->cols;
    int newRow = result->rows;
    int newCol = result->cols;
    int i, j;
    double **dataOriginal = mat->data;
    double **dataNew = result->data;

    if (!(originalRow == newRow && originalCol == newCol)) {
        PyErr_SetString(PyExc_ValueError, "Not valid dimensions");
    	return -1;
    }

    for (i = 0; i < newRow; i++) {
    	for (j = 0; j < newCol; j++) {
    		double beforeAbs = dataOriginal[i][j];
    		if (beforeAbs >= 0) {
    			dataNew[i][j] = beforeAbs;
    		} else {
    			dataNew[i][j] = -beforeAbs;
    		}
    	}
    }

    return 0;

}

