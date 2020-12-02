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

    if (rows <= 0 || cols <= 0) {
        PyErr_SetString(PyExc_ValueError, "Given dimensions are not valid!");
        return -1;
    }
    *mat = (matrix *) malloc(sizeof(matrix));

    if (mat == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Could not allocate, not enough space!");
        return -1;
    }
    (*mat)->rows = rows;
    (*mat)->cols = cols;
    (*mat)->parent = NULL;
    (*mat)->ref_cnt = 1;
    if (rows == 1 || cols == 1) {
        (*mat)->is_1d = 1;
    } else {
        (*mat)->is_1d = 0;
    }
    int i;
    double **outer = (double **)calloc(rows, sizeof(double*));
    double *inner = (double *)calloc(cols * rows, sizeof(double));
    if (outer == NULL || inner == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Could not allocate, not enough space!");
        return -1;
    }
    for (i = 0; i < rows; i++) {  
        outer[i] = inner + (i * cols);
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

    if (rows <= 0 || cols <= 0 
    || rows > from->rows || cols > from->cols 
    || row_offset < 0 || col_offset < 0) {
        PyErr_SetString(PyExc_ValueError, "Given dimensions are not valid!");
        return -1;
    }

    *mat = (matrix *) malloc(sizeof(matrix));

    if (mat == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Could not allocate, not enough space!");
        return -1;
    }

    (*mat)->rows = rows;
    (*mat)->cols = cols;
    if (rows == 1 || cols == 1) {
        (*mat)->is_1d = 1;
    } else {
        (*mat)->is_1d = 0;
    }

    (*mat)->parent = from;
    (*mat)->ref_cnt = 1;
    from->ref_cnt += 1;
    double **outer = (double **)calloc(rows, sizeof(double*));
    double *inner = (double *)calloc(cols * rows, sizeof(double));
    if (outer == NULL || inner == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Could not allocate, not enough space!");
        return -1;
    }  
    int w,i,j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            inner[(i * cols) + j] = from->data[i + row_offset][j + col_offset];
        }
    }
    for (w = 0; w < rows; w++) {  
        outer[w] = inner + (w * cols);
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
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
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

    if (!(rowMat1 == rowMat2 && colMat1 == colMat2
            && rowMat1 == rowMatRes && colMat1 == colMatRes)) {
        PyErr_SetString(PyExc_ValueError, "Not valid dimensions");
        return -1;
    }

    
    #pragma omp parallel for if(rowMat1 * colMat1 > 100000)
    for (int i = 0; i < result->rows * result->cols; i ++) {
        result->data[0][i] = mat1->data[0][i] + mat2->data[0][i];
    }
    // for (int i = 0; i < result->rows; i++) {
    //     for (int j = 0; j < result->cols; j++) {
    //         result->data[i][j] = mat1->data[i][j] + mat2->data[i][j];
    //     }
    // }

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

    if (!(rowMat1 == rowMat2 && colMat1 == colMat2
            && rowMat1 == rowMatRes && colMat1 == colMatRes)) {
        PyErr_SetString(PyExc_ValueError, "Not valid dimensions");
        return -1;
    }

    #pragma omp parallel for if(rowMat1 * colMat1 > 100000)
    for (int i = 0; i < result->rows * result->cols; i ++) {
        result->data[0][i] = mat1->data[0][i] - mat2->data[0][i];
    }
    // for (int i = 0; i < result->rows; i++) {
    //     for (int j = 0; j < result->cols; j++) {
    //         result->data[i][j] = mat1->data[i][j] - mat2->data[i][j];
    //     }
    // }
    return 0;
}

/*
 * Store the result of multiplying mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that matrix multiplication is not the same as multiplying individual elements.
 */
int mul_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */

    //1. cache block transpose
    //2. simd
    //3. unroll 
    //4. openMP
    int matrix1Rows = mat1->rows;
    int matrix1Cols = mat1->cols;

    int matrix2Rows = mat2->rows;
    int matrix2Cols = mat2->cols;

    int matrixResRows = result->rows;
    int matrixResCols = result->cols;

    //double **mat1Data = mat1->data;
    //double **mat2Data = mat2->data;
    //double **resData = result->data;
    //double sum;


    if (matrix1Cols != matrix2Rows || matrixResCols != matrix2Cols
            || matrixResRows != matrix1Rows) {
        return -1;
    }

    matrix *copy = NULL;
    allocate_matrix(&copy, mat2->cols, mat2->rows);

    #pragma omp parallel for if(matrix2Cols * matrix2Rows > 100000)
    for (int i = 0; i < (matrix2Cols / 4) * 4 ; i += 4) {
        for (int j = 0; j < (matrix2Rows / 4) * 4; j += 4) {
            set(copy, i, j, mat2->data[j][i]);
            set(copy, i, j+1, mat2->data[j+1][i]);
            set(copy, i, j+2, mat2->data[j+2][i]);
            set(copy, i, j+3, mat2->data[j+3][i]);
        }

        for (int j = matrix2Rows - (matrix2Rows % 4); j < matrix2Rows; j++) {
            set(copy, i, j, mat2->data[j][i]);
        }

        for (int j = 0; j < (matrix2Rows / 4) * 4; j += 4) {
            set(copy, i+1, j, mat2->data[j][i+1]);
            set(copy, i+1, j+1, mat2->data[j+1][i+1]);
            set(copy, i+1, j+2, mat2->data[j+2][i+1]);
            set(copy, i+1, j+3, mat2->data[j+3][i+1]);
        }

        for (int j = matrix2Rows - (matrix2Rows % 4); j < matrix2Rows; j++) {
            set(copy, i+1, j, mat2->data[j][i+1]);
        }

        for (int j = 0; j < (matrix2Rows / 4) * 4; j += 4) {
            set(copy, i+2, j, mat2->data[j][i+2]);
            set(copy, i+2, j+1, mat2->data[j+1][i+2]);
            set(copy, i+2, j+2, mat2->data[j+2][i+2]);
            set(copy, i+2, j+3, mat2->data[j+3][i+2]);
        }

        for (int j = matrix2Rows - (matrix2Rows % 4); j < matrix2Rows; j++) {
            set(copy, i+2, j, mat2->data[j][i+2]);
        }

        for (int j = 0; j < (matrix2Rows / 4) * 4; j += 4) {
            set(copy, i+3, j, mat2->data[j][i+3]);
            set(copy, i+3, j+1, mat2->data[j+1][i+3]);
            set(copy, i+3, j+2, mat2->data[j+2][i+3]);
            set(copy, i+3, j+3, mat2->data[j+3][i+3]);
        }

        for (int j = matrix2Rows - (matrix2Rows % 4); j < matrix2Rows; j++) {
            set(copy, i+3, j, mat2->data[j][i+3]);
            }
        }
    

    #pragma omp parallel for if((matrix2Cols - ((matrix2Cols / 4) * 4)) * matrix2Rows > 100000)    
    for (int i = matrix2Cols - (matrix2Cols % 8); i < matrix2Cols; i++) {
        for (int j = 0; j < (matrix2Rows / 4) * 4; j += 4) {
            set(copy, i, j, mat2->data[j][i]);
            set(copy, i, j+1, mat2->data[j+1][i]);
            set(copy, i, j+2, mat2->data[j+2][i]);
            set(copy, i, j+3, mat2->data[j+3][i]);
        }

        for (int j = matrix2Rows - (matrix2Rows % 4); j < matrix2Rows; j++) {
            set(copy, i, j, mat2->data[j][i]);
        }
    }



    //The following until line until 443 is an unrolled SIMD loop that does the matrix multiplication
    #pragma omp parallel for if(matrix1Rows * matrix2Cols > 100000)
    for (int i = 0; i < mat1->rows; i++) {
        for (int j = 0; j < copy->rows; j++) {
            set(result, i, j, 0);
            double p[4] = {0,0,0,0};
            __m256d z = _mm256_set1_pd(0);
            for (int k = 0; k < (copy->cols / 4) * 4; k += 4) {
                __m256d x = _mm256_loadu_pd(&mat1->data[i][k]);
                __m256d y = _mm256_loadu_pd(&copy->data[j][k]);
                z = _mm256_fmadd_pd(x, y, z);
            }
            _mm256_storeu_pd(p,z);
            result->data[i][j] += p[0] + p[1] + p[2] + p[3];
            //printf("%f\n", result->data[i][j]);
        }
    }

    //This is the tail case of the above, your actual matrix multiplication
    //for a small matrix will more than likely take place here 
    #pragma omp parallel for if((copy->cols - ((copy->cols / 4) * 4))*copy->rows > 100000)
    for (int i = 0; i < mat1->rows; i++) {
        for (int j = 0; j < copy->rows; j++) {
            for (int k = copy->cols - (copy->cols % 4); k < copy->cols; k++) {
                result->data[i][j] += mat1->data[i][k] * copy->data[j][k];
            }
            //printf("%f\n", result->data[i][j]);
        }
    }

    deallocate_matrix(copy);
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
    int colMatRes = result->cols;
    int rowMatRes = result->rows;
    int i, j;

    for (i = 0; i < rowMatRes; i++) {
    	for (j = 0; j < colMatRes; j++) {
            if (i == j) {
                set(result,i,j,1);
            } else {
                set(result,i,j,0);
            }
    	}
    }

    if (pow == 0) {
        return 0;
    }

    if (pow == 1) {
        result = mat;
        return 0;
    }

    matrix *product = NULL;
    allocate_matrix(&product, rowMatRes, colMatRes);

    matrix *mat_copy = NULL;
    allocate_matrix(&mat_copy, rowMatRes, colMatRes);
    for (i = 0; i < result->rows; i++) {
		    for (j = 0; j < result->cols; j++) {
			    mat_copy->data[i][j] = mat->data[i][j];
		}
	}


    while (pow > 1) {
        if (pow % 2 != 0) {
            mul_matrix(product, result, mat_copy);
            mat = result;
            result = product;
            product = mat;
            pow = pow/2;
        } else {
            mul_matrix(result, mat_copy ,result);
            mul_matrix(product, mat_copy, mat_copy);
            mat = mat_copy;
            mat_copy = product;
            product = mat_copy;
            pow = (pow-1)/2;
        }
    }

    mul_matrix(result, mat_copy, result);

    deallocate_matrix(product);
    deallocate_matrix(mat_copy);

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

    if (!(originalRow == newRow && originalCol == newCol)) {
        PyErr_SetString(PyExc_ValueError, "Not valid dimensions");
        return -1;
    }

     
    #pragma omp parallel for if(originalRow * originalCol > 100000)
    for (int i = 0; i < result->rows * result->cols; i++) {
        result->data[0][i] = -mat->data[0][i];
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

    if (!(originalRow == newRow && originalCol == newCol)) {
        PyErr_SetString(PyExc_ValueError, "Not valid dimensions");
        return -1;
    }
    
    #pragma omp parallel for if(originalRow * originalCol > 100000)
    for (int i = 0; i < result->rows * result->cols; i++) {
        if (mat->data[0][i] >= 0) {
                result->data[0][i] = mat->data[0][i];
            }else {
                result->data[0][i] = -mat->data[0][i];
            } 
    }

        

    return 0;

}

