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
        free(mat);
	return -1;
    }

    #pragma omp parallel for if(rows * cols > 100000)
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
        free(mat);
	return -1;
    } 


    #pragma omp parallel for if(rows * cols > 100000)
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            inner[(i * cols) + j] = from->data[i + row_offset][j + col_offset];
        }
    }
    
    #pragma omp parallel for if(rows > 100000)
    for (int w = 0; w < rows; w++) {  
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
    } else {
        if (mat->parent == NULL) {
            if (mat->ref_cnt == 1) {
                free(*(mat->data));
                free(mat->data);
                free(mat);
            } else {
                (mat->ref_cnt)--;
            }
        } else {
            if (mat->parent->ref_cnt == 1) {
                free(*(mat->parent->data));
                free(mat->parent->data);
                free(mat->parent);
            }
            (mat->parent->ref_cnt)--;
            free(*(mat->data));
            free(mat->data);
            free(mat);
        }
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
    double * mat1Inner = *(mat1->data);
    double * mat2Inner = *(mat2->data);
    double * resInner = *(result->data);
    int size = rowMat1 * colMat1;

    if (!(rowMat1 == rowMat2 && colMat1 == colMat2
            && rowMat1 == rowMatRes && colMat1 == colMatRes)) {
        PyErr_SetString(PyExc_ValueError, "Not valid dimensions");
        return -1;
    }

    __m256d mat1block = _mm256_set1_pd(0);
    __m256d mat2block = _mm256_set1_pd(0); 
    __m256d resultblock =_mm256_set1_pd(0);

    #pragma omp parallel for if(size > 10000)
    for (int i = 0 ; i < (size / 4) * 4; i+=4){
	    mat1block =_mm256_loadu_pd(&mat1Inner[i]);
	    mat2block = _mm256_loadu_pd(&mat2Inner[i]);
	    resultblock =_mm256_add_pd(mat1block,mat2block);
	    _mm256_storeu_pd(&resInner[i], resultblock);
    }

    for (int i = size - (size % 4); i < size; i++) {
        resInner[i] = mat1Inner[i] + mat2Inner[i];
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
    double * mat1Inner = *(mat1->data);
    double * mat2Inner = *(mat2->data);
    double * resInner = *(result->data);
    int size = rowMat1 * colMat1;

    if (!(rowMat1 == rowMat2 && colMat1 == colMat2
            && rowMat1 == rowMatRes && colMat1 == colMatRes)) {
        PyErr_SetString(PyExc_ValueError, "Not valid dimensions");
        return -1;
    }

    __m256d mat1block = _mm256_set1_pd(0);
    __m256d mat2block = _mm256_set1_pd(0); 
    __m256d resultblock =_mm256_set1_pd(0);

    #pragma omp parallel for if(size > 10000)
    for (int i = 0 ; i < (size / 4) * 4; i+=4){
	    mat1block =_mm256_loadu_pd(&mat1Inner[i]);
	    mat2block = _mm256_loadu_pd(&mat2Inner[i]);
	    resultblock =_mm256_sub_pd(mat1block,mat2block);
	    _mm256_storeu_pd(&resInner[i], resultblock);
    }

    for (int i = size - (size % 4); i < size; i++) {
        resInner[i] = mat1Inner[i] - mat2Inner[i];
    }

    return 0;
}

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

    if (matrix1Cols != matrix2Rows || matrixResCols != matrix2Cols
            || matrixResRows != matrix1Rows) {
        return -1;
    }

    #pragma omp parallel for if(matrix1Rows * matrix2Cols > 25000)
    for (int i = 0; i < matrixResCols * matrixResRows; i++) {
        result->data[0][i] = 0;
    }


    #pragma omp parallel for if(matrix1Rows * matrix2Cols > 25000)
    for (int i = 0; i < matrix1Rows; i++) {
        for (int j = 0; j < matrix1Cols; j++) {
            for (int w = 0; w < matrix2Cols; w++) {
                resData[i][w] += mat1Data[i][j] * mat2Data[j][w];
            }
        }
    }

    return 0;
}

int pow_mul(matrix *result, matrix *mat1, matrix *mat2) {
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

    if (matrix1Cols != matrix2Rows || matrixResCols != matrix2Cols
            || matrixResRows != matrix1Rows) {
        return -1;
    }

    #pragma omp parallel for if(matrix1Rows * matrix2Cols > 6000)
    for (int i = 0; i < matrixResCols * matrixResRows; i++) {
        result->data[0][i] = 0;
    }


    #pragma omp parallel for if(matrix1Rows * matrix2Cols > 6000)
    for (int i = 0; i < matrix1Rows; i++) {
        for (int j = 0; j < matrix1Cols; j++) {
            for (int w = 0; w < matrix2Cols; w++) {
                resData[i][w] += mat1Data[i][j] * mat2Data[j][w];
            }
        }
    }

    return 0;
}

int pow_mul2(matrix *result, matrix *mat1, matrix *mat2) {
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

    if (matrix1Cols != matrix2Rows || matrixResCols != matrix2Cols
            || matrixResRows != matrix1Rows) {
        return -1;
    }

    #pragma omp parallel for if(matrix1Rows * matrix2Cols > 6000)
    for (int i = 0; i < matrixResCols * matrixResRows; i++) {
        result->data[0][i] = 0;
    }


    #pragma omp parallel for if(matrix1Rows * matrix2Cols > 6000)
    for (int i = 0; i < matrix1Rows; i++) {
        for (int j = 0; j < (matrix1Cols / 4) * 4; j+=4) {
            for (int w = 0; w < matrix2Cols; w++) {
                resData[i][w] += mat1Data[i][j] * mat2Data[j][w];
                resData[i][w] += mat1Data[i][j+1] * mat2Data[j+1][w];
                resData[i][w] += mat1Data[i][j+2] * mat2Data[j+2][w];
                resData[i][w] += mat1Data[i][j+3] * mat2Data[j+3][w];
            }
        }
        for (int j = matrix1Cols - (matrix1Cols % 4); j < matrix1Cols; j++) {
            for (int w = 0; w < matrix2Cols; w++) {
                resData[i][w] += mat1Data[i][j] * mat2Data[j][w];
            }
        }

    }
    return 0;
}



/*
 * Store the result of multiplying mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that matrix multiplication is not the same as multiplying individual elements.
 */
int mul_matrix2(matrix *result, matrix *mat1, matrix *mat2) {
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

    double **mat1Data = mat1->data;
    double **mat2Data = mat2->data;
    double **resData = result->data;
    //double sum;


    if (matrix1Cols != matrix2Rows || matrixResCols != matrix2Cols
            || matrixResRows != matrix1Rows) {
        return -1;
    }

    matrix *copy = NULL;
    allocate_matrix(&copy, mat2->cols, mat2->rows);
    int matrixCopyRows = copy->rows;
    int matrixCopyCols = copy->cols;

    double **copyData = copy->data;

    double *p1;

    #pragma omp parallel for if(matrix2Rows * matrix2Cols > 30000)
    for (int i = 0; i < (matrix2Rows / 8) * 8; i++) {
        for (int j = 0; j < (matrix2Cols / 8) * 8; j+=8) {
            __m256d mat2DataBlock = _mm256_loadu_pd(&mat2Data[i][j]);
            p1 = (double *) &mat2DataBlock;
            copyData[j][i] = p1[0];
            copyData[j+1][i] = p1[1];
            copyData[j+2][i] = p1[2];
            copyData[j+3][i] = p1[3];

            mat2DataBlock = _mm256_loadu_pd(&mat2Data[i][j + 4]);
            p1 = (double *) &mat2DataBlock;
            copyData[j+4][i] = p1[0];
            copyData[j+1+4][i] = p1[1];
            copyData[j+2+4][i] = p1[2];
            copyData[j+3+4][i] = p1[3];
        }
    }

    #pragma omp parallel for if(matrix2Cols > 30000)
    for (int i = matrix2Rows - (matrix2Rows % 8); i < matrix2Rows; i++) {
        for (int j = 0; j < matrix2Cols; j ++) {
            copyData[j][i] = mat2Data[i][j];
        }
    }

    #pragma omp parallel for if(matrix2Rows > 30000)
    for (int i = 0; i < (matrix2Rows / 8) * 8; i++) {
        for (int j = matrix2Cols - (matrix2Cols % 8); j < matrix2Cols; j++) {
            copyData[j][i] = mat2Data[i][j];
        }
    }

    //The following until line until 443 is an unrolled SIMD loop that does the matrix multiplication
    double *p;

    #pragma omp parallel for if(matrix1Rows * matrix2Cols > 30000)
    for (int i = 0; i < mat1->rows; i++) {
        double * ptr = mat1Data[i];
        double * resRowBlock = resData[i];
        for (int j = 0; j < matrixCopyRows; j++) {
            double * ptrc = copyData[j];
            resRowBlock[j] = 0;
            __m256d z = _mm256_set1_pd(0);
            for (int k = 0; k < (matrixCopyCols / 8) * 8; k += 8) {
                __m256d x = _mm256_loadu_pd(&mat1Data[i][k]);
                __m256d y = _mm256_loadu_pd(&copyData[j][k]);
                z = _mm256_fmadd_pd(x, y, z);

                x = _mm256_loadu_pd(&mat1Data[i][k+4]);
                y = _mm256_loadu_pd(&copyData[j][k+4]);
                z = _mm256_fmadd_pd(x, y, z);
            }

            p = (double *) &z;
            resRowBlock[j] += p[0] + p[1] + p[2] + p[3];

            for (int k = matrixCopyCols - (matrixCopyCols % 8); k < matrixCopyCols; k++) {
                resRowBlock[j] += ptr[k] * ptrc[k];
            }
            //printf("%f\n", resData[i][j]);
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
int pow_matrix2(matrix *result, matrix *mat, int pow) {
    /* TODO: YOUR CODE HERE */
        if (pow < 0 || result->rows != mat->rows 
        || result->cols != mat->cols || mat->rows != mat->cols) {
        PyErr_SetString(PyExc_ValueError, "Not valid dimensions");
        return -1;
    }
    int res_swap = 0;
    int colMatRes = result->cols;
    int rowMatRes = result->rows;

    double **resData = result->data;

    #pragma omp parallel for if(rowMatRes * colMatRes > 20000)
    for (int i = 0; i < rowMatRes; i++) {
        double * currRow = resData[i];
    	for (int j = 0; j < colMatRes; j++) {
            if (i == j) {
                currRow[j] = 1;
            } else {
                currRow[j] = 0;
            }
    	}
    }

    if (pow == 0) {
        return 0;
    }

    double * resDataRow = resData[0];
    double * matRow = mat->data[0];

    if (pow == 1) {
        #pragma omp parallel for if(rowMatRes * colMatRes > 20000)
        for (int i = 0; i < rowMatRes * colMatRes; i++) {
            resDataRow[i] = matRow[i];
        }
        return 0;   
    }

    if (pow == 2) {
        mul_matrix(result, mat, mat);
        return 0;
    }

    matrix *mat_copy_temp = NULL;
    allocate_matrix(&mat_copy_temp, rowMatRes, colMatRes);

    matrix *result_temp = NULL;
    allocate_matrix(&result_temp, rowMatRes, colMatRes);

    matrix *mat_copy = NULL;
    allocate_matrix(&mat_copy, rowMatRes, colMatRes);
    double * matCopyRow = mat_copy->data[0];

    #pragma omp parallel for if(rowMatRes * colMatRes > 20000)
    for (int i = 0; i < rowMatRes * colMatRes; i++) {
        matCopyRow[i] = matRow[i];
    }

    matrix * temp_pointer;

    while (pow > 1) {
        if (pow % 2 == 0) {
            temp_pointer = mat_copy;
            mat_copy = mat_copy_temp;
            mat_copy_temp = temp_pointer;
            mul_matrix(mat_copy, mat_copy_temp, mat_copy_temp);
            pow = pow / 2;

        } else {
            temp_pointer = result;
            result = result_temp;
            result_temp = temp_pointer;
            mul_matrix(result, result_temp, mat_copy);
            res_swap++;

            temp_pointer = mat_copy;
            mat_copy = mat_copy_temp;
            mat_copy_temp = temp_pointer;
            mul_matrix(mat_copy, mat_copy_temp, mat_copy_temp);
            pow = (pow - 1) / 2;
        }
    }

    temp_pointer = result;
    result = result_temp;
    result_temp = temp_pointer;
    mul_matrix(result, result_temp, mat_copy);

    double * currResTempRow = result_temp->data[0]; 
    double * currResRow = result->data[0];

    if (res_swap % 2 == 0) {
        #pragma omp parallel for if(rowMatRes * colMatRes > 20000)
        for (int i = 0; i < rowMatRes * colMatRes; i++) {
            currResTempRow[i] = currResRow[i];
        }

        temp_pointer = result;
        result = result_temp;
        result_temp = temp_pointer;
    }

    deallocate_matrix(mat_copy_temp);
    deallocate_matrix(result_temp);
    deallocate_matrix(mat_copy);
    return 0;
}

int pow_matrix(matrix *result, matrix *mat, int pow) {
    /* TODO: YOUR CODE HERE */
        if (pow < 0 || result->rows != mat->rows 
        || result->cols != mat->cols || mat->rows != mat->cols) {
        PyErr_SetString(PyExc_ValueError, "Not valid dimensions");
        return -1;
    }
    int res_swap = 0;
    int colMatRes = result->cols;
    int rowMatRes = result->rows;

    double **resData = result->data;

    #pragma omp parallel for if(rowMatRes * colMatRes > 20000)
    for (int i = 0; i < rowMatRes; i++) {
        double * currRow = resData[i];
    	for (int j = 0; j < colMatRes; j++) {
            if (i == j) {
                currRow[j] = 1;
            } else {
                currRow[j] = 0;
            }
    	}
    }

    if (pow == 0) {
        return 0;
    }

    double * resDataRow = resData[0];
    double * matRow = mat->data[0];

    if (pow == 1) {
        #pragma omp parallel for if(rowMatRes * colMatRes > 20000)
        for (int i = 0; i < rowMatRes * colMatRes; i++) {
            resDataRow[i] = matRow[i];
        }
        return 0;   
    }

    if (pow == 2) {
        pow_mul2(result, mat, mat);
        return 0;
    }

    matrix *mat_copy_temp = NULL;
    allocate_matrix(&mat_copy_temp, rowMatRes, colMatRes);

    matrix *result_temp = NULL;
    allocate_matrix(&result_temp, rowMatRes, colMatRes);

    matrix *mat_copy = NULL;
    allocate_matrix(&mat_copy, rowMatRes, colMatRes);
    double * matCopyRow = mat_copy->data[0];

    #pragma omp parallel for if(rowMatRes * colMatRes > 20000)
    for (int i = 0; i < rowMatRes * colMatRes; i++) {
        matCopyRow[i] = matRow[i];
    }

    matrix * temp_pointer;
    int first_odd = 1;

    while (pow > 0) {
        if (pow % 2 != 0) {

            if (first_odd == 1) {
                for (int i = 0; i < rowMatRes * colMatRes; i++) {
                    resDataRow[i] = mat_copy->data[0][i];
                }
                first_odd = 0;
            } else {
                temp_pointer = result;
                result = result_temp;
                result_temp = temp_pointer;
                pow_mul2(result, result_temp, mat_copy);
                res_swap++;
            }
            pow--;

        } else {
            temp_pointer = mat_copy;
            mat_copy = mat_copy_temp;
            mat_copy_temp = temp_pointer;
            pow_mul2(mat_copy, mat_copy_temp, mat_copy_temp);
            pow = pow / 2;
        }
    }

    double * currResTempRow = result_temp->data[0]; 
    double * currResRow = result->data[0];

    if (res_swap % 2 != 0) {
        #pragma omp parallel for if(rowMatRes * colMatRes > 20000)
        for (int i = 0; i < rowMatRes * colMatRes; i++) {
            currResTempRow[i] = currResRow[i];
        }

        temp_pointer = result;
        result = result_temp;
        result_temp = temp_pointer;
    }

    deallocate_matrix(mat_copy_temp);
    deallocate_matrix(result_temp);
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
    int size = newRow * newCol;

    if (!(originalRow == newRow && originalCol == newCol)) {
        PyErr_SetString(PyExc_ValueError, "Not valid dimensions");
        return -1;
    }

    double * matInner = *(mat->data);
    double * resInner = *(result->data);


    __m256d negblock = _mm256_set1_pd(-1);     
    __m256d matblock = _mm256_set1_pd(0);
    __m256d resultblock =_mm256_set1_pd(0);

    #pragma omp parallel for if(size > 10000)
    for (int i = 0 ; i < (size / 4) * 4; i+=4){
	    matblock =_mm256_loadu_pd(&matInner[i]);
	    resultblock =_mm256_mul_pd(matblock,negblock);
	    _mm256_storeu_pd(&resInner[i], resultblock);
    }

    for (int i = size - (size % 4); i < size; i++) {
        resInner[i] = -matInner[i];
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
    int size = newRow * newCol;

    if (!(originalRow == newRow && originalCol == newCol)) {
        PyErr_SetString(PyExc_ValueError, "Not valid dimensions");
        return -1;
    }

    double * matInner = *(mat->data);
    double * resInner = *(result->data);


    __m256d zeroblock = _mm256_set1_pd(0);     
    __m256d matblock = _mm256_set1_pd(0);
    __m256d resultblock =_mm256_set1_pd(0);

    #pragma omp parallel for if(size > 10000)
    for (int i = 0 ; i < (size / 4) * 4; i+=4){
	    matblock =_mm256_loadu_pd(&matInner[i]);
	    resultblock =_mm256_sub_pd(zeroblock,matblock);
        resultblock = _mm256_max_pd(resultblock,matblock);
	    _mm256_storeu_pd(&resInner[i], resultblock);
    }
    
    for (int i = size - (size % 4); i < size; i++) {
        if (matInner[i] >= 0) {
                resInner[i] = matInner[i];
            }else {
                resInner[i] = -matInner[i];
            } 
    }
    return 0;
}
