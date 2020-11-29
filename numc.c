#include "numc.h"
#include <structmember.h>

PyTypeObject Matrix61cType;

/* Helper functions for initalization of matrices and vectors */

/*
 * Return a tuple given rows and cols
 */
PyObject *get_shape(int rows, int cols) {
    if (rows == 1 || cols == 1) {
        return PyTuple_Pack(1, PyLong_FromLong(rows * cols));
    } else {
        return PyTuple_Pack(2, PyLong_FromLong(rows), PyLong_FromLong(cols));
    }
}
/*
 * Matrix(rows, cols, low, high). Fill a matrix random double values
 */
int init_rand(PyObject *self, int rows, int cols, unsigned int seed, double low,
              double high) {
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed) return alloc_failed;
    rand_matrix(new_mat, seed, low, high);
    ((Matrix61c *)self)->mat = new_mat;
    ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    return 0;
}

/*
 * Matrix(rows, cols, val). Fill a matrix of dimension rows * cols with val
 */
int init_fill(PyObject *self, int rows, int cols, double val) {
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed)
        return alloc_failed;
    else {
        fill_matrix(new_mat, val);
        ((Matrix61c *)self)->mat = new_mat;
        ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    }
    return 0;
}

/*
 * Matrix(rows, cols, 1d_list). Fill a matrix with dimension rows * cols with 1d_list values
 */
int init_1d(PyObject *self, int rows, int cols, PyObject *lst) {
    if (rows * cols != PyList_Size(lst)) {
        PyErr_SetString(PyExc_ValueError, "Incorrect number of elements in list");
        return -1;
    }
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed) return alloc_failed;
    int count = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            set(new_mat, i, j, PyFloat_AsDouble(PyList_GetItem(lst, count)));
            count++;
        }
    }
    ((Matrix61c *)self)->mat = new_mat;
    ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    return 0;
}

/*
 * Matrix(2d_list). Fill a matrix with dimension len(2d_list) * len(2d_list[0])
 */
int init_2d(PyObject *self, PyObject *lst) {
    int rows = PyList_Size(lst);
    if (rows == 0) {
        PyErr_SetString(PyExc_ValueError,
                        "Cannot initialize numc.Matrix with an empty list");
        return -1;
    }
    int cols;
    if (!PyList_Check(PyList_GetItem(lst, 0))) {
        PyErr_SetString(PyExc_ValueError, "List values not valid");
        return -1;
    } else {
        cols = PyList_Size(PyList_GetItem(lst, 0));
    }
    for (int i = 0; i < rows; i++) {
        if (!PyList_Check(PyList_GetItem(lst, i)) ||
                PyList_Size(PyList_GetItem(lst, i)) != cols) {
            PyErr_SetString(PyExc_ValueError, "List values not valid");
            return -1;
        }
    }
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed) return alloc_failed;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            set(new_mat, i, j,
                PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(lst, i), j)));
        }
    }
    ((Matrix61c *)self)->mat = new_mat;
    ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    return 0;
}

/*
 * This deallocation function is called when reference count is 0
 */
void Matrix61c_dealloc(Matrix61c *self) {
    deallocate_matrix(self->mat);
    Py_TYPE(self)->tp_free(self);
}

/* For immutable types all initializations should take place in tp_new */
PyObject *Matrix61c_new(PyTypeObject *type, PyObject *args,
                        PyObject *kwds) {
    /* size of allocated memory is tp_basicsize + nitems*tp_itemsize*/
    Matrix61c *self = (Matrix61c *)type->tp_alloc(type, 0);
    return (PyObject *)self;
}

/*
 * This matrix61c type is mutable, so needs init function. Return 0 on success otherwise -1
 */
int Matrix61c_init(PyObject *self, PyObject *args, PyObject *kwds) {
    /* Generate random matrices */
    if (kwds != NULL) {
        PyObject *rand = PyDict_GetItemString(kwds, "rand");
        if (!rand) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
        if (!PyBool_Check(rand)) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
        if (rand != Py_True) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }

        PyObject *low = PyDict_GetItemString(kwds, "low");
        PyObject *high = PyDict_GetItemString(kwds, "high");
        PyObject *seed = PyDict_GetItemString(kwds, "seed");
        double double_low = 0;
        double double_high = 1;
        unsigned int unsigned_seed = 0;

        if (low) {
            if (PyFloat_Check(low)) {
                double_low = PyFloat_AsDouble(low);
            } else if (PyLong_Check(low)) {
                double_low = PyLong_AsLong(low);
            }
        }

        if (high) {
            if (PyFloat_Check(high)) {
                double_high = PyFloat_AsDouble(high);
            } else if (PyLong_Check(high)) {
                double_high = PyLong_AsLong(high);
            }
        }

        if (double_low >= double_high) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }

        // Set seed if argument exists
        if (seed) {
            if (PyLong_Check(seed)) {
                unsigned_seed = PyLong_AsUnsignedLong(seed);
            }
        }

        PyObject *rows = NULL;
        PyObject *cols = NULL;
        if (PyArg_UnpackTuple(args, "args", 2, 2, &rows, &cols)) {
            if (rows && cols && PyLong_Check(rows) && PyLong_Check(cols)) {
                return init_rand(self, PyLong_AsLong(rows), PyLong_AsLong(cols), unsigned_seed, double_low,
                                 double_high);
            }
        } else {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
    }
    PyObject *arg1 = NULL;
    PyObject *arg2 = NULL;
    PyObject *arg3 = NULL;
    if (PyArg_UnpackTuple(args, "args", 1, 3, &arg1, &arg2, &arg3)) {
        /* arguments are (rows, cols, val) */
        if (arg1 && arg2 && arg3 && PyLong_Check(arg1) && PyLong_Check(arg2) && (PyLong_Check(arg3)
                || PyFloat_Check(arg3))) {
            if (PyLong_Check(arg3)) {
                return init_fill(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), PyLong_AsLong(arg3));
            } else
                return init_fill(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), PyFloat_AsDouble(arg3));
        } else if (arg1 && arg2 && arg3 && PyLong_Check(arg1) && PyLong_Check(arg2) && PyList_Check(arg3)) {
            /* Matrix(rows, cols, 1D list) */
            return init_1d(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), arg3);
        } else if (arg1 && PyList_Check(arg1) && arg2 == NULL && arg3 == NULL) {
            /* Matrix(rows, cols, 1D list) */
            return init_2d(self, arg1);
        } else if (arg1 && arg2 && PyLong_Check(arg1) && PyLong_Check(arg2) && arg3 == NULL) {
            /* Matrix(rows, cols, 1D list) */
            return init_fill(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), 0);
        } else {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return -1;
    }
}

/*
 * List of lists representations for matrices
 */
PyObject *Matrix61c_to_list(Matrix61c *self) {
    int rows = self->mat->rows;
    int cols = self->mat->cols;
    PyObject *py_lst = NULL;
    if (self->mat->is_1d) {  // If 1D matrix, print as a single list
        py_lst = PyList_New(rows * cols);
        int count = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                PyList_SetItem(py_lst, count, PyFloat_FromDouble(get(self->mat, i, j)));
                count++;
            }
        }
    } else {  // if 2D, print as nested list
        py_lst = PyList_New(rows);
        for (int i = 0; i < rows; i++) {
            PyList_SetItem(py_lst, i, PyList_New(cols));
            PyObject *curr_row = PyList_GetItem(py_lst, i);
            for (int j = 0; j < cols; j++) {
                PyList_SetItem(curr_row, j, PyFloat_FromDouble(get(self->mat, i, j)));
            }
        }
    }
    return py_lst;
}

PyObject *Matrix61c_class_to_list(Matrix61c *self, PyObject *args) {
    PyObject *mat = NULL;
    if (PyArg_UnpackTuple(args, "args", 1, 1, &mat)) {
        if (!PyObject_TypeCheck(mat, &Matrix61cType)) {
            PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
            return NULL;
        }
        Matrix61c* mat61c = (Matrix61c*)mat;
        return Matrix61c_to_list(mat61c);
    } else {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return NULL;
    }
}

/*
 * Add class methods
 */
PyMethodDef Matrix61c_class_methods[] = {
    {"to_list", (PyCFunction)Matrix61c_class_to_list, METH_VARARGS, "Returns a list representation of numc.Matrix"},
    {NULL, NULL, 0, NULL}
};

/*
 * Matrix61c string representation. For printing purposes.
 */
PyObject *Matrix61c_repr(PyObject *self) {
    PyObject *py_lst = Matrix61c_to_list((Matrix61c *)self);
    return PyObject_Repr(py_lst);
}

/* NUMBER METHODS */

/*
 * Add the second numc.Matrix (Matrix61c) object to the first one. The first operand is
 * self, and the second operand can be obtained by casting `args`.
 */
PyObject *Matrix61c_add(Matrix61c* self, PyObject* args) {
    /* TODO: YOUR CODE HERE */

    //if (PyArg_UnpackTuple(args, "args", 1, 1, &mat)) {
    if (PyObject_TypeCheck(args, &Matrix61cType)) {

        Matrix61c* other = (Matrix61c *) args;
        matrix **newMat = (matrix **) malloc(sizeof(matrix*));
        int allocateSuccess = allocate_matrix(newMat, self->mat->rows, self->mat->cols);
        if (allocateSuccess == 0) {
            Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);

            rv->mat = *newMat;
            rv->shape = get_shape(rv->mat->rows, rv->mat->cols);
            // if (other->mat->rows != self->mat->rows || other->mat->cols != self->mat->cols) {
            // //throw ValueError
            // }
            add_matrix(rv->mat, self->mat, other->mat);
            return (PyObject *)rv;
        } else {
            deallocate_matrix(*newMat);
            return NULL;
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
        return NULL;
    }
}

/*
 * Substract the second numc.Matrix (Matrix61c) object from the first one. The first operand is
 * self, and the second operand can be obtained by casting `args`.
 */
PyObject *Matrix61c_sub(Matrix61c* self, PyObject* args) {
    /* TODO: YOUR CODE HERE */
    if (PyObject_TypeCheck(args, &Matrix61cType)) {

        Matrix61c* other = (Matrix61c *) args;
        matrix **newMat = (matrix **) malloc(sizeof(matrix*));
        int allocateSuccess = allocate_matrix(newMat, self->mat->rows, self->mat->cols);
        if (allocateSuccess == 0) {
            Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
            rv->mat = *newMat;
            rv->shape = get_shape(rv->mat->rows, rv->mat->cols);
            sub_matrix(rv->mat, self->mat, other->mat);
            return (PyObject *)rv;
        } else {
            deallocate_matrix(*newMat);
            return NULL;
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
        return NULL;
    }
}

/*
 * NOT element-wise multiplication. The first operand is self, and the second operand
 * can be obtained by casting `args`.
 */
PyObject *Matrix61c_multiply(Matrix61c* self, PyObject *args) {
    /* TODO: YOUR CODE HERE */
    if (PyObject_TypeCheck(args, &Matrix61cType)) {

        Matrix61c* other = (Matrix61c *) args;
        matrix **newMat = (matrix **) malloc(sizeof(matrix*));
        int allocateSuccess = allocate_matrix(newMat, self->mat->rows, self->mat->cols);
        if (allocateSuccess == 0) {
            Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
            rv->mat = *newMat;
            rv->shape = get_shape(rv->mat->rows, rv->mat->cols);
            mul_matrix(rv->mat, self->mat, other->mat);
            return (PyObject *)rv;
        } else {
            deallocate_matrix(*newMat);
            return NULL;
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
        return NULL;
    }
}

/*
 * Negates the given numc.Matrix.
 */
PyObject *Matrix61c_neg(Matrix61c* self) {
    /* TODO: YOUR CODE HERE */

    matrix **newMat = (matrix **) malloc(sizeof(matrix*));
    int allocateSuccess = allocate_matrix(newMat, self->mat->rows, self->mat->cols);
    if (allocateSuccess == 0) {
        Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
        rv->mat = *newMat;
        rv->shape = get_shape(rv->mat->rows, rv->mat->cols);
        neg_matrix(rv->mat, self->mat);
        return (PyObject *)rv;
    } else {
        deallocate_matrix(*newMat);
        return NULL;
    }


}

/*
 * Take the element-wise absolute value of this numc.Matrix.
 */
PyObject *Matrix61c_abs(Matrix61c *self) {
    /* TODO: YOUR CODE HERE */

    matrix **newMat = (matrix **) malloc(sizeof(matrix*));
    int allocateSuccess = allocate_matrix(newMat, self->mat->rows, self->mat->cols);
    if (allocateSuccess == 0) {
        Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
        rv->mat = *newMat;
        rv->shape = get_shape(rv->mat->rows, rv->mat->cols);
        abs_matrix(rv->mat, self->mat);
        return (PyObject *)rv;
    } else {
        deallocate_matrix(*newMat);
        return NULL;
    }
}

/*
 * Raise numc.Matrix (Matrix61c) to the `pow`th power. You can ignore the argument `optional`.
 */
PyObject *Matrix61c_pow(Matrix61c *self, PyObject *pow, PyObject *optional) {
    /* TODO: YOUR CODE HERE */

    if (PyLong_Check(pow)) {

        int toPow = (int)PyLong_AsLong(pow);
        matrix **newMat = (matrix **) malloc(sizeof(matrix*));
        int allocateSuccess = allocate_matrix(newMat, self->mat->rows, self->mat->cols);
        if (allocateSuccess == 0) {
            Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
            rv->mat = *newMat;
            rv->shape = get_shape(rv->mat->rows, rv->mat->cols);
            pow_matrix(rv->mat, self->mat, toPow);
            return (PyObject *)rv;
        } else {
            deallocate_matrix(*newMat);
            return NULL;
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "Argument must be an integer!");
        return NULL;
    }
}

/*
 * Create a PyNumberMethods struct for overloading operators with all the number methods you have
 * define. Ybinaryfuncou might find this link helpful: https://docs.python.org/3.6/c-api/typeobj.html
 */
PyNumberMethods Matrix61c_as_number = {
    /* TODO: YOUR CODE HERE */
    .nb_add = (binaryfunc) Matrix61c_add,
    .nb_subtract = (binaryfunc) Matrix61c_sub,
    .nb_multiply = (binaryfunc) Matrix61c_multiply,
    .nb_absolute = (unaryfunc) Matrix61c_abs,
    .nb_negative = (unaryfunc) Matrix61c_neg,
    .nb_power = (ternaryfunc) Matrix61c_pow

};


/* INSTANCE METHODS */

/*
 * Given a numc.Matrix self, parse `args` to (int) row, (int) col, and (double/int) val.
 * Return None in Python (this is different from returning null).
 */
PyObject *Matrix61c_set_value(Matrix61c *self, PyObject* args) {
    /* TODO: YOUR CODE HERE */
    if (PyTuple_Size(args) == 3) {
        PyObject *rows = NULL;
        PyObject *cols = NULL;
        PyObject *val = NULL;
        PyArg_UnpackTuple(args, "args", 3, 3, &rows, &cols, &val);
        if (rows && cols && val && PyLong_Check(rows) && PyLong_Check(cols) && (PyLong_Check(val)
                || PyFloat_Check(val))) {
            int toRows = (int)PyLong_AsLong(rows);
            int toCols = (int)PyLong_AsLong(cols);
            double toVal;
            if (PyLong_Check(val)) {
                toVal = (double)PyLong_AsDouble(val);
            } else {
                toVal = (double)PyFloat_AsDouble(val);
            }

            if (!(toRows >= self->mat->rows || toCols >= self->mat->cols)) {
                set(self->mat, toRows, toCols, toVal);
                Py_RETURN_NONE;
            } else {
                PyErr_SetString(PyExc_IndexError, "Specified row or column is out of range!");
                Py_RETURN_NONE;
            }

        } else {
            PyErr_SetString(PyExc_TypeError, "Wrong types for arguments!");
            Py_RETURN_NONE;
        }

    } else {
        PyErr_SetString(PyExc_TypeError, "Wrong number of arguments!");
        Py_RETURN_NONE;
    }
}

/*
 * Given a numc.Matrix `self`, parse `args` to (int) row and (int) col.
 * Return the value at the `row`th row and `col`th column, which is a Python
 * float/int.
 */
PyObject *Matrix61c_get_value(Matrix61c *self, PyObject* args) {
    /* TODO: YOUR CODE HERE */
    if (PyTuple_Size(args) == 2) {
        PyObject *rows = NULL;
        PyObject *cols = NULL;
        PyArg_UnpackTuple(args, "args", 2, 2, &rows, &cols);
        if (rows && cols && PyLong_Check(rows) && PyLong_Check(cols)) {
            int toRows = (int)PyLong_AsLong(rows);
            int toCols = (int)PyLong_AsLong(cols);
            //printf("%d ----- %d\n", toRows, self->mat->rows);
            if (!(toRows >= self->mat->rows || toCols >= self->mat->cols)) {

                double val = get(self->mat, toRows, toCols);

                PyObject *result = PyFloat_FromDouble(val);
                return result;

            } else {
                PyErr_SetString(PyExc_IndexError, "Specified row or column is out of range!");
                return NULL;
            }

        } else {
            PyErr_SetString(PyExc_TypeError, "Wrong types for arguments!");
            return NULL;
        }

    } else {
        PyErr_SetString(PyExc_TypeError, "Wrong number of arguments!");
        return NULL;
    }
}

/*
 * Create an array of PyMethodDef structs to hold the instance methods.
 * Name the python function corresponding to Matrix61c_get_value as "get" and Matrix61c_set_value
 * as "set"
 * You might find this link helpful: https://docs.python.org/3.6/c-api/structures.html
 */
PyMethodDef Matrix61c_methods[] = {
    /* TODO: YOUR CODE HERE */
    //{"set", (PyCFunction)(*Matrix61c_set_value), 4, "docstring"},
    {"get", (PyCFunction)(&Matrix61c_get_value), METH_VARARGS, "matrix getter"},
    {"set", (PyCFunction)(&Matrix61c_set_value), METH_VARARGS, "matrix setter"},
    {NULL, NULL, 0, NULL}
};

/* INDEXING */

/*
 * Given a numc.Matrix `self`, index into it with `key`. Return the indexed result.
 */
PyObject *Matrix61c_subscript(Matrix61c* self, PyObject* key) {
    /* TODO: YOUR CODE HERE */
    if (self->mat->rows == 1 || self->mat->cols == 1) {
        //1d matrix
        int length;
        if (self->mat->rows == 1) {
            length = self->mat->cols;
        } else {
            length = self->mat->rows;
        }
        if (PyLong_Check(key)) {
            int index = (int)PyLong_AsLong(key);
            if (self->mat->rows == 1) {
                if (index >= self->mat->cols || index < 0) {
                    PyErr_SetString(PyExc_IndexError, "Index out of range!");
                    return NULL;
                }
                return PyFloat_FromDouble(self->mat->data[0][index]);
            } else {
                if (index >= self->mat->rows || index < 0) {
                    PyErr_SetString(PyExc_IndexError, "Index out of range!");
                    return NULL;
                }
                return PyFloat_FromDouble(self->mat->data[index][0]);
            }
        } else if (PySlice_Check(key)) {
            Py_ssize_t start = 0;
            Py_ssize_t end = 0;
            Py_ssize_t step = 0;
            Py_ssize_t sliceLength = 0;
            int success = PySlice_GetIndicesEx(key, length, &start, &end, &step, &sliceLength);
            if (end - start == 0 || step != 1 || start > end) {
                PyErr_SetString(PyExc_ValueError, "Slice info not valid!");
                return NULL;
            }
            if (success == 0) {
                if (end - start == 1) {
                    if (self->mat->rows == 1) {
                        return PyFloat_FromDouble(self->mat->data[0][start]);
                    } else {
                        return PyFloat_FromDouble(self->mat->data[start][0]);
                    }
                }
                matrix **newMat = (matrix **) malloc(sizeof(matrix*));
                Matrix61c *rv = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
                int allRefSuccess;
                if (self->mat->rows == 1) {
                    allRefSuccess = allocate_matrix_ref(newMat, self->mat, 0, start, 1, end - start);
                } else {
                    allRefSuccess = allocate_matrix_ref(newMat, self->mat, start, 0, end - start, 1);
                }
                if (allRefSuccess != 0) {
                    //RuntimeError
                    Matrix61c_dealloc(rv);
                    deallocate_matrix(*newMat);
                    return NULL;
                }
                rv->mat = *newMat;
                rv->shape = get_shape(rv->mat->rows, rv->mat->cols);
                return (PyObject *)rv;
            } else {
                PyErr_SetString(PyExc_TypeError, "PySlice_GetIndicesEx could not parse key!");
                return NULL;
            }
        } else {
            PyErr_SetString(PyExc_TypeError, "Key is not an integer or a slice for 1D matrix!");
            return NULL;
        }
    } else {
        int length = self->mat->rows;
        if (PyLong_Check(key)) {
            int selectedRow = (int)PyLong_AsLong(key);
            if (selectedRow >= self->mat->rows || selectedRow < 0) {
                PyErr_SetString(PyExc_IndexError, "Index out of range!");
                return NULL;
            }
            matrix **newMat2D = (matrix **) malloc(sizeof(matrix*));
            Matrix61c *rv2 = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
            int integer2DSuccess = allocate_matrix_ref(newMat2D, self->mat, selectedRow, 0, 1, self->mat->cols);
            if (integer2DSuccess != 0) {
                //runtime error
                Matrix61c_dealloc(rv2);
                deallocate_matrix(*newMat2D);
                return NULL;
            }
            rv2->mat = *newMat2D;
            rv2->shape = get_shape(rv2->mat->rows, rv2->mat->cols);
            return (PyObject *)rv2;
        } else if (PySlice_Check(key)) {
            Py_ssize_t start2Dslice = 0;
            Py_ssize_t end2Dslice = 0;
            Py_ssize_t step2Dslice = 0;
            Py_ssize_t sliceLength2Dslice = 0;
            int success2Dslice = PySlice_GetIndicesEx(key, length, &start2Dslice, &end2Dslice, &step2Dslice, &sliceLength2Dslice);
            if (end2Dslice - start2Dslice == 0 || step2Dslice != 1 || start2Dslice > end2Dslice) {
                PyErr_SetString(PyExc_ValueError, "Slice info not valid!");
                return NULL;
            }
            if (success2Dslice == 0) {
                matrix **newMat2D = (matrix **) malloc(sizeof(matrix*));
                Matrix61c *rv2 = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
                int allRefSuccess2D;
                allRefSuccess2D = allocate_matrix_ref(newMat2D, self->mat, start2Dslice, 0, end2Dslice - start2Dslice, self->mat->cols);
                if (allRefSuccess2D != 0) {
                    //RuntimeError
                    Matrix61c_dealloc(rv2);
                    deallocate_matrix(*newMat2D);
                    return NULL;
                }
                rv2->mat = *newMat2D;
                rv2->shape = get_shape(rv2->mat->rows, rv2->mat->cols);
                return (PyObject *)rv2;
            } else {

                PyErr_SetString(PyExc_TypeError, "PySlice_GetIndicesEx could not parse key!");
                return NULL;
            }

        } else if (PyTuple_Check(key)) {
            if (PyTuple_Size(key) == 2) {
                PyObject *rows = NULL;
                PyObject *cols = NULL;

                if (PyArg_UnpackTuple(key, "args", 2, 2, &rows, &cols)) {
                    if (PyLong_Check(rows)) {
                        int rowIndex = (int)PyLong_AsLong(rows);
                        if (rowIndex >= self->mat->rows || rowIndex < 0) {
                            PyErr_SetString(PyExc_IndexError, "Index out of range!");
                            return NULL;
                        }
                        if (PyLong_Check(cols)) {
                            int colIndex = (int)PyLong_AsLong(cols);

                            if (colIndex >= self->mat->cols || colIndex < 0) {
                                PyErr_SetString(PyExc_IndexError, "Index out of range!");
                                return NULL;
                            }
                            return PyFloat_FromDouble(self->mat->data[rowIndex][colIndex]);
                        } else if (PySlice_Check(cols)) {
                            Py_ssize_t start2Dtuple = 0;
                            Py_ssize_t end2Dtuple = 0;
                            Py_ssize_t step2Dtuple = 0;
                            Py_ssize_t sliceLength2Dtuple = 0;
                            int success2Dtuple = PySlice_GetIndicesEx(cols, length, &start2Dtuple, &end2Dtuple, &step2Dtuple, &sliceLength2Dtuple);
                            if (end2Dtuple - start2Dtuple == 0 || step2Dtuple != 1 || start2Dtuple > end2Dtuple) {
                                PyErr_SetString(PyExc_ValueError, "Slice info not valid!");
                                return NULL;
                            }
                            if (success2Dtuple == 0) {
                                if (end2Dtuple - start2Dtuple == 1) {
                                    return PyFloat_FromDouble(self->mat->data[rowIndex][start2Dtuple]);
                                }
                                matrix **newMat2D = (matrix **) malloc(sizeof(matrix*));
                                Matrix61c *rv2 = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
                                int allRefSuccess2D;
                                allRefSuccess2D = allocate_matrix_ref(newMat2D, self->mat, rowIndex, start2Dtuple, 1, end2Dtuple - start2Dtuple);
                                if (allRefSuccess2D != 0) {
                                    //RuntimeError
                                    Matrix61c_dealloc(rv2);
                                    deallocate_matrix(*newMat2D);
                                    return NULL;
                                }
                                rv2->mat = *newMat2D;
                                rv2->shape = get_shape(rv2->mat->rows, rv2->mat->cols);
                                return (PyObject *)rv2;
                            } else {
                                PyErr_SetString(PyExc_TypeError, "PySlice_GetIndicesEx could not parse key!");
                                return NULL;
                            }
                        } else {
                            PyErr_SetString(PyExc_TypeError, "Key is not an integer, slice or a tuple for 2D matrix!");
                            return NULL;
                        }
                    } else if (PyLong_Check(cols)) {
                        int colIndex = (int)PyLong_AsLong(cols);
                        if (colIndex >= self->mat->cols || colIndex < 0) {
                            PyErr_SetString(PyExc_IndexError, "Index out of range!");
                            return NULL;
                        }
                        if (PySlice_Check(rows)) {
                            Py_ssize_t start2Dtuple = 0;
                            Py_ssize_t end2Dtuple = 0;
                            Py_ssize_t step2Dtuple = 0;
                            Py_ssize_t sliceLength2Dtuple = 0;
                            int success2Dtuple = PySlice_GetIndicesEx(rows, length, &start2Dtuple, &end2Dtuple, &step2Dtuple, &sliceLength2Dtuple);
                            if (end2Dtuple - start2Dtuple == 0 || step2Dtuple != 1 || start2Dtuple > end2Dtuple) {
                                PyErr_SetString(PyExc_TypeError, "slice info not valid!");
                                return NULL;
                            }
                            if (success2Dtuple == 0) {
                                if (end2Dtuple - start2Dtuple == 1) {
                                    return PyFloat_FromDouble(self->mat->data[start2Dtuple][colIndex]);
                                }
                                matrix **newMat2D = (matrix **) malloc(sizeof(matrix*));
                                Matrix61c *rv2 = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
                                int allRefSuccess2D;
                                allRefSuccess2D = allocate_matrix_ref(newMat2D, self->mat, start2Dtuple, colIndex, end2Dtuple - start2Dtuple, 1);
                                if (allRefSuccess2D != 0) {
                                    //RuntimeError
                                    Matrix61c_dealloc(rv2);
                                    deallocate_matrix(*newMat2D);
                                    return NULL;
                                }
                                rv2->mat = *newMat2D;
                                rv2->shape = get_shape(rv2->mat->rows, rv2->mat->cols);
                                return (PyObject *)rv2;
                            } else {
                                PyErr_SetString(PyExc_TypeError, "PySlice_GetIndicesEx could not parse key!");
                                return NULL;
                            }
                        } else {
                            PyErr_SetString(PyExc_TypeError, "Key is not an integer, slice or a tuple for 2D matrix!");
                            return NULL;
                        }
                    } else if (PySlice_Check(cols) && PySlice_Check(rows)) {
                        Py_ssize_t start2Dtuple1 = 0;
                        Py_ssize_t end2Dtuple1 = 0;
                        Py_ssize_t step2Dtuple1 = 0;
                        Py_ssize_t sliceLength2Dtuple1 = 0;
                        int success2Dtuple1 = PySlice_GetIndicesEx(rows, length, &start2Dtuple1, &end2Dtuple1, &step2Dtuple1, &sliceLength2Dtuple1);
                        Py_ssize_t start2Dtuple2 = 0;
                        Py_ssize_t end2Dtuple2 = 0;
                        Py_ssize_t step2Dtuple2 = 0;
                        Py_ssize_t sliceLength2Dtuple2 = 0;
                        int success2Dtuple2 = PySlice_GetIndicesEx(cols, length, &start2Dtuple2, &end2Dtuple2, &step2Dtuple2, &sliceLength2Dtuple2);
                        if (end2Dtuple1 - start2Dtuple1 == 0 || end2Dtuple2 - start2Dtuple2 == 0 || step2Dtuple1 != 1 || step2Dtuple2 != 1 
                            || start2Dtuple1 > end2Dtuple1 || start2Dtuple2 > end2Dtuple2) {
                            PyErr_SetString(PyExc_ValueError, "Slice info not valid!");
                            return NULL;
                        }
                        if (end2Dtuple1 - start2Dtuple1 == 1 && end2Dtuple2 - start2Dtuple2 == 1) {
                            return PyFloat_FromDouble(self->mat->data[start2Dtuple1][start2Dtuple2]);
                        }
                        if (success2Dtuple1 == 0 && success2Dtuple2 == 0) {

                            matrix **newMat2D = (matrix **) malloc(sizeof(matrix*));
                            Matrix61c *rv2 = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
                            int allRefSuccess2D;
                            allRefSuccess2D = allocate_matrix_ref(newMat2D, self->mat, start2Dtuple1, start2Dtuple2, end2Dtuple1 - start2Dtuple1, end2Dtuple2 - start2Dtuple2);
                            if (allRefSuccess2D != 0) {
                                //RuntimeError
                                Matrix61c_dealloc(rv2);
                                deallocate_matrix(*newMat2D);
                                return NULL;
                            }
                            rv2->mat = *newMat2D;
                            rv2->shape = get_shape(rv2->mat->rows, rv2->mat->cols);
                            return (PyObject *)rv2;
                        } else {
                            PyErr_SetString(PyExc_TypeError, "PySlice_GetIndicesEx could not parse key!");
                            return NULL;
                        }
                    } else {
                        PyErr_SetString(PyExc_TypeError, "Key is not an integer, slice or a tuple for 2D matrix!");
                        return NULL;
                    }
                } else {
                    PyErr_SetString(PyExc_TypeError, "Invalid arguments for 2D tuple!");
                    return NULL;
                }
            } else {
                PyErr_SetString(PyExc_TypeError, "Tuple is not of size 2!");
                return NULL;
            }
        } else {
            PyErr_SetString(PyExc_TypeError, "Key is not an integer, slice or a tuple for 2D matrix!");
            return NULL;
        }
        //2d matrix
    }
}



/*
 * Given a numc.Matrix `self`, index into it with `key`, and set the indexed result to `v`.
 */
int Matrix61c_set_subscript(Matrix61c * self, PyObject * key, PyObject * v) {
    /* TODO: YOUR CODE HERE */

    PyObject *subscripted = Matrix61c_subscript(self, key);
    if (subscripted != NULL) {
        //Matrix61c *rv = (Matrix61c*)subscripted;
        //int rows = rv->mat->rows;
        //int cols = rv->mat->cols;
        //int i, j, count;


        if (PyFloat_Check(subscripted)) {
            if (PyLong_Check(key)) {
                int index = (int)PyLong_AsLong(key);
                double value;
                if (!(PyLong_Check(v) || PyFloat_Check(v))) {
                    PyErr_SetString(PyExc_ValueError, "The Value is not of type integer or float!");
                    return -1;
                }
                if (PyLong_Check(v)) {
                    value = (double)PyLong_AsDouble(v);
                } else {
                    value = (double)PyFloat_AsDouble(v);
                }

                if (self->mat->rows == 1) {
                    set(self->mat, 0, index, value);
                    return 0;
                } else {
                    set(self->mat, index, 0, value);
                    return 0;
                }
            } else {
                PyErr_SetString(PyExc_TypeError, "The Key is not valid!");
                return -1;
            }
        // } else if ((Matrix61c*)) {
        //    return init_1d(subscripted, rows, cols, v);
        } else {
            //if (!(PyList_Check(v)) || PyList_Size(v) != rv->mat->rows * rv->mat->cols) {
	    //    return -1;
	    //}
	    Matrix61c *rv = (Matrix61c*)subscripted;
	    int rows = rv->mat->rows;
	    int cols = rv->mat->cols;
	    int size = PyList_Size(v);
	    int i, j;

//	    if (!(PyList_Check(v)) || PyList_Size(v) != rows * cols) {
//		    return -1;
//	    }
	    
	    if (rows == 1 || cols == 1) {
		
	   	if (!(PyList_Check(v)) || PyList_Size(v) != rows * cols) {
		       return -1;
		}	       
		    
		    
		int w;	
		for (w = 0; w < size; w++) {
			//int x = PyList_GetItem(v ,w);
			if (!PyLong_Check(PyList_GetItem(v, w)) && !PyFloat_Check(PyList_GetItem(v, w))) {
				return -1;
			}
		}
		int count = 0;
	    	for (i = 0; i < rows; i++) {
		    for (j = 0; j < cols; j++) {
			    set(rv->mat, i, j, PyFloat_AsDouble(PyList_GetItem(v, count)));
		    	    count++;
		    }
		}
		return 0;
	    } else {
		
		int y;
		for (y = 0; y < size; y++) {
			if (!PyList_Check(PyList_GetItem(v, y))) {
					return -1;
			} else {
				int x;
				int innerSize = PyList_Size(PyList_GetItem(v, y));
				for (x = 0; x < innerSize; x++) {
				       if (!(PyLong_Check(PyList_GetItem(PyList_GetItem(v, y), x))) && !(PyFloat_Check(PyList_GetItem(PyList_GetItem(v, y), x)))) {
						       return -1;
						       }	       
				}
			}
		}

		int outer, inner;
		for (outer = 0; outer < rows; outer++) {
			for (inner = 0; inner < cols; inner++) {
				set(rv->mat, outer, inner, PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(v, outer), inner)));
			



			}
		}	
		return 0;
	    }
           // return init_2d(subscripted, v);
        }
    } else {
        return -1;
    }

}

PyMappingMethods Matrix61c_mapping = {
    NULL,
    (binaryfunc) Matrix61c_subscript,
    (objobjargproc) Matrix61c_set_subscript,
};

/* INSTANCE ATTRIBUTES*/
PyMemberDef Matrix61c_members[] = {
    {
        "shape", T_OBJECT_EX, offsetof(Matrix61c, shape), 0,
        "(rows, cols)"
    },
    {NULL}  /* Sentinel */
};

PyTypeObject Matrix61cType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "numc.Matrix",
    .tp_basicsize = sizeof(Matrix61c),
    .tp_dealloc = (destructor)Matrix61c_dealloc,
    .tp_repr = (reprfunc)Matrix61c_repr,
    .tp_as_number = &Matrix61c_as_number,
    .tp_flags = Py_TPFLAGS_DEFAULT |
    Py_TPFLAGS_BASETYPE,
    .tp_doc = "numc.Matrix objects",
    .tp_methods = Matrix61c_methods,
    .tp_members = Matrix61c_members,
    .tp_as_mapping = &Matrix61c_mapping,
    .tp_init = (initproc)Matrix61c_init,
    .tp_new = Matrix61c_new
};


struct PyModuleDef numcmodule = {
    PyModuleDef_HEAD_INIT,
    "numc",
    "Numc matrix operations",
    -1,
    Matrix61c_class_methods
};

/* Initialize the numc module */
PyMODINIT_FUNC PyInit_numc(void) {
    PyObject* m;

    if (PyType_Ready(&Matrix61cType) < 0)
        return NULL;

    m = PyModule_Create(&numcmodule);
    if (m == NULL)
        return NULL;

    Py_INCREF(&Matrix61cType);
    PyModule_AddObject(m, "Matrix", (PyObject *)&Matrix61cType);
    printf("CS61C Fall 2020 Project 4: numc imported!\n");
    fflush(stdout);
    return m;
}
