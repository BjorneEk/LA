/*==========================================================*
 *
 * @author Gustaf Franzén :: https://github.com/BjorneEk;
 *
 * matrix header file for my linear algebra library
 *
 * functions with a name ending in _r  aswell as the row
 * and column swap functions alter the input matrix
 * and other functions return a new matrix with the result.
 *
 *==========================================================*/

#ifndef _LA_MATRIX_H_
#define _LA_MATRIX_H_


#include <UL/types.h>


#define LA_MIN_COEF 0.000000000000001
#define MIN_MAT_THRESHOLD 40

typedef struct LA_matrix_s {
        /**
         * the height and width of the matrix;
         **/
        u32_t rows, cols;

        /**
        * the matrix data stored in linear continous memory for faster acces
        **/
        f64_t * data;
} LA_matrix;


/*==========================================================*
 *  allocates a new matrix of dimension width x height and
 *  setts all elements to zero.
 *==========================================================*/
LA_matrix * LA_empty(u32_t cols, u32_t rows);

/*==========================================================*
 *  allocates a new matrix of dimension width x height and
 *  setts the elements to the suplied values
 *==========================================================*/
LA_matrix * LA_new(u32_t cols, u32_t rows, ...);

/*==========================================================*
 *  creates a new copy of the input matrix
 *==========================================================*/
LA_matrix * LA_copy(const LA_matrix * const m);

/*==========================================================*
 *  frees the memory allocated to the matrix m
 *==========================================================*/
void        LA_free(LA_matrix * m);

/*==========================================================*
 *  prints m to the console
 *==========================================================*/
void        LA_printmat(const LA_matrix * m);

/*==========================================================*
 *  prints m to the console verbously
 *==========================================================*/
void        LA_debugmat(const LA_matrix * m);

/*==========================================================*
 *  reads a matrix from a file according to the following
 *  format: width space height newline then the matrix
 *  with cols separated by spaces and rows by newlines
 *==========================================================*/
LA_matrix * LA_readmat(const char * filename);

/*==========================================================*
 *  reads a array of matricies from a file
 *==========================================================*/
LA_matrix * LA_read_matricies(const char * filename, i32_t * len);

/*==========================================================*
 *  writes a matrix to a file
 *==========================================================*/
void        LA_writemat(const LA_matrix * m,
                        const char * filename);

/*==========================================================*
 *  sets the diagonal in m to d
 *==========================================================*/
void LA_set_diagonal(LA_matrix *m, f64_t d);

/*==========================================================*
 *  creates a nxn identity matrix
 *==========================================================*/
LA_matrix * LA_eye(u32_t n);

/*==========================================================*
 *  creates a nxm random matrix
 *==========================================================*/
LA_matrix * LA_random(u32_t width, u32_t height,
                      f64_t min,   f64_t max);

/*==========================================================*
 *  creates a transpose the matrix m
 *==========================================================*/
LA_matrix * LA_transpose(const LA_matrix * m);

/*==========================================================*
 *  multiplies each element in m with s
 *==========================================================*/
void        LA_scale_r(LA_matrix * m, f64_t s);
/*==========================================================*
 *  creates a new matrix with the values of m scaled by s
 *==========================================================*/
LA_matrix * LA_scale(const LA_matrix * m, f64_t s);

/*==========================================================*
 *  add the matrix b to the matrix a
 *==========================================================*/
void        LA_add_r(LA_matrix * a, const LA_matrix * b);
/*==========================================================*
 *  creates a new matrix that is equal to a + b
 *==========================================================*/
LA_matrix * LA_add(const LA_matrix * a, const LA_matrix * b);

/*==========================================================*
 *  subtract the matrix b to the matrix a
 *==========================================================*/
void        LA_sub_r(LA_matrix * a, const LA_matrix * b);
/*==========================================================*
 *  creates a new matrix that is equal to a - b
 *==========================================================*/
LA_matrix * LA_sub(const LA_matrix * a, const LA_matrix * b);

/*==========================================================*
 *  multiply matrix a by matrix b, store result in a
 *==========================================================*/
void        LA_mult_r(LA_matrix * a, const LA_matrix * b);

LA_matrix * LA_mult_naive(const LA_matrix * a, const LA_matrix * b);
/*==========================================================*
 *  creates a new matrix that is equal to a * b
 *==========================================================*/
LA_matrix * LA_mult(const LA_matrix * a, const LA_matrix * b);

/*==========================================================*
 *  creates a new matrix that is equal to a * b using
 *  the BLAS function cblas_dgemm
 *==========================================================*/
LA_matrix * LA_multBLAS(const LA_matrix * a, const LA_matrix * b);

/*==========================================================*
 *  returns true is a is identical to b
 *==========================================================*/
bool        LA_equals(const LA_matrix * a, const LA_matrix * b);

/*==========================================================*
 * returns true is the width and height of a is equal to
 * the width and height of b
 *==========================================================*/
bool        LA_equal_dim(const LA_matrix * a, const LA_matrix * b);

/*==========================================================*
 *  swap row _p with row _q in the matrix m
 *==========================================================*/
void        LA_rswap(LA_matrix * m, i32_t _p, i32_t _q);

/*==========================================================*
 *  swap column _p with column _q in the matrix m
 *==========================================================*/
void        LA_cswap(LA_matrix * m, i32_t _p, i32_t _q);

/*==========================================================*
 *  returns the row r of a matrix m
 *==========================================================*/
LA_matrix * LA_rowat(const LA_matrix * m, u32_t r);

/*==========================================================*
 *  returns the col c of a matrix m
 *==========================================================*/
LA_matrix * LA_colat(const LA_matrix * m, u32_t c);

/*==========================================================*
 *  wultiply al elements in row r with s
 *==========================================================*/
void        LA_mult_row_r(LA_matrix * m, u32_t r, f64_t s);

/*==========================================================*
 *  create a new copy of a where row r has been
 *  multiplied by s
 *==========================================================*/
LA_matrix * LA_mult_row(const LA_matrix * m, u32_t r, f64_t s);

/*==========================================================*
 *  wultiply al elements in col c with s
 *==========================================================*/
void        LA_mult_col_r(LA_matrix * m, u32_t c, f64_t s);

/*==========================================================*
 *  create a new copy of a where col c has been
 *  multiplied by s
 *==========================================================*/
LA_matrix * LA_mult_col(const LA_matrix * m, u32_t c, f64_t s);

/*==========================================================*
 * add row r multiplied with the scalar s to the row d
 *==========================================================*/
void        LA_row_add_row_r(LA_matrix * m, u32_t d, u32_t r, f64_t s);

/*==========================================================*
 * create a new copy of m where the row r multiplied with
 * the scalar s have been addet to the row d
 *==========================================================*/
LA_matrix * LA_row_add_row(LA_matrix * m, u32_t d, u32_t r, f64_t s);

u32_t LA_lup_decompose(LA_matrix *m, LA_matrix **P, LA_matrix **L, LA_matrix **U);

/*==========================================================*
 * linear system solve forward
 *
 * L is a lower triangular matrix
 *
 * b is a column vector ie. nx1 matrix
 * returns the solution in a column vector
 *==========================================================*/
LA_matrix * LA_ls_solvefwd(LA_matrix *L, LA_matrix *b);

/*==========================================================*
 * linear system solve backwards
 *
 * U is a upper triangular matrix
 *
 * b is a column vector ie. nx1 matrix
 * returns the solution in a column vector
 *==========================================================*/
LA_matrix * LA_ls_solvebck(LA_matrix *U, LA_matrix *b);

/*==========================================================*
 * linear equation system solver that solves the equation
 *
 * AX = Y
 *
 * where A is a nxn matrix, and X and Y are column vectros
 *
 * returns X in a LA_matrix with width 1
 *==========================================================*/
LA_matrix * LA_solve(LA_matrix *A, LA_matrix * Y);

/*==========================================================*
 * linear equation system solver that solves the equation
 *
 * AX = Y
 *
 * using precalculated lower and upper triangupar matricies
 *==========================================================*/
LA_matrix * LA_solve_LU(LA_matrix *L, LA_matrix *U, LA_matrix *P, LA_matrix * Y);


/*==========================================================*
 * calculates the inverse of the matrix m
 *==========================================================*/
LA_matrix * LA_inverse(LA_matrix *m);

/*==========================================================*
 * calculates the determinant of the matrix m
 *==========================================================*/
f64_t LA_det(LA_matrix * m);

/*==========================================================*
 * calculates dot product of two vectors
 *==========================================================*/
f64_t LA_dot(const LA_matrix * a, const LA_matrix * b);
/*==========================================================*
 * calculates dot product of two columns in two matricies
 *==========================================================*/
f64_t LA_dot_col(const LA_matrix * a, u32_t ac, const LA_matrix * b, u32_t bc);

/*==========================================================*
 * calculates absolute value (length) of a vector
 *==========================================================*/
f64_t LA_abs(const LA_matrix * a);
/*==========================================================*
 * calculates absolute value (length) of a column
 *==========================================================*/
f64_t LA_abs_col(const LA_matrix * a, u32_t c);
/*==========================================================*
 * calculates absolute value (length) of each column
 * in a matrix and returns a 1 x 3 matrix with the result
 *==========================================================*/
LA_matrix * LA_abs_mat(const LA_matrix * a);


 #endif /*_LA_MATRIX_H_*/
