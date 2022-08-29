#include "matrix.h"
#include "util.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/*==========================================================*
 *
 * @author Gustaf Franzén :: https://github.com/BjorneEk;
 *
 *  functions related to and using the
*   LU decomposition of matricies
 *
 * this library uses LU Factorisation with partial pivoting
 * and solves the equation:
 *
 * P * A = L * U
 *
 * where L and U are lower and upper triangular matricieas,
 * P is a permutation matrix and A is the input matrix.
 *==========================================================*/


/*==========================================================*
 * Finds the maxid on the column (starting from k -> num_rows)
 * This method is used for pivoting in LUP decomposition
 *==========================================================*/
int _LA_absmaxr(LA_matrix *m, unsigned int k) {
  int i, max_idx;
  double max;

  max = m->data[k + k*m->cols];
  max_idx = k;
  for(i = k+1; i < m->rows; i++) {
    if (fabs(m->data[k+i*m->rows]) > max) {
      max = fabs(m->data[k+i*m->rows]);
      max_idx = i;
    }
  }
  return max_idx;
}

/*==========================================================*
 * LU Factorisation function with partial pivoting
 * that solves the equation:
 *
 * P * A = L * U
 *
 * using Gauss-jordan elimination
 * where L and U are lower and upper triangular matricieas,
 * P is a permutation matrix and A is the input matrix.
 * returns the number of permutations
 *==========================================================*/
unsigned int LA_lup_decompose(LA_matrix *m, LA_matrix **P, LA_matrix **L, LA_matrix **U) {
  unsigned int num_permutations;
  int i,j, pivot;
  double mult;

  assert(m->cols == m->rows, "Matrix must be square (LUP_decompose)");

  *L = LA_new_mat(m->cols, m->rows);
  *U = LA_copy(m);
  *P = LA_identity_mat(m->cols);
  num_permutations = 0;

  for(j = 0; j < (*U)->cols; j++) {
    // Retrieves the row with the biggest element for column (j)
    pivot = _LA_absmaxr(*U, j);
    assert(fabs((*U)->data[j+pivot*m->cols]) > LA_MIN_COEF, "cannot LU matrix degenerate");
    if (pivot!=j) {
      //Pivots LU and P accordingly to the rule
      LA_rswap(*U, j, pivot);
      LA_rswap(*L, j, pivot);
      LA_rswap(*P, j, pivot);
      // Keep the number of permutations to easily calculate the
      // determinant sign afterwards
      num_permutations++;
    }
    for(i = j+1; i < (*U)->rows; i++) {
      mult = (*U)->data[j+i*m->cols] / (*U)->data[j+j*m->cols];
      // Building the U upper rows
      LA_row_add_row_r(*U, i, j, -mult);
      // Store the multiplier in L
      (*L)->data[j+i*m->cols] = mult;
    }
  }
  LA_set_diagonal(*L, 1.0);
  return num_permutations;
}

/*==========================================================*
 * linear system solve forward
 *
 * L is a lower triangular matrix
 *
 * b is a column vector ie. nx1 matrix
 * returns the solution in a column vector
 *==========================================================*/
LA_matrix * LA_ls_solvefwd(LA_matrix *L, LA_matrix *b) {
  LA_matrix * x;
  int i,j;
  double tmp;

  assert(b->cols <= 1,"LA_ls_solvefwd input is not a column vector");
  x = LA_new_mat(1, L->cols);

  for(i = 0; i < L->cols; i++) {
    tmp = b->data[i*b->cols];
    for(j = 0; j < i ; j++)
      tmp -= L->data[j+i*L->cols] * x->data[j*x->cols];
    x->data[i*x->cols] = tmp / L->data[i+i*L->cols];
  }
  return x;
}

/*==========================================================*
 * linear system solve backwards
 *
 * U is a upper triangular matrix
 *
 * b is a column vector ie. nx1 matrix
 * returns the solution in a column vector
 *==========================================================*/
LA_matrix * LA_ls_solvebck(LA_matrix *U, LA_matrix *b) {
  LA_matrix * x;
  int i, j;
  double tmp;
  i = U->cols;
  x = LA_new_mat(1, U->cols);
  while(i-- > 0) {
    tmp = b->data[i*b->cols];
    for(j = i; j < U->cols; j++) {
      tmp -= U->data[j+i*U->cols] * x->data[j*x->cols];
    }
    x->data[i*x->cols] = tmp / U->data[i+i*U->cols];
  }
  return x;
}

/*==========================================================*
 * linear equation system solver that solves the equation
 *
 * AX = Y
 *
 * where A is a nxn matrix, and X and Y are column vectros
 *
 * returns X in a LA_matrix with cols 1
 *==========================================================*/
LA_matrix * LA_solve(LA_matrix *A, LA_matrix * Y) {
  LA_matrix *L, *U, *P, *PY, *b, *x;
  LA_lup_decompose(A, &P, &L, &U);
  assert(!(U->rows != Y->rows || Y->cols != 1),"Cannot solve system");
  PY = LA_mult(P, Y);

  // We solve L*y = P*b using forward substition
  b = LA_ls_solvefwd(L, PY);
  // We solve U*x=y
  x = LA_ls_solvebck(U, b);

  LA_free(b);
  LA_free(PY);
  LA_free(L);
  LA_free(U);
  LA_free(P);
  return x;
}

/*==========================================================*
 * linear equation system solver that solves the equation
 *
 * AX = Y
 *
 * using precalculated lower and upper triangupar matricies
 *==========================================================*/
LA_matrix * LA_solve_LU(LA_matrix *L, LA_matrix *U, LA_matrix *P, LA_matrix * Y) {
  LA_matrix *PY, *b, *x;
  assert(!(U->rows != Y->rows || Y->cols != 1),"Cannot solve system");
  PY = LA_mult(P, Y);
  // We solve L*y = P*b using forward substition
  b = LA_ls_solvefwd(L, PY);

  // We solve U*x=y
  x = LA_ls_solvebck(U, b);

  LA_free(b);
  LA_free(PY);
  return x;
}

/*==========================================================*
 * calculates the inverse of the matrix m
 *==========================================================*/
LA_matrix * LA_inverse(LA_matrix *m) {
  LA_matrix *r, *I, *invx, *Ix, *P, *U, *L;
  int i,j, n;
  LA_lup_decompose(m, &P, &L, &U);
  n = m->cols;
  r = LA_new_mat(n,n);
  I = LA_identity_mat(m->rows);

  for(j =0; j < n; j++) {
    Ix = LA_colat(I, j);
    invx = LA_solve_LU(L, U, P, Ix);
    for(i = 0; i < invx->rows; i++) {
      r->data[j+i*r->cols] = invx->data[i*invx->cols];
    }
    LA_free(invx);
    LA_free(Ix);
  }
  LA_free(I);
  LA_free(P);
  LA_free(U);
  LA_free(L);
  return r;
}

/*==========================================================*
 * calculates the determinant of the matrix m
 *==========================================================*/
double LA_det(LA_matrix * m) {
  LA_matrix *P, *U, *L;
  int sz, i;
  double res, s;

  assert(m->cols == m->rows, "Matrix must be square (calc determinant)");
  s = ((LA_lup_decompose(m, &P, &L, &U) % 2) == 0) ? 1.0 : -1.0;
  res = 1.0f;

  for(i = 0; i < U->rows; i++)
    res *= U->data[i+i*U->cols];
  LA_free(P);
  LA_free(L);
  LA_free(U);
  return res * s;
}
