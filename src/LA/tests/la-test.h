

/*==========================================================*
 *
 * @author Gustaf Franzén :: https://github.com/BjorneEk;
 *
 * header file for tests to test the LA matrix and linear
 * algebra library
 *
 *==========================================================*/

#ifndef _LA_TEST_H_
#define _LA_TEST_H_

#include "../matrix.h"

void test(bool success, const char * msg);
void LA_test_copy(LA_matrix * mats);
void LA_test_equality(LA_matrix * mats);
void LA_test_transpose(LA_matrix * mats);
void LA_test_mult(LA_matrix * mats);
void LA_test_mult_advanced(int max);
void LA_test_gaus_jordan_LUP(LA_matrix * mats);
void LA_test_solve(LA_matrix * mats);
void LA_test_det(LA_matrix * mats);
#endif /* _LA_TEST_H_ */
