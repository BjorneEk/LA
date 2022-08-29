#include "tests/la-test.h"
#include "util.h"
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>


int main(int argc, char const *argv[]) {
  LA_matrix * mats, *inv;
  int len, i, det;
  mats = LA_read_matricies("matricies.txt", &len);
  printf("read %i matricies\n", len);
  for(i = 0; i < len; i++) {
    printf("Matrix %i:\n",i);
    LA_printmat(&mats[i]);
  }
  LA_test_copy(mats);
  LA_test_equality(mats);
  LA_test_transpose(mats);
  LA_test_mult(mats);
  LA_test_gaus_jordan_LUP(mats);
  LA_test_solve(mats);
  LA_test_det(mats);
  inv = LA_inverse(&mats[4]);
  LA_printmat(inv);

  double dot = LA_dot(&mats[0], &mats[1]);

  printf("det = %f", dot);
  for(i = 0; i < len; i++) {
    free(mats[i].data);
  }
  free(mats);

  LA_test_mult_advanced(700);
return 0;
}
