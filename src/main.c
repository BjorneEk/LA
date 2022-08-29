
////////////////////////////////////////////////////////////////////////////
///        @author Gustaf Franzén :: https://github.com/BjorneEk;        ///
////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include "LA/matrix.h"

int main(int argc, char const *argv[]) {
  LA_matrix *m1, *m2, *p, *l, *u, *inv, *c;

  m1 = LA_readmat("matrix3.txt");
  m2 = LA_readmat("matrix2.txt");

  printf("read matricies.\n");
  printf("matrix 1:\n");
  LA_printmat(m1);
  printf("matrix 2:\n");
  LA_printmat(m2);

  printf("taking the second column of m2\n");
  c = LA_colat(m2, 1);
  LA_printmat(c);
  LA_printmat(m2);


  printf("attempting lu decomposition on m2\n");
  int res;
  res = LA_lup_decompose(m2, &p, &l, &u);
  printf("p:\n");
  LA_printmat(p);
  printf("\nl:\n");
  LA_printmat(l);
  printf("\nu:\n");
  LA_printmat(u);

  printf("attempting invert m2\n");
  inv = LA_inverse(m2);
  LA_printmat(inv);
  LA_printmat(m2);

  LA_free(m1);
  LA_free(m2);
  LA_free(p);
  LA_free(l);
  LA_free(u);
  LA_free(inv);
  LA_free(c);
  return 0;
}
