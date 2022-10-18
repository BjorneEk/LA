
////////////////////////////////////////////////////////////////////////////
///        @author Gustaf Franzén :: https://github.com/BjorneEk;        ///
////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include "../../LA/matrix.h"
#include "la-test.h"

void old_test()
{
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
        i32_t res;
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
}

int main(int argc, char const *argv[])
{
        LA_matrix * mats, *inv;
        i32_t len, i, det;
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

        f64_t dot = LA_dot(&mats[0], &mats[1]);

        printf("det = %f", dot);
        for(i = 0; i < len; i++) {
                free(mats[i].data);
        }
        free(mats);

        LA_test_mult_advanced(70);
        return 0;

}
