/*==========================================================*
 *
 * @author Gustaf Franzén :: https://github.com/BjorneEk;
 *
 * functions for testing the LA matrix and linear
 * algebra library, function needs the matricies.txt file to
 * be read to a LA_matrix array and passed as the single
 * argument
 *==========================================================*/


#include "la-test.h"
#include <stdbool.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <UL/assert.h>
#include <UL/math.h>

void test(bool success, const char * msg) {
        if (!success)
                printf("Test: %s [\033[31;1;4mTest failed\033[0m]\n", msg);
        else
                printf("Test: %s [\033[32;1;4mTest succeded\033[0m]\n", msg);
}

void LA_test_copy(LA_matrix * mats) {
        LA_matrix * copy;
        copy = LA_copy(&mats[3]);
        test(LA_equals(copy, &mats[3]), "copy matrix");
        LA_free(copy);
}

void LA_test_equality(LA_matrix * mats) {
        LA_matrix * copy;

        copy = LA_copy(&mats[4]);
        test(LA_equals(&mats[3], &mats[3]), "equals itself       ");
        test(LA_equals(&mats[4], copy), "equal comparison    ");
        test(!LA_equals(&mats[4], &mats[3]), "inequal comparison 1");
        test(!LA_equals(&mats[0], &mats[2]), "inequal comparison 2");
        LA_free(copy);
}
void LA_test_transpose(LA_matrix * mats) {
        LA_matrix *t1, *t2;

        t1 = LA_transpose(&mats[4]);
        t2 = LA_transpose(&mats[0]);
        test(LA_equals(t1, &mats[5]), "test transpose 1");
        test(LA_equals(t2, &mats[2]), "test transpose 2");
        LA_free(t1);
        LA_free(t2);
}
void LA_test_mult(LA_matrix * mats) {
        LA_matrix *m1, *m2;

        m1 = LA_mult(&mats[4], &mats[0]);
        m2 = LA_mult(&mats[4], &mats[5]);
        LA_printmat(m1);
        LA_printmat(m2);
        test(LA_equals(m1, &mats[7]), "test multiplication 1");
        test(LA_equals(m2, &mats[8]), "test multiplication 2");
        LA_free(m1);
        LA_free(m2);
}

void LA_test_mult_advanced(i32_t max) {
        LA_matrix *m1, *m2, *res;
        FILE *fp_1, *fp_2, *fp_3;
        i64_t * times, * times_naive, *times_blas;
        clock_t begin, end;
        i32_t i, dim;

        times = malloc(max * sizeof(long));
        times_naive = malloc(max * sizeof(long));
        times_blas = malloc(max * sizeof(long));
        assertf(times != NULL && times_naive != NULL, "Out of memory (%s)", __func__);

        for (i = 1; i <= max; i++) {
                if(randf(0, 1) > 0.5) dim = i;
                else dim = (int) (randf(1, 400));
                m1 = LA_random(i,i,-500,500);
                m2 = LA_random(i,i,-500,500);

                begin = clock();
                res = LA_mult(m1, m2);
                end = clock();

                times[i-1] = (long)(1000000.0*(f64_t)(end - begin) / CLOCKS_PER_SEC);
                LA_free(res);

                begin = clock();
                res = LA_multBLAS(m1, m2);
                end = clock();

                times_blas[i-1] = (long)(1000000.0*(f64_t)(end - begin) / CLOCKS_PER_SEC);
                LA_free(res);

                begin = clock();
                res = LA_mult_naive(m1, m2);
                end = clock();

                times_naive[i-1] = (long)(1000000.0*(f64_t)(end - begin) / CLOCKS_PER_SEC);

                LA_free(m1);
                LA_free(m2);
                LA_free(res);
        }

        fp_1 = fopen("matmult_times_fast.csv", "w");
        fp_2 = fopen("matmult_times_naive.csv", "w");
        fp_3 = fopen("matmult_times_blas.csv", "w");
        for(i = 0; i < max; i++){
                if(i!= max-1){
                        fprintf(fp_1, "%ld,", times[i]);
                        fprintf(fp_2, "%ld,", times_naive[i]);
                        fprintf(fp_3, "%ld,", times_blas[i]);
                }
                else{
                        fprintf(fp_1, "%ld\n", times[i]);
                        fprintf(fp_2, "%ld\n", times_naive[i]);
                        fprintf(fp_3, "%ld\n", times_blas[i]);
                }
        }
        fclose(fp_1);
        fclose(fp_2);
        fclose(fp_3);
        free(times);
        free(times_naive);
        free(times_blas);
}

void LA_test_gaus_jordan_LUP(LA_matrix * mats) {
        LA_matrix *P, *L, *U;
        LA_lup_decompose(&mats[4], &P, &L, &U);
        printf("P:\n");
        LA_printmat(P);
        printf("L:\n");
        LA_printmat(L);
        printf("U:\n");
        LA_printmat(U);
        printf("\n");
}
void LA_test_solve(LA_matrix * mats) {
        LA_matrix *x;

        x = LA_solve(&mats[4], &mats[9]);
        LA_printmat(x);
        LA_printmat(&mats[10]);
        test(LA_equals(x, &mats[10]), "test equation solve");
        LA_free(x);
}

void LA_test_det(LA_matrix * mats) {

        test(LA_det(&mats[4]) == 64.0f, "test determinant");
}
