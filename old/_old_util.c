
#include "util.h"
#include <stdlib.h>
#include <stdio.h>


void assert(int assertion, char* message) {
        if (assertion) return;
        fprintf(stderr, "[\033[31;1;4mAssert error\033[0m]: %s\n", message);
        exit(1);
}

double random_d(double min, double max) {
        double d;
        d = (double) rand() / ((double) RAND_MAX + 1);
        return (min + d * (max - min));
}

char * read_file(const char * filename) {

        char *source = NULL;
        FILE *fp;

        fp = fopen(filename, "r");
        assert(fp != NULL, "could not read from file");
        /* Go to the end of the file. */
        if (fseek(fp, 0L, SEEK_END) == 0) {
                /* Get the size of the file. */
                long bufsize = ftell(fp);
                assert(bufsize != -1, "Error reading file");
                /* Allocate our buffer to that size. */
                source = malloc(sizeof(char) * (bufsize + 1));
                assert(source != NULL, "Out of memory");

                /* Go back to the start of the file. */
                assert(fseek(fp, 0L, SEEK_SET) == 0, "Error reading file");

                /* Read the entire file into memory. */
                size_t newLen = fread(source, sizeof(char), bufsize, fp);
                if ( ferror( fp ) != 0 ) {
                        fputs("Error reading file", stderr);
                } else {
                        source[newLen++] = '\0'; /* Just to be safe. */
                }
        }
        fclose(fp);
        return source;
}
