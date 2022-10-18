
#ifndef _UTIL_H
#define _UTIL_H

/*==========================================================*
 *
 * @author Gustaf Franzén :: https://github.com/BjorneEk;
 *
 * header for various utilities
 *
 *==========================================================*/

#define RAND_MAX 0x7fffffff

void assert(int assertion, char* message);


double random_d(double min, double max);

char * read_file(const char * filename);

 #endif /*_UTIL_H*/
