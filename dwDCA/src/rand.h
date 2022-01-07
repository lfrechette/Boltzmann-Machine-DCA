#ifndef RAND_H
#define RAND_H

#include <gsl/gsl_rng.h>

extern gsl_rng *rg;

void init_rng(long unsigned int seed);

#endif
