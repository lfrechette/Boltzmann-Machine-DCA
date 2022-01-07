#ifndef MC_IO_H
#define MC_IO_H

#include <stdlib.h>
#include <armadillo>

extern int N;
extern int q;

void read_dims(const char *name, const char *method);
void read_params(const char *name, const char *method, arma::mat &h, arma::cube &J);
void print_traj_scalar(double *s, int n, const char *dir, const char *name);
void pconf(int *seq, int time, const char *dir, const char *name);

#endif
