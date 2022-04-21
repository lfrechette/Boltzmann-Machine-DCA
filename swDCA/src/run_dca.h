#ifndef RUN_DCA_H
#define RUN_DCA_H

#include <armadillo>

extern gsl_rng *rg;

void init_default(model &mymodel, arma::mat &msa_freq);

void fit(model &mymodel, arma::mat &msa_freq, arma::cube &msa_corr, int nrep);

void read_inputs(std::string in_name);

void read_params(model &mymodel, std::string in_name);

void print_params(model &mymodel, std::string out_name);

void print_seqs(model &mymodel, std::string out_name);

arma::mat get_norm(model &mymodel);

#endif
