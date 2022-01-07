/* mc.cpp Created by Layne Frechette on Feb. 17, 2021
 *
 * Implements Monte Carlo moves for sampling under given parameters.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <armadillo>
#include <sstream>

//User-defined headers
#include "mc.h"
#include "msa.h"
#include "model.h"
#include "io.h"
#include "rand.h"

int nequil=1000;
std::string ga_str = "MEAVDANSLAQAKEAAIKELKQYGIGDYYIKLINNAKTVEGVESLKNEILKALPTE";
std::string gb_str = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE";

/***************/
/*** Methods ***/
/***************/

void run_mc_traj(model &model, int n){

  int dump_freq = n/model.num_seqs;
  int nseed=1;
  int nblock=10;
  model.avg_ene=0;

  model.mom1.zeros();
  model.mom2.zeros();
  model.mom1_ene1.zeros();
  model.mom1_ene2.zeros();
  model.mom2_ene1.zeros();
  model.mom2_ene2.zeros();
  model.mom1_err=0;
  model.mom2_err=0;

  std::vector<arma::mat> mom1_temp(nblock);
  std::vector<arma::cube> mom2_temp(nblock);
  for(int i=0; i<nblock; i++){
    mom1_temp[i] = arma::zeros(model.N,model.q);
    mom2_temp[i] = arma::zeros(model.q,model.q,model.N*(model.N-1)/2);
  }

  double *mom1_block = new double[nblock];
  double *mom2_block = new double[nblock];
  for(int i=0; i<nblock; i++){
    mom1_block[i]=0;
    mom2_block[i]=0;
  }
  int block_cnt=0;

  std::vector<int> seq(model.N,0.0);
  std::vector<int> seq2(model.N,0.0);
  std::vector<int> ga(model.N,0.0);
  std::vector<int> gb(model.N,0.0);

  for(int i=0; i<model.N; i++){
    ga[i] = letter_to_number(ga_str.at(i),model.q);
    gb[i] = letter_to_number(gb_str.at(i),model.q);
  }
  //Compute 1st and 2nd moments by sampling w/ Metropolis algorithm.
  //Could run the different seeds in parallel...
  for(int s=0; s<nseed; s++){
    
    if(model.mc_init=="random"){
      //Initialize a random sequence
      for(int i=0; i<model.N; i++){
        seq[i] = (int)gsl_rng_uniform_int(rg, (unsigned long int)model.q);
      }
    }
    else if(model.mc_init=="gagb"){
      for(int i=0; i<model.N; i++) seq[i] = ga[i];
    }

    //Equilibrate
    if(model.mc_init=="random"){
    for(int t=0;  t<nequil; t++){
      //std::cout << "Equilibrating...Pass " << t << ", energy=" << model.get_energy(seq) << std::endl;
      do_sweep(model, seq);
    }
    }

    //Production
    int dump_cnt=0;
    for(int t=0; t<n/nseed; t++){
      //std::cout << "Pass number " << t << std::endl;
      if(t==n/nseed/2 && model.mc_init=="gagb"){
        for(int i=0; i<model.N; i++) seq[i] = gb[i];
      }
      if(t%dump_freq==0){
        std::stringstream seqstr;
        for(int j=0; j<seq.size(); j++){
          seqstr << number_to_letter(seq[j],model.q);
        }
        model.seqs[dump_cnt] = seqstr.str();
        dump_cnt++;
      }
      do_sweep(model, seq);
      model.avg_ene += model.get_energy(seq); 
      //Update moments
      for(int i=0; i<model.N; i++){
        model.mom1(i,seq[i])++;
        mom1_temp[block_cnt](i,seq[i])++;
      }
      for(int i=0; i<(model.N-1); i++){
        for(int j=(i+1); j<model.N; j++){
          int index = (model.N-1)*i-i*(i+1)/2+j-1;
          model.mom2(seq[i],seq[j],index)++;
          mom2_temp[block_cnt](seq[i],seq[j],index)++;
          if(seq[i]!=seq[j]){
            model.mom2(seq[j],seq[i],index) += 0.5;//symmetrize
            mom2_temp[block_cnt](seq[j],seq[i],index) += 0.5;
            model.mom2(seq[i],seq[j],index) -= 0.5;
            mom2_temp[block_cnt](seq[i],seq[j],index) -= 0.5;
          }
        }
      }

      //Get energy-weighted frequencies if doing double-well
      if(model.nwell==2){
        double E1 = model.get_energy_single(seq, 1);
        double E2 = model.get_energy_single(seq, 2);
        double exp1 = exp(-E1/model.Tmix);
        double exp2 = exp(-E2/model.Tmix);
        double q1 = exp1/(exp1+exp2);
        double q2 = exp2/(exp1+exp2);
        for(int i=0; i<model.N; i++){
          model.mom1_ene1(i,seq[i]) += q1;
          model.mom1_ene2(i,seq[i]) += q2;
        }
        for(int i=0; i<(model.N-1); i++){
          for(int j=(i+1); j<model.N; j++){
            int index = (model.N-1)*i-i*(i+1)/2+j-1;
            model.mom2_ene1(seq[i],seq[j],index) += q1;
            model.mom2_ene2(seq[i],seq[j],index) += q2;
/*
            if(seq[i]!=seq[j]){
              model.mom2_ene1(seq[j],seq[i],index) += q1;//symmetrize
              model.mom2_ene2(seq[j],seq[i],index) += q2;
            }
*/
          }
        }
      }

      if((t+1)%(n/nseed/nblock)==0){
        block_cnt++;
      }
    }
  } 
 
  model.mom1 /= n;
  model.mom2 /= n; 
  model.mom1_ene1 /= n;
  model.mom1_ene2 /= n;
  model.mom2_ene1 /= n;
  model.mom2_ene2 /= n;
  model.avg_ene /= n;

  for(int i=0; i<nblock; i++){
    mom1_temp[i] /= (n/nblock);
    mom2_temp[i] /= (n/nblock);

    model.mom1_err += arma::accu(arma::square(mom1_temp[i]-model.mom1))/model.mom1.n_elem;
    model.mom2_err += arma::accu(arma::square(mom2_temp[i]-model.mom2))/model.mom2.n_elem;
  }
  model.mom1_err /= (nblock-1);
  model.mom2_err /= (nblock-1);
  model.mom1_err = sqrt(model.mom1_err);
  model.mom2_err = sqrt(model.mom2_err); 

  
}

void do_sweep(model &model, std::vector<int> &seq){

  for(int i=0; i<model.N; i++){
    unsigned long int s = gsl_rng_uniform_int(rg, (unsigned long int)model.N);
    unsigned long int r = gsl_rng_uniform_int(rg, (unsigned long int)model.q);
    while((int)r==seq[s]) r = gsl_rng_uniform_int(rg, (unsigned long int)model.q);
    double prob = std::min(1.0, get_boltzmann(model, seq, (int)s, (int)r));
    double xsi = gsl_rng_uniform(rg);
    if(prob>xsi) seq[s] = (int)r; 
  }
}

double get_boltzmann(model &model, std::vector<int> &seq, int m, int r){

  double E1 = model.get_energy(seq);  
  int old = seq[m];
  seq[m] = r;
  double E2 = model.get_energy(seq); 
  seq[m] = old;
  
  return exp(-(E2-E1));
}

