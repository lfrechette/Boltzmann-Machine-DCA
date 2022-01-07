/* model.cpp Created by Layne Frechette on Feb. 17, 2021
 *
 * Implements DCA model class.
 */

#include "model.h"
#include <stdlib.h>
#include <vector>


model::model(){
  lambda = 0.01;
}

model::model(int N1, int q1, double lam, bool symon, std::string mc_init1){
  N = N1;
  q = q1;
  lambda = lam;
  avg_ene = 0;
  num_seqs = 100;
  seqs.resize(num_seqs);
  h.zeros(N,q);
  J.zeros(q,q,N*(N-1)/2);
  mom1.zeros(N,q);
  mom2.zeros(q,q,N*(N-1)/2);
  mom1_err=0;
  mom2_err=0;
  mc_init = mc_init1;
  symmetrize_on=symon;
}

model::~model(){
}

double model::get_energy(std::vector<int> &seq){

  double energy = 0;

  for(int i=0; i<N; i++){
    energy -= h(i,seq[i]);
  }
  for(int i=0; i<(N-1); i++){
    for(int j=(i+1); j<N; j++){
      int index = (N-1)*i-i*(i+1)/2+j-1;
      if(symmetrize_on) energy -= (J(seq[i], seq[j], index) + J(seq[j], seq[i], index))/2;
      else energy -= J(seq[i], seq[j], index);
    }
  }

  return energy;
}


double model::get_Z(){

  //Uses code from https://www.geeksforgeeks.org/print-all-sequences-of-given-length/

  double Z=0;
  std::vector<int> seq(N);

  while(true){
    double ene = get_energy(seq); 
    Z += exp(-ene);
    for(int i=0; i<N; i++){
      printf("%d ", seq[i]);
    }
    printf("\n");
    int p=N-1;
    while(seq[p]==q-1) p--;
    if(p<0) break;
    seq[p]=seq[p]+1;
    for(int i=p+1; i<N; i++) seq[i]=0;
  }
  
  return Z;
}

void model::convert_to_zero_sum(){
  
  for(int i=0; i<N; i++){
    double sum=0;
    for(int a=0; a<q; a++) sum += h(i,a);
    for(int a=0; a<q; a++) h(i,a) -= sum/q;
  }
  for(int i=0; i<N-1; i++){
    for(int j=i+1; j<N; j++){
      int index = (N-1)*i-i*(i+1)/2+j-1;
      double dsum1=0;
      std::vector<double> sum1_a(q);
      std::vector<double> sum1_b(q);
      for(int a=0; a<q; a++){
        sum1_a[a]=0;
        sum1_b[a]=0;
      }
      for(int b=0; b<q; b++){
        for(int a=0; a<q; a++){
          sum1_b[b] += J(a,b,index)/q;
          dsum1 += J(a,b,index)/(q*q);
          sum1_a[a] += J(a,b,index)/q;
        }
      }
      for(int a=0; a<q; a++){
        for(int b=0; b<q; b++){
          J(a,b,index) = J(a,b,index) - sum1_a[a] - sum1_b[b] + dsum1;
        }
      }
    }
  }
}
