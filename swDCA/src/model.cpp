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

double model::get_delta_energy(std::vector<int> &seq, int m, int r){

  double dE = 0.;

  dE -= h(m,r) - h(m,seq[m]);

  for (int i=0; i<m; i++) {
      int index = (N-1)*i-i*(i+1)/2+m-1;
      if(symmetrize_on) {
	      dE -= (J(seq[i], r, index) + J(r, seq[i], index))/2;
	      dE += (J(seq[i], seq[m], index) + J(seq[m], seq[i], index))/2;
      } else { 
	      dE -= J(seq[i], r, index);
	      dE += J(seq[i], seq[m], index);
      }
  }

  for(int j=(m+1); j<N; j++){
      int index = (N-1)*m-m*(m+1)/2+j-1;
      if(symmetrize_on) {
	      dE -= (J(r, seq[j], index) + J(seq[j], r, index))/2;
	      dE += (J(seq[m], seq[j], index) + J(seq[j], seq[m], index))/2;
      } else {
	      dE -= J(r, seq[j], index);
	      dE += J(seq[m], seq[j], index);
      }
  }

  return dE;
}


/*
 * this version is not any faster; presumably the compiler takes 
 * care of this already
 *
double model::get_delta_energy(std::vector<int> &seq, int m, int r){

  double dE = 0.;

  dE -= h(m,r) - h(m,seq[m]);

  int seqm = seq[m];
  int save1 = m-1;
  int save2 = N-1;
  for (int i=0; i<m; i++) {
      int seqi = seq[i];
      int index = save2*i-i*(i+1)/2+save1;
      if(symmetrize_on) {
	      dE -= (J(seqi, r, index) + J(r, seqi, index))/2;
	      dE += (J(seqi, seqm, index) + J(seqm, seqi, index))/2;
      } else { 
	      dE -= J(seqi, r, index);
	      dE += J(seqi, seqm, index);
      }
  }

  int save = (N-1)*m-m*(m+1)/2-1;
  for(int j=(m+1); j<N; j++){
      int seqj = seq[j];
      //int index = (N-1)*m-m*(m+1)/2+j-1;
      int index = save+j;
      if(symmetrize_on) {
	      dE -= (J(r, seqj, index) + J(seqj, r, index))/2.;
	      dE += (J(seqm, seqj, index) + J(seqj, seqm, index))/2.;
      } else {
	      dE -= J(r, seqj, index);
	      dE += J(seqm, seqj, index);
      }
  }

  return dE;
}
*/

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

  std::cout << "Converting to zero sum gauge..." << std::endl;
  
  for(int i=0; i<N; i++){
    double sum=0;
    for(int a=0; a<q; a++) sum += h(i,a);
    for(int a=0; a<q; a++) h(i,a) -= sum/q;
    for(int a=0; a<q; a++){
      for(int q1=0; q1<q; q1++){
        for(int j=0; j<N; j++){
          if(j!=i){
            if(j>i){
              int index = (N-1)*i-i*(i+1)/2+j-1;
              h(i,a) += J(a,q1,index)/q;
              for(int q2=0; q2<q; q2++){
                h(i,a) -= J(q1,q2,index)/(q*q);
              }
            }
            else{
              int index = (N-1)*j-j*(j+1)/2+i-1;
              h(i,a) += J(q1,a,index)/q;
              for(int q2=0; q2<q; q2++){
                h(i,a) -= J(q1,q2,index)/(q*q);
              }
            }
          }
        }
      }
    }
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
