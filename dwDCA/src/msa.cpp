#include "msa.h"
#include "model.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <numeric>
#include <string>
#include <iostream>
#include <fstream>
#include <armadillo>

int letter_to_number(char a, int q){

  int n = -1;
  if(q==21){
    if(a=='-') n=0;
    else if(a=='A') n=1;
    else if(a=='C') n=2;
    else if(a=='D') n=3;
    else if(a=='E') n=4;
    else if(a=='F') n=5;
    else if(a=='G') n=6;
    else if(a=='H') n=7;
    else if(a=='I') n=8;
    else if(a=='K') n=9;
    else if(a=='L') n=10;
    else if(a=='M') n=11;
    else if(a=='N') n=12;
    else if(a=='P') n=13;
    else if(a=='Q') n=14;
    else if(a=='R') n=15;
    else if(a=='S') n=16;
    else if(a=='T') n=17;
    else if(a=='V') n=18;
    else if(a=='W') n=19;
    else if(a=='Y') n=20;
    else{
      n=0;
    }
  }
  else if(q==2){
    if(a=='H') n=0;
    else n=1;
  }
  else{
    printf("Error: unsupported number of Potts states. Exiting.\n");
    exit(-1);
  }
  return n;
}

std::string number_to_letter(int a, int q){

  std::string s;
  if(q==21){
    if(a==0) s='-';
    else if(a==1) s='A';
    else if(a==2) s='C';
    else if(a==3) s='D';
    else if(a==4) s='E';
    else if(a==5) s='F';
    else if(a==6) s='G';
    else if(a==7) s='H';
    else if(a==8) s='I';
    else if(a==9) s='K';
    else if(a==10) s='L';
    else if(a==11) s='M';
    else if(a==12) s='N';
    else if(a==13) s='P';
    else if(a==14) s='Q';
    else if(a==15) s='R';
    else if(a==16) s='S';
    else if(a==17) s='T';
    else if(a==18) s='V';
    else if(a==19) s='W';
    else if(a==20) s='Y';
    else{
      s='?';
    }
  }
  else if(q==2){
    if(a==0) s='H';
    else s='P';
  }
  else{
    printf("Error: unsupported number of Potts states. Exiting.\n");
    exit(-1);
  }
  return s;
}

std::vector<int> read_msa(std::string msafile, int *N, int*q){

  int nseq = 0;

  //Get number of Potts states
  int index = msafile.find(".");
  std::string extension = msafile.substr(index+1);
  if(extension=="fas") *q=21;
  else if(extension=="hp") *q=2;
  else{
    std::cout << "Error: file extension not recognized. Exiting." << std::endl;
    exit(-1);
  }

  //Count number of lines, divide by 2 to get # of seqs
  std::string line;
  std::ifstream myfile(msafile);
  while(std::getline(myfile,line)) nseq++;
  nseq = nseq/2;
  std::cout << "No. of sequences: " << nseq << std::endl;

  //Get sequence length
  std::ifstream stream(msafile);
  std::string dummyLine;
  std::string firstSeq;
  std::getline(stream,dummyLine);
  std::getline(stream,firstSeq);
  *N = firstSeq.length();
  std::cout << "Sequence length: " << *N << std::endl;

  std::vector<int> msa_seqs(nseq*(*N));

  int linecount=0;
  std::ifstream file(msafile);
  std::string str;
  while(std::getline(file, str)){
    if(linecount%2==1){
      for(int i=0; i<(*N); i++){
        msa_seqs[(linecount-1)/2*(*N)+i] = letter_to_number(str.at(i),*q);
      }
    }
    linecount++;
  }

  return msa_seqs;
}

std::vector<double> get_weights(std::vector<int> &msa_seqs, int nseq, double delta){

  std::cout << "Constructing sequence weights..." << std::endl;
  std::vector<double> weights(nseq, 0.0);
  int N = msa_seqs.size()/nseq;
  double Meff=0;
  for(int m1=0; m1<nseq; m1++){
    weights[m1] += 1.0;
    if(delta>0.0){
    std::vector<int> seq1(msa_seqs.begin()+m1*N, msa_seqs.begin()+(m1+1)*N);
    for(int m2=0; m2<nseq; m2++){
      if(m1!=m2){
        std::vector<int> seq2(msa_seqs.begin()+m2*N, msa_seqs.begin()+(m2+1)*N);
        double dist = get_hamming_dist(seq1, seq2);
        if(dist<delta*N) weights[m1] += 1.0;
      }
    }
    }
    weights[m1] = 1.0/weights[m1];
  }
  
  for(int m1=0; m1<nseq; m1++) Meff += weights[m1];
  std::cout << "Computed weights. Meff=" << Meff << std::endl;
  return weights;
}

void fill_freq(std::vector<int> &msa_seqs, std::vector<double> &weights,
		arma::mat &msa_freq, arma::cube &msa_corr,int nseq, bool symon){

  //double Meff = std::accumulate(weights.begin(), weights.end(), 0);
  std::cout << nseq << std::endl;
  double Meff=0;
  for(int i=0; i<nseq; i++) Meff+=weights[i];
  std::cout << "Checking Meff: " << Meff << std::endl;
  int N = msa_seqs.size()/nseq;
  
  //Set freqs and corrs to zero
  msa_freq.zeros();
  msa_corr.zeros();

  //Compute frequencies
  for(int m=0; m<nseq; m++){
    for(int i=0; i<N; i++){
      msa_freq(i,msa_seqs[m*N+i]) += weights[m];
    }
    for(int i=0; i<(N-1); i++){
      for(int j=(i+1); j<N; j++){
        int index = (N-1)*i-i*(i+1)/2+j-1;
        int ind1 = msa_seqs[m*N+i];
        int ind2 = msa_seqs[m*N+j];
        msa_corr(ind1,ind2,index) += weights[m];
        if(symon){
          if(ind1!=ind2){
            msa_corr(ind2,ind1,index) += 0.5*weights[m];
            msa_corr(ind1,ind2,index) -= 0.5*weights[m];
          }
        }
      }
    } 
  }

  //Normalize
  msa_freq /= Meff;
  msa_corr /= Meff;
}

void fill_ene_weighted_freq(std::vector<int> &msa_seqs, std::vector<double> &weights,
		arma::mat &freq1, arma::mat &freq2, arma::cube &corr1, arma::cube &corr2, model &mymodel, int nseq){

  double Meff = std::accumulate(weights.begin(), weights.end(), 0);
  int N = mymodel.N;
  double Tmix = mymodel.Tmix;
  
  //Set freqs and corrs to zero
  freq1.zeros();
  freq2.zeros();
  corr1.zeros();
  corr2.zeros();

  //Compute energy weights q1, q2
  std::vector<double> q1(nseq, 0.0);
  std::vector<double> q2(nseq, 0.0);
  for(int m=0; m<nseq; m++){
    std::vector<int> seq_curr(msa_seqs.begin()+m*N, msa_seqs.begin()+(m+1)*N);
    double E1 = mymodel.get_energy_single(seq_curr, 1);
    double E2 = mymodel.get_energy_single(seq_curr, 2);
    double exp1 = exp(-E1/Tmix);
    double exp2 = exp(-E2/Tmix);
    q1[m] = exp1/(exp1+exp2);
    q2[m] = exp2/(exp1+exp2);
  } 

  //Compute frequencies
  for(int m=0; m<nseq; m++){
    for(int i=0; i<N; i++){
      freq1(i,msa_seqs[m*N+i]) += weights[m]*q1[m];
      freq2(i,msa_seqs[m*N+i]) += weights[m]*q2[m];
    }
    for(int i=0; i<(N-1); i++){
      for(int j=(i+1); j<N; j++){
        int index = (N-1)*i-i*(i+1)/2+j-1;
        int ind1 = msa_seqs[m*N+i];
        int ind2 = msa_seqs[m*N+j];
        corr1(ind1,ind2,index) += weights[m]*q1[m];
        corr2(ind1,ind2,index) += weights[m]*q2[m];

/*
        if(ind1!=ind2){
	  corr1(ind2,ind1,index) += 0.5*weights[m]*q1[m];
	  corr2(ind2,ind1,index) += 0.5*weights[m]*q2[m];
	  corr1(ind1,ind2,index) -= 0.5*weights[m]*q1[m];
	  corr2(ind1,ind2,index) -= 0.5*weights[m]*q2[m];
        }
*/
      }
    } 
  }

  //Normalize
  freq1 /= Meff;
  freq2 /= Meff;
  corr1 /= Meff;
  corr2 /= Meff;
}
double get_hamming_dist(std::vector<int> v1, std::vector<int> v2){
  
  if(v1.size()!=v2.size()){
    std::cout << "ERROR: vectors need the same size!" << std::endl;
    exit(-1);
  }
  double dist=0;
  for(int i=0; i<v1.size(); i++){
    if(v1[i]!=v2[i]) dist+=1;
  }
  return dist;
}
