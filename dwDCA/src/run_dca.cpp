/* run_dca.cpp Created by Layne Frechette on Feb. 17, 2021
 *
 * Run DCA to fit MSA to single- or double-well Potts model.
 */

//Headers
#include "mc.h"
#include "io.h"
#include "rand.h"
#include "model.h"
#include "msa.h"
#include "run_dca.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
using std::istringstream;

//Global variables

std::string msa_name;
std::string out_name;
std::string input_name;
std::string conf_name;
std::string mc_init;
std::string param_init;

std::string scratch_dir;
std::string output_dir;
std::string folder_name;

gsl_rng *rg;
gsl_rng **rg_replica;

int max_iter;
int mc_steps;
int nrep; 
double eps0_h;
double eps0_J;
double eps_inc;
double eps_dec;
double eps_min_h;
double eps_max_h;
double eps_min_J;
double eps_max_J_N;
double lambda;
double gamma_mom;
double cutoff_freq;
bool adaptive_stepsize_on;
bool symmetrize_on;
bool cd_on;
bool verbose;
int nwell;
double Tmix=1.0;
double delta=0.2; //reweighting threshold
long unsigned int myseed;

/*** Main ***/
int main(int argc, char *argv[]){

  //Read command line arguments
  if(argc == 2){
    //Try reading in parameters from conf file
    conf_name = argv[1];
    if(fs::exists(conf_name)){
      read_inputs(conf_name);
    }
    else{ 
      printf("ERROR: conf file does not exist! Exiting.\n");
      exit(-1);
    }
  }
  else{
    printf("ERROR: conf file required (1 argument).\n");
    exit(-1);
  }

  //Set up random number generator
  init_rng(myseed);
  
  //Memory allocation
  int N;
  int q;
  int nseq;

  std::cout << msa_name << std::endl;
  std::vector<int> msa_seqs = read_msa(msa_name, &N, &q);
  nseq = msa_seqs.size()/N;
  std::cout << nseq << " sequences" << std::endl;
  std::cout << "Tmix: " << Tmix << std::endl;
  model model(N, q, nwell, lambda, symmetrize_on, Tmix, mc_init);
  std::cout << model.mc_init << std::endl;
  arma::mat msa_freq(model.N,model.q,arma::fill::ones);
  arma::cube msa_corr(model.q, model.q, model.N*(model.N-1)/2,arma::fill::zeros);
  //msa_freq = msa_freq/21.0;

  //Define data output directory
  output_dir = "data/";

  //scratch_dir = "/data/frechettelb/protein_evolution/dwDCA/data/";
  scratch_dir = "data/";

  output_dir += folder_name + "/";
  scratch_dir += folder_name + "/";

  fs::create_directories(output_dir);
  fs::create_directories(scratch_dir);

  //Get frequencies and pair correlations from MSA
  std::vector<double> weights = get_weights(msa_seqs, nseq, delta);
  fill_freq(msa_seqs, weights, msa_freq, msa_corr, nseq, symmetrize_on);
  std::string msa_stats_2 = scratch_dir + "/stat_align_2p.txt";
  msa_freq.save(scratch_dir + "stat_align_1p.txt", arma::arma_ascii);
  msa_corr.save(scratch_dir + "stat_align_2p.txt", arma::arma_ascii);

  if(N==2){
    std::cout << "msa frequencies:" << std::endl;
    std::cout << msa_freq << std::endl;
    std::cout << "msa correlations" << std::endl;
    std::cout << msa_corr << std::endl;
  }

  //Initialize parameters
  if(input_name=="none"){
    if(model.nwell==1){
      //Initialize h's to be consistent with single-site frequencies
      for(int i=0; i<model.N; i++){
        for(int a=0; a<model.q; a++){
          model.h1(i,a) = log(msa_freq(i,a)+0.0001); //Small regularization to avoid divergence
        }
      }
    }
    if(model.nwell==2){
      //Make initial parameters different for the two wells
      //by assigning half the fields to one well and half to the other.
      if(param_init=="half"){
        std::cout << "HALF" << std::endl;
        for(int i=0; i<model.N/2; i++){
          for(int a=0; a<model.q; a++){
            model.h1(i,a) = log(msa_freq(i,a)+0.0001);
          }
        }
        for(int i=model.N/2; i<model.N; i++){
          for(int a=0; a<model.q; a++){
            model.h2(i,a) = log(msa_freq(i,a)+0.0001);
          }
        }
      }
      if(param_init=="all_in_one"){
        std::cout << "all_in_one" << std::endl;
        for(int i=0; i<model.N; i++){
          for(int a=0; a<model.q; a++){
            model.h1(i,a) = log(msa_freq(i,a)+0.0001);
          }
        }
      }
    }
  }
  else{
    read_params(model,input_name);
  }

  //Ensure parameters are in zero-sum gauge
  //model.convert_to_zero_sum();

  //Run DCA
  fit(model,msa_freq,msa_corr,msa_seqs,weights,nrep);

  if(model.N==2){
    std::cout << "h1:" << std::endl;
    std::cout << model.h1 << std::endl;
    std::cout << "h2:" << std::endl;
    std::cout << model.h2 << std::endl;
    std::cout << "J1:" << std::endl;
    std::cout << model.J1 << std::endl;
    std::cout << "J2:" << std::endl;
    std::cout << model.J2 << std::endl;
  }

  //Print parameters to file
  print_params(model,output_dir+out_name);
  
  return 1;
}

void fit(model &mymodel, arma::mat &msa_freq, arma::cube &msa_corr, std::vector<int> msa_seqs, std::vector<double> weights, int nr){

  bool converged = false;
  int niter=0;

  //Initialize adaptive learning rates
  arma::mat alpha_h1(mymodel.N,mymodel.q,arma::fill::ones);
  arma::cube alpha_J1(mymodel.q, mymodel.q, mymodel.N*(mymodel.N-1)/2,arma::fill::ones);
  alpha_h1 = eps0_h*alpha_h1;
  alpha_J1 = eps0_J*alpha_J1;

  arma::mat alpha_h2(mymodel.N,mymodel.q,arma::fill::ones);
  arma::cube alpha_J2(mymodel.q, mymodel.q, mymodel.N*(mymodel.N-1)/2,arma::fill::ones);
  alpha_h2 = eps0_h*alpha_h2;
  alpha_J2 = eps0_J*alpha_J2;

  //Derivatives
  arma::mat dh1(mymodel.N,mymodel.q,arma::fill::zeros);
  arma::cube dJ1(mymodel.q, mymodel.q, mymodel.N*(mymodel.N-1)/2,arma::fill::zeros);

  arma::mat dh2(mymodel.N,mymodel.q,arma::fill::zeros);
  arma::cube dJ2(mymodel.q, mymodel.q, mymodel.N*(mymodel.N-1)/2,arma::fill::zeros);

  //Stored derivatives of previous
  arma::mat dh1_prev(mymodel.N,mymodel.q,arma::fill::ones);
  arma::cube dJ1_prev(mymodel.q, mymodel.q, mymodel.N*(mymodel.N-1)/2,arma::fill::ones);

  arma::mat dh2_prev(mymodel.N,mymodel.q,arma::fill::ones);
  arma::cube dJ2_prev(mymodel.q, mymodel.q, mymodel.N*(mymodel.N-1)/2,arma::fill::ones);

  //Changes in parameters (for momentum)
  arma::mat change_h1(mymodel.N,mymodel.q,arma::fill::zeros);
  arma::cube change_J1(mymodel.q,mymodel.q,mymodel.N*(mymodel.N-1)/2,arma::fill::zeros);

  arma::mat change_h2(mymodel.N,mymodel.q,arma::fill::zeros);
  arma::cube change_J2(mymodel.q,mymodel.q,mymodel.N*(mymodel.N-1)/2,arma::fill::zeros);

  //Average energy
  arma::vec avg_energies(max_iter, arma::fill::zeros);

  //Cross entropy
  arma::vec entropy(max_iter, arma::fill::zeros);

  //Fit
  std::cout << "Fitting model." << std::endl;
  while(!converged){

    std::cout << std::endl << "iteration " << niter << std::endl;
    run_mc_traj(mymodel,mc_steps,nr);
    avg_energies[niter] = mymodel.avg_ene;
    std::cout << "Avg. energy: " << avg_energies[niter] << std::endl;
    print_seqs(mymodel, scratch_dir + "seqs_" + std::to_string(niter) + ".txt");

    //Compute RMS difference between MSA and MC
    double diff1 = sqrt(arma::accu(arma::square(dh1_prev))/mymodel.mom1.n_elem);
    double diff2 = sqrt(arma::accu(arma::square(dJ1_prev))/mymodel.mom2.n_elem);
    double max1 = arma::abs(dh1_prev).max();
    double max2 = arma::abs(dJ1_prev).max();
    double diff1_2 = sqrt(arma::accu(arma::square(dh2_prev))/mymodel.mom1.n_elem);
    double diff2_2 = sqrt(arma::accu(arma::square(dJ2_prev))/mymodel.mom2.n_elem);
    double max1_2 = arma::abs(dh2_prev).max();
    double max2_2 = arma::abs(dJ2_prev).max();
    std::cout << "RMS of dh1, dJ1: " << diff1 << " " << diff2 << std::endl;
    std::cout << "Max of dh1, dJ1: " << max1 << " " << max2 << std::endl;
    if(mymodel.nwell==2){
      std::cout << "RMS of dh2, dJ2: " << diff1_2 << " " << diff2_2 << std::endl;
      std::cout << "Max of dh2, dJ2: " << max1_2 << " " << max2_2 << std::endl;
    }
    if(mymodel.N==2){
      std::cout << "h1:" << std::endl;
      std::cout << mymodel.h1 << std::endl;
      std::cout << "h2:" << std::endl;
      std::cout << mymodel.h2 << std::endl;
      std::cout << "J1:" << std::endl;
      std::cout << mymodel.J1 << std::endl;
      std::cout << "J2:" << std::endl;
      std::cout << mymodel.J2 << std::endl;
      std::cout << "model 1p freq:" << std::endl;
      std::cout << mymodel.mom1 << std::endl;
      std::cout << "model 2p freq:" << std::endl;
      std::cout << mymodel.mom2 << std::endl;
    }
    if(diff1<mymodel.mom1_err && diff2<mymodel.mom2_err){
      std::cout << "Std. errors: " << mymodel.mom1_err << " " << mymodel.mom2_err << std::endl;
      std::cout << "Model difference is within sampling error. Increasing sampling..." << std::endl;
      mc_steps *= 2;
      std::cout << "MC steps per iteration: " << mc_steps << std::endl;
    }
/*
    if(mymodel.q==2){
//      double Z = mymodel.get_Z();
//      entropy[niter] = log(Z)-arma::accu(mymodel.mom1%mymodel.h1)-arma::accu(mymodel.mom2%mymodel.J1); 
      std::cout << "Cross-entropy: " << entropy[niter] << std::endl;
    }
*/
    if (verbose) {
    //Dump current statistics and parameters to file
    mymodel.mom1.save(scratch_dir + "stat_MC_1p_" + std::to_string(niter) + ".txt", arma::arma_ascii);
    mymodel.mom2.save(scratch_dir + "stat_MC_2p_" + std::to_string(niter) + ".txt",arma::arma_ascii);
    mymodel.h1.save(scratch_dir + "h1_" + std::to_string(niter) + ".txt", arma::arma_ascii);
    mymodel.J1.save(scratch_dir + "J1_" + std::to_string(niter) + ".txt", arma::arma_ascii);
    if(mymodel.nwell==2){
      mymodel.h2.save(scratch_dir + "h2_" + std::to_string(niter) + ".txt", arma::arma_ascii);
      mymodel.J2.save(scratch_dir + "J2_" + std::to_string(niter) + ".txt",arma::arma_ascii);
    }
    }
    
    /*** Compute derivatives ***/
    if(mymodel.nwell==1){ 
      dh1 = (msa_freq-mymodel.mom1 - 2*mymodel.lambda*mymodel.h1);
      dJ1 = (msa_corr-mymodel.mom2 - 2*mymodel.lambda*mymodel.J1);

      change_h1 = gamma_mom*change_h1 + alpha_h1%dh1;
      change_J1 = gamma_mom*change_J1 + alpha_J1%dJ1;
    }
    else{
      //Get energy-weighted frequencies
      arma::mat msa_freq1(mymodel.N, mymodel.q, arma::fill::zeros);
      arma::mat msa_freq2(mymodel.N, mymodel.q, arma::fill::zeros);
      arma::cube msa_corr1(mymodel.q, mymodel.q, mymodel.N*(mymodel.N-1)/2,arma::fill::zeros);
      arma::cube msa_corr2(mymodel.q, mymodel.q, mymodel.N*(mymodel.N-1)/2,arma::fill::zeros);

      fill_ene_weighted_freq(msa_seqs, weights, msa_freq1, msa_freq2, msa_corr1, msa_corr2, mymodel, weights.size());

      if(mymodel.N==2){
        std::cout << "msa freqs:" << std::endl;
        std::cout << msa_freq1 << std::endl;
        std::cout << msa_freq2 << std::endl;
        std::cout << "msa corrs:" << std::endl;
        std::cout << msa_corr1 << std::endl;
        std::cout << msa_corr2 << std::endl;
      }

      dh1 = (msa_freq1-mymodel.mom1_ene1 - 2*mymodel.lambda*mymodel.h1);
      dJ1 = (msa_corr1-mymodel.mom2_ene1 - 2*mymodel.lambda*mymodel.J1);
      dh2 = (msa_freq2-mymodel.mom1_ene2 - 2*mymodel.lambda*mymodel.h2);
      dJ2 = (msa_corr2-mymodel.mom2_ene2 - 2*mymodel.lambda*mymodel.J2);

      if(mymodel.N==2){
        std::cout << "dJ1:" << std::endl;
        std::cout << dJ1 << std::endl;
        std::cout << "dJ2:" << std::endl;
        std::cout << dJ2 << std::endl;
      }

      change_h1 = gamma_mom*change_h1 + alpha_h1%dh1;
      change_J1 = gamma_mom*change_J1 + alpha_J1%dJ1;
      change_h2 = gamma_mom*change_h2 + alpha_h2%dh2;
      change_J2 = gamma_mom*change_J2 + alpha_J2%dJ2;
    }

    //Update parameters
    mymodel.h1 += change_h1;
    mymodel.J1 += change_J1;

    if(mymodel.nwell==2){
      mymodel.h2 += change_h2;
      mymodel.J2 += change_J2;
    }

    if (verbose) {
   //Dump derivatives to file
    dh1.save(scratch_dir + "dh1_" + std::to_string(niter) + ".txt", arma::arma_ascii);
    dJ1.save(scratch_dir + "dJ1_" + std::to_string(niter) + ".txt", arma::arma_ascii);

    if(mymodel.nwell==2){
      dh2.save(scratch_dir + "dh2_" + std::to_string(niter) + ".txt", arma::arma_ascii);
      dJ2.save(scratch_dir + "dJ2_" + std::to_string(niter) + ".txt", arma::arma_ascii);
    }
    }

    //Update learning rates
    if(adaptive_stepsize_on){
      for(int i=0; i<mymodel.N; i++){
        for(int a=0; a<mymodel.q; a++){
          double prod1 = dh1_prev(i,a)*dh1(i,a);
          if(prod1>0 && alpha_h1(i,a)<eps_max_h) alpha_h1(i,a) = eps_inc*alpha_h1(i,a);
          else if(prod1<0 && alpha_h1(i,a)>eps_min_h) alpha_h1(i,a) = eps_dec*alpha_h1(i,a);
          else {};
          if(alpha_h1(i,a)>eps_max_h) alpha_h1(i,a)=eps_max_h;
          if(alpha_h1(i,a)<eps_min_h) alpha_h1(i,a)=eps_min_h;
          if(mymodel.nwell==2){
            double prod2 = dh2_prev(i,a)*dh2(i,a);
            if(prod2>0 && alpha_h2(i,a)<eps_max_h) alpha_h2(i,a) = eps_inc*alpha_h2(i,a);
            else if(prod2<0 && alpha_h2(i,a)>eps_min_h) alpha_h2(i,a) = eps_dec*alpha_h2(i,a);
            else {};
            if(alpha_h2(i,a)>eps_max_h) alpha_h2(i,a)=eps_max_h;
            if(alpha_h2(i,a)<eps_min_h) alpha_h2(i,a)=eps_min_h;
          }
        }
      }
      for(int i=0; i<mymodel.N-1; i++){
        for(int j=i+1; j<mymodel.N; j++){ 
          int index = (mymodel.N-1)*i-i*(i+1)/2+j-1;
          for(int a=0; a<mymodel.q; a++){
            for(int b=0; b<mymodel.q; b++){
              double prod1 = dJ1_prev(a,b,index)*dJ1(a,b,index);
              if(prod1>0 && alpha_J1(a,b,index)<eps_max_J_N/mymodel.N) alpha_J1(a,b,index) = eps_inc*alpha_J1(a,b,index);
              else if(prod1<0 && alpha_J1(a,b,index)>eps_min_J) alpha_J1(a,b,index) = eps_dec*alpha_J1(a,b,index);
              else {};
              if(alpha_J1(a,b,index)>eps_max_J_N/mymodel.N) alpha_J1(a,b,index)=eps_max_J_N/mymodel.N;
              if(alpha_J1(a,b,index)<eps_min_J) alpha_J1(a,b,index)=eps_min_J;
              if(mymodel.nwell==2){
                double prod2 = dJ2_prev(a,b,index)*dJ2(a,b,index);
                if(prod2>0 && alpha_J2(a,b,index)<eps_max_J_N/mymodel.N) alpha_J2(a,b,index) = eps_inc*alpha_J2(a,b,index);
                else if(prod2<0 && alpha_J2(a,b,index)>eps_min_J) alpha_J2(a,b,index) = eps_dec*alpha_J2(a,b,index);
                else {};
                if(alpha_J2(a,b,index)>eps_max_J_N/mymodel.N) alpha_J2(a,b,index)=eps_max_J_N/mymodel.N;
                if(alpha_J2(a,b,index)<eps_min_J) alpha_J2(a,b,index)=eps_min_J;
              }
            }
          }
        }
      }
    }

    //Update previous step gradient
    dh1_prev = dh1;
    dJ1_prev = dJ1;
    if(mymodel.nwell==2){
      dh2_prev = dh2;
      dJ2_prev = dJ2;
    }
    
    niter++;
    if(mymodel.nwell==1){
      if(niter>max_iter || (max1<cutoff_freq && max2<cutoff_freq)){
        std::cout << "Converged!" << std::endl;
        std::cout << "RMS of dh, dJ: " << diff1 << " " << diff2 << std::endl;
        std::cout << "Max of |dh|, |dJ|: " << max1 << " " << max2 << std::endl;
        converged=true;
      }
    }
    if(mymodel.nwell==2){
      if(niter>max_iter || (max1<cutoff_freq && max2<cutoff_freq && max1_2<cutoff_freq && max2_2<cutoff_freq)){
        std::cout << "Converged!" << std::endl;
        std::cout << "RMS of dh1, dJ1, dh2, dJ2: " << diff1 << " " << diff2 << " " << diff1_2 << " " << diff2_2 << std::endl;
        std::cout << "Max of |dh1|, |dJ1|, |dh2|, |dJ2|: " << max1 << " " << max2 << " " << max1_2 << " " << max2_2 << std::endl;
        converged=true;
      }
    }
  }

  //Dump energies to file
  avg_energies.save(output_dir + "avg_ene.txt", arma::arma_ascii);
  entropy.save(output_dir + "entropy.txt", arma::arma_ascii);
  //mymodel.convert_to_zero_sum();

}

void read_inputs(std::string name){

  nrep=1;
  verbose=false;
  //Adapted from https://stackoverflow.com/questions/7868936/read-file-line-by-line-using-ifstream-in-c
  std::ifstream file(name);
  if(file.is_open()){
    std::string line;
    std::cout << "INPUT PARAMETERS:" << std::endl;
    while(std::getline(file,line)){
      int index = line.find(":");
      std::string value = line.substr(index+1);
      std::string key = line.substr(0,index);
      std::cout << key << " " << value << std::endl;
      if(key.compare("msa_name")==0) msa_name = value;
      if(key.compare("out_name")==0) out_name = value;
      if(key.compare("input_name")==0) input_name = value;
      if(key.compare("folder_name")==0) folder_name = value;
      if(key.compare("lambda")==0) lambda = std::stof(value);
      if(key.compare("myseed")==0) myseed = (long unsigned)std::stoi(value);
      if(key.compare("max_iter")==0) max_iter = std::stoi(value);
      if(key.compare("mc_steps")==0) mc_steps = std::stoi(value);
      if(key.compare("eps0_h")==0) eps0_h = std::stof(value);
      if(key.compare("eps0_J")==0) eps0_J = std::stof(value);
      if(key.compare("adaptive_stepsize_on")==0) adaptive_stepsize_on = std::stoi(value);
      if(key.compare("cd_on")==0) cd_on = std::stoi(value);
      if(key.compare("eps_inc")==0) eps_inc = std::stof(value);
      if(key.compare("eps_dec")==0) eps_dec = std::stof(value);
      if(key.compare("eps_min_h")==0) eps_min_h = std::stof(value);
      if(key.compare("eps_max_h")==0) eps_max_h = std::stof(value);
      if(key.compare("eps_min_J")==0) eps_min_J = std::stof(value);
      if(key.compare("eps_max_J_N")==0) eps_max_J_N = std::stof(value);
      if(key.compare("gamma")==0) gamma_mom = std::stof(value);
      if(key.compare("cd_on")==0) cd_on = std::stoi(value);
      if(key.compare("cutoff_freq")==0) cutoff_freq = std::stof(value);
      if(key.compare("nwell")==0) nwell = std::stoi(value);
      if(key.compare("Tmix")==0) Tmix = std::stof(value);
      if(key.compare("delta")==0) delta = std::stof(value);
      if(key.compare("mc_init")==0) mc_init = value;
      if(key.compare("param_init")==0) param_init = value;
      if(key.compare("symmetrize_on")==0) symmetrize_on = std::stoi(value);
      if(key.compare("nrep")==0) nrep = std::stoi(value);
      if(key.compare("verbose")==0) verbose = std::stoi(value);
    }
    std::cout << std::endl;
    file.close();
  }  
}

void read_params(model &mymodel, std::string in_name){

  std::ifstream file(in_name);
  if(file.is_open()){
    std::string line;
    while(std::getline(file,line)){
      std::string::size_type pos;
      pos = line.find(' ',0);
      std::string data = line.substr(pos+1);
      std::string id = line.substr(0,pos);
      istringstream my_stream(data);
      if(id=="J1"){
        int i;
        int j;
        int a;
        int b;
        my_stream >> i;
        my_stream >> j;
        my_stream >> a;
        my_stream >> b;
        int index = (mymodel.N-1)*i-i*(i+1)/2+j-1; 
        my_stream >> mymodel.J1(a,b,index);
      }
      else if(id=="J2"){
        int i;
        int j;
        int a;
        int b;
        my_stream >> i;
        my_stream >> j;
        my_stream >> a;
        my_stream >> b;
        int index = (mymodel.N-1)*i-i*(i+1)/2+j-1; 
        my_stream >> mymodel.J2(a,b,index);
      }
      else if(id=="h1"){
        int i;
        int a;
        my_stream >> i;
        my_stream >> a;
        my_stream >> mymodel.h1(i,a);
      }
      else if(id=="h2"){
        int i;
        int a;
        my_stream >> i;
        my_stream >> a;
        my_stream >> mymodel.h2(i,a);
      }
      else{
        std::cout << "Error: unrecognized specifier. Exiting." << std::endl;
        exit(-1);
      }
    }
  }
}

void print_params(model &mymodel, std::string name){

  std::ofstream ofile;
  ofile.open(name);
  for(int i=0; i<mymodel.N-1; i++){
    for(int j=i+1; j<mymodel.N; j++){
      int index = (mymodel.N-1)*i-i*(i+1)/2+j-1;
      for(int a=0; a<mymodel.q; a++){
        for(int b=0; b<mymodel.q; b++){
          ofile << "J1 " << i << " " << j << " " << a << " " << b << " " << mymodel.J1(a,b,index) << std::endl;
        }
      }
    }
  }
  if(mymodel.nwell==2){
    for(int i=0; i<mymodel.N-1; i++){
      for(int j=i+1; j<mymodel.N; j++){
        int index = (mymodel.N-1)*i-i*(i+1)/2+j-1;
        for(int a=0; a<mymodel.q; a++){
          for(int b=0; b<mymodel.q; b++){
            ofile << "J2 " << i << " " << j << " " << a << " " << b << " " << mymodel.J2(a,b,index) << std::endl;
          }
        }
      }
    }
  }
  for(int i=0; i<mymodel.N; i++){
    for(int a=0; a<mymodel.q; a++){
      ofile << "h1 " << i << " " << a << " " << mymodel.h1(i,a) << std::endl;
    }
  }
  if(mymodel.nwell==2){
    for(int i=0; i<mymodel.N; i++){
      for(int a=0; a<mymodel.q; a++){
        ofile << "h2 " << i << " " << a << " " << mymodel.h2(i,a) << std::endl;
      }
    }
  }
  ofile.close();
}

void print_seqs(model &mymodel, std::string name){

  std::ofstream ofile;
  ofile.open(name);
  for(int i=0; i<mymodel.num_seqs; i++){
    ofile << mymodel.seqs[i] << std::endl;
  }
  ofile.close();
}
