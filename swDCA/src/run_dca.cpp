/* run_dca.cpp Created by Layne Frechette on Feb. 17, 2021
 *
 * Run DCA to fit MSA to Potts model.
 */

//Headers
#include "mc.h"
#include "rand.h"
#include "model.h"
#include "msa.h"
#include "run_dca.h"
#include "omp.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <filesystem>
//#include <experimental/filesystem>
//namespace fs = std::experimental::filesystem; //use this for C++14
namespace fs = std::filesystem;
using std::istringstream;

//Global variables

std::string msa_name;
std::string out_name;
std::string input_name;
std::string conf_name;
std::string mc_init;

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
bool adaptive_sampling_on;
bool symmetrize_on; 
bool verbose;
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
  model model(N, q, lambda, symmetrize_on, mc_init);
  std::cout << model.mc_init << std::endl;
  arma::mat msa_freq(model.N,model.q,arma::fill::ones);
  arma::cube msa_corr(model.q, model.q, model.N*(model.N-1)/2,arma::fill::zeros);
  msa_freq = msa_freq/21.0;

  //Define data output directory
  output_dir = "data/";
  folder_name += "/lambda=" + std::to_string(lambda);
  if(adaptive_sampling_on) folder_name += "_as";
  if(adaptive_stepsize_on) folder_name += "_al";
  if(gamma_mom>0.0) folder_name += "_gamma=" + std::to_string(gamma_mom);

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

  //Initialize parameters
  if(input_name=="none") init_default(model, msa_freq);
  else read_params(model,input_name); 

  /***TEST FOR SW TOY MODEL---MOVE TO A TEST FILE***/

/*
  msa_freq(0,0)=0.5;
  msa_freq(0,1)=0.5;
  msa_freq(1,0)=0.5;
  msa_freq(0,1)=0.5;
  msa_corr(0,0,0)=0.440399;
  msa_corr(0,1,0)=0.059601;
  msa_corr(1,0,0)=0.059601;
  msa_corr(1,1,0)=0.440399;
*/
  model.convert_to_zero_sum();

  if(model.N==2){
    std::cout << "h:" << std::endl;
    std::cout << model.h << std::endl;
    std::cout << "J:" << std::endl;
    std::cout << model.J << std::endl;
  }

  //Run DCA
  fit(model,msa_freq,msa_corr,nrep);

  if(model.N==2){
    std::cout << "h:" << std::endl;
    std::cout << model.h << std::endl;
    std::cout << "J:" << std::endl;
    std::cout << model.J << std::endl;
  }

  //Ensure parameters are in zero-sum gauge
  model.convert_to_zero_sum();

  if(model.N==2){
    std::cout << "After zero-sum conversion:" << std::endl;
    std::cout << "h:" << std::endl;
    std::cout << model.h << std::endl;
    std::cout << "J:" << std::endl;
    std::cout << model.J << std::endl;
  }

  //Print parameters to file
  print_params(model,output_dir+out_name);
  
  return 1;
}

/***********
 * METHODS *
 ***********/

void init_default(model &mymodel, arma::mat &msa_freq){
  //Initialize h's to be consistent with single-site frequencies
  for(int i=0; i<mymodel.N; i++){
    for(int a=0; a<mymodel.q; a++){
      mymodel.h(i,a) = log(msa_freq(i,a)+0.0001); //Small regularization to avoid divergence
    }
  }
}

void fit(model &mymodel, arma::mat &msa_freq, arma::cube &msa_corr, int nr){

  bool converged = false;
  int niter=0;

  //Initialize adaptive learning rates
  arma::mat alpha_h(mymodel.N,mymodel.q,arma::fill::ones);
  arma::cube alpha_J(mymodel.q, mymodel.q, mymodel.N*(mymodel.N-1)/2,arma::fill::ones);
  alpha_h = eps0_h*alpha_h;
  alpha_J = eps0_J*alpha_J;

  //Derivatives
  arma::mat dh(mymodel.N,mymodel.q,arma::fill::zeros);
  arma::cube dJ(mymodel.q, mymodel.q, mymodel.N*(mymodel.N-1)/2,arma::fill::zeros);

  //Stored derivatives of previous
  arma::mat dh_prev(mymodel.N,mymodel.q,arma::fill::ones);
  arma::cube dJ_prev(mymodel.q, mymodel.q, mymodel.N*(mymodel.N-1)/2,arma::fill::ones);

  //Changes in parameters (for momentum)
  arma::mat change_h(mymodel.N,mymodel.q,arma::fill::zeros);
  arma::cube change_J(mymodel.q,mymodel.q,mymodel.N*(mymodel.N-1)/2,arma::fill::zeros);

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
    double diff1 = sqrt(arma::accu(arma::square(dh_prev))/mymodel.mom1.n_elem);
    double diff2 = sqrt(arma::accu(arma::square(dJ_prev))/mymodel.mom2.n_elem);
    double max1 = arma::abs(dh_prev).max();
    double max2 = arma::abs(dJ_prev).max();
    std::cout << "RMS of dh, dJ: " << diff1 << " " << diff2 << std::endl;
    std::cout << "Max of dh, dJ: " << max1 << " " << max2 << std::endl;
    if(mymodel.N==2){
      std::cout << "h:" << std::endl;
      std::cout << mymodel.h << std::endl;
      std::cout << "J:" << std::endl;
      std::cout << mymodel.J << std::endl;
      std::cout << "dh:" << std::endl;
      std::cout << dh << std::endl;
      std::cout << "dJ:" << std::endl;
      std::cout << dJ << std::endl;
      std::cout << "model 1p freq:" << std::endl;
      std::cout << mymodel.mom1 << std::endl;
      std::cout << "model 2p freq:" << std::endl;
      std::cout << mymodel.mom2 << std::endl;
    }

    if(adaptive_sampling_on){
      if(diff1<mymodel.mom1_err && diff2<mymodel.mom2_err){
        std::cout << "Std. errors: " << mymodel.mom1_err << " " << mymodel.mom2_err << std::endl;
        std::cout << "Model difference is within sampling error. Increasing sampling..." << std::endl;
        mc_steps *= 2;
        std::cout << "MC steps per iteration: " << mc_steps << std::endl;
      }
    }
/*
    if(mymodel.q==2){
//      double Z = mymodel.get_Z();
//      entropy[niter] = log(Z)-arma::accu(mymodel.mom1%mymodel.h)-arma::accu(mymodel.mom2%mymodel.J); 
      std::cout << "Cross-entropy: " << entropy[niter] << std::endl;
    }
*/
    //Dump current statistics and parameters to file
    if (verbose) {
	    mymodel.mom1.save(scratch_dir + "stat_MC_1p_" + std::to_string(niter) + ".txt", arma::arma_ascii);
	    mymodel.mom2.save(scratch_dir + "stat_MC_2p_" + std::to_string(niter) + ".txt",arma::arma_ascii);
	    mymodel.h.save(scratch_dir + "h_" + std::to_string(niter) + ".txt", arma::arma_ascii);
	    mymodel.J.save(scratch_dir + "J_" + std::to_string(niter) + ".txt", arma::arma_ascii);
    }
    

    /*** Compute derivatives ***/
    dh = (msa_freq-mymodel.mom1 - 2*mymodel.lambda*mymodel.h);
    dJ = (msa_corr-mymodel.mom2 - 2*mymodel.lambda*mymodel.J);

    change_h = gamma_mom*change_h + alpha_h%dh;
    change_J = gamma_mom*change_J + alpha_J%dJ;

    //Update parameters
    mymodel.h += change_h;
    mymodel.J += change_J;


   //Dump derivatives to file
    if (verbose) {
	    dh.save(scratch_dir + "dh_" + std::to_string(niter) + ".txt", arma::arma_ascii);
	    dJ.save(scratch_dir + "dJ_" + std::to_string(niter) + ".txt", arma::arma_ascii);
    }


    //Update learning rates **MAKE THIS INTO A SEPARATE FUNCTION**
    if(adaptive_stepsize_on){
      for(int i=0; i<mymodel.N; i++){
        for(int a=0; a<mymodel.q; a++){
          double prod = dh_prev(i,a)*dh(i,a);
          if(prod>0 && alpha_h(i,a)<eps_max_h) alpha_h(i,a) = eps_inc*alpha_h(i,a);
          else if(prod<0 && alpha_h(i,a)>eps_min_h) alpha_h(i,a) = eps_dec*alpha_h(i,a);
          else {};
          if(alpha_h(i,a)>eps_max_h) alpha_h(i,a)=eps_max_h;
          if(alpha_h(i,a)<eps_min_h) alpha_h(i,a)=eps_min_h;
        }
      }

      for(int i=0; i<mymodel.N-1; i++){
        for(int j=i+1; j<mymodel.N; j++){ 
          int index = (mymodel.N-1)*i-i*(i+1)/2+j-1;
          for(int a=0; a<mymodel.q; a++){
            for(int b=0; b<mymodel.q; b++){
              double prod = dJ_prev(a,b,index)*dJ(a,b,index);
              if(prod>0 && alpha_J(a,b,index)<eps_max_J_N/mymodel.N) alpha_J(a,b,index) = eps_inc*alpha_J(a,b,index);
              else if(prod<0 && alpha_J(a,b,index)>eps_min_J) alpha_J(a,b,index) = eps_dec*alpha_J(a,b,index);
              else {};
              if(alpha_J(a,b,index)>eps_max_J_N/mymodel.N) alpha_J(a,b,index)=eps_max_J_N/mymodel.N;
              if(alpha_J(a,b,index)<eps_min_J) alpha_J(a,b,index)=eps_min_J;
            }
          }
        }
      }
    }


    //Update previous step gradient
    dh_prev = dh;
    dJ_prev = dJ;
    
    niter++;
    if(niter>max_iter || (max1<cutoff_freq && max2<cutoff_freq)){
      std::cout << "Converged!" << std::endl;
      std::cout << "RMS of dh, dJ: " << diff1 << " " << diff2 << std::endl;
      std::cout << "Max of |dh|, |dJ|: " << max1 << " " << max2 << std::endl;
      converged=true;
    }

  }

  //Dump energies to file
  avg_energies.save(output_dir + "avg_ene.txt", arma::arma_ascii);
  //entropy.save(output_dir + "entropy.txt", arma::arma_ascii);

}

void read_inputs(std::string name){

  // default values...
  nrep = 1;
  verbose = 0;

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
      if(key.compare("scratch_dir")==0) scratch_dir = value;
      if(key.compare("lambda")==0) lambda = std::stof(value);
      if(key.compare("myseed")==0) myseed = (long unsigned)std::stoi(value);
      if(key.compare("max_iter")==0) max_iter = std::stoi(value);
      if(key.compare("mc_steps")==0) mc_steps = std::stoi(value);
      if(key.compare("eps0_h")==0) eps0_h = std::stof(value);
      if(key.compare("eps0_J")==0) eps0_J = std::stof(value);
      if(key.compare("adaptive_stepsize_on")==0) adaptive_stepsize_on = std::stoi(value);
      if(key.compare("adaptive_sampling_on")==0) adaptive_sampling_on = std::stoi(value);
      if(key.compare("eps_inc")==0) eps_inc = std::stof(value);
      if(key.compare("eps_dec")==0) eps_dec = std::stof(value);
      if(key.compare("eps_min_h")==0) eps_min_h = std::stof(value);
      if(key.compare("eps_max_h")==0) eps_max_h = std::stof(value);
      if(key.compare("eps_min_J")==0) eps_min_J = std::stof(value);
      if(key.compare("eps_max_J_N")==0) eps_max_J_N = std::stof(value);
      if(key.compare("gamma")==0) gamma_mom = std::stof(value);
      if(key.compare("cutoff_freq")==0) cutoff_freq = std::stof(value);
      if(key.compare("delta")==0) delta = std::stof(value);
      if(key.compare("mc_init")==0) mc_init = value;
      if(key.compare("symmetrize_on")==0) symmetrize_on = std::stoi(value);
      if(key.compare("nrep")==0) nrep = std::stoi(value);
      if(key.compare("verbose")==0) verbose=std::stoi(value);
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
      if(id=="J"){
        int i;
        int j;
        int a;
        int b;
        my_stream >> i;
        my_stream >> j;
        my_stream >> a;
        my_stream >> b;
        int index = (mymodel.N-1)*i-i*(i+1)/2+j-1; 
        my_stream >> mymodel.J(a,b,index);
      }
      else if(id=="h"){
        int i;
        int a;
        my_stream >> i;
        my_stream >> a;
        my_stream >> mymodel.h(i,a);
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
          ofile << "J " << i << " " << j << " " << a << " " << b << " " << mymodel.J(a,b,index) << std::endl;
        }
      }
    }
  }

  for(int i=0; i<mymodel.N; i++){
    for(int a=0; a<mymodel.q; a++){
      ofile << "h " << i << " " << a << " " << mymodel.h(i,a) << std::endl;
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
