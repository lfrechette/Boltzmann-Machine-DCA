#include "mc_io.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>

/*************************/
/*** Reading Functions ***/
/*************************/

void read_dims(const char *name, const char *method){

  char myfile[200];
  struct stat fileStat;
  sprintf(myfile, "../fitness/%s_params_%s.txt", name, method);
  if(stat(myfile, &fileStat)<0){
    printf("Error: parameter file does not exist. Exiting.\n");
    exit(-1);
  }
  FILE *f;
  f = fopen(myfile, "r");
  char *ignore = (char*)malloc(200); 
  if(strcmp(method, "plm")==0){
    if(fscanf(f, "%s %d\n", ignore, &N)!=2) printf("Error: require integer (sequence size.)\n");
    if(fscanf(f, "%s %d\n", ignore, &q)!=2) printf("Error: require integer (number of possible spin states.)\n");
  }
  else{
    int t1=0, t2=0, t3=0, t4=0;
    double d=0;
    q=0;
    while(q<=t4){
      q=t4;
      if(fscanf(f, "%s %d %d %d %d %lf\n", ignore, &t1, &t2, &t3, &t4, &d)!=6) printf("Error: 6 entries expected.\n");
    }
    q++;
    while(N<=t2){
      N=t2;
      if(fscanf(f, "%s %d %d %d %d %lf\n", ignore, &t1, &t2, &t3, &t4, &d)!=6) printf("Error: 6 entries expected.\n");
    }
    N++;
  }
  fclose(f);
  f=NULL;
  printf("N: %d\n", N);
  printf("q: %d\n", q);
}

void read_params(const char *name, const char *method, arma::mat &h, arma::cube &J){

  char myfile[200];
  struct stat fileStat;
  sprintf(myfile, "../fitness/%s_params_%s.txt", name, method);
  if(stat(myfile, &fileStat)<0){
    printf("Error: parameter file does not exist. Exiting.\n");
    exit(-1);
  }
  FILE *f;
  f = fopen(myfile, "r");
  char *ignore = (char*)malloc(1000); 
  if(strcmp(method, "plm")==0){
    if(fgets(ignore, (int)sizeof(ignore), f)==NULL) printf("nothing\n"); //skip line
    if(fgets(ignore, (int)sizeof(ignore), f)==NULL) printf("nothing\n"); //skip line
    if(fscanf(f, "%s\n", ignore)!=1) printf("Error: field header expected, not found.\n");
    for(int i=0;  i<N; i++){
      for(int j=0; j<(q-1); j++){
        if(fscanf(f, "%lf ", &h(i,j))!=1) printf("Error: expecting float.\n"); 
      }
      if(fscanf(f, "%lf\n", &h(i,q-1))!=1) printf("Error: expecting float.\n");
    }
    if(fscanf(f, "%s\n", ignore)!=1) printf("Error: coupling header expected, not found.\n");
    for(int i=0; i<(N*(N-1)/2); i++){
/*
      if(fgets(ignore, (int)sizeof(ignore), f)==NULL) printf("nothing\n"); //skip "pair" line
      std::cout << ignore << std::endl;
      if(fgets(ignore, (int)sizeof(ignore), f)==NULL) printf("nothing\n"); //skip "pair" line (need 2 for comma)
      std::cout << ignore << std::endl;
*/
      if(fscanf(f, "%*[^\n]\n")!=0) printf("not zero");
      /*
      if(strcmp(name, "rfah")==0){
        if(fgets(ignore, (int)sizeof(ignore), f)==NULL) printf("nothing\n");
        std::cout << ignore << std::endl;
      } 
      */
      for(int j=0; j<q; j++){
        for(int k=0; k<(q-1); k++){
          if(fscanf(f, "%lf ", &J(i,j,k))!=1){
            printf("Error: expecting float.\n");
            printf("%d\n", i);
            printf("%s\n", fgets(ignore, (int)sizeof(ignore), f));
            exit(0);
          }
        }
        if(fscanf(f, "%lf\n", &J(i,j,q-1))!=1) printf("Error: expecting float.n");
      } 
    }
  }
  else{
    int t1=0, t2=0, t3=0, t4=0;
    for(int i=0; i<(N-1); i++){
      for(int j=(i+1); j<N; j++){
        for(int q1=0; q1<q; q1++){
          for(int q2=0; q2<q; q2++){
            int l = (N-1)*i-i*(i+1)/2+j-1;
            if(fscanf(f, "%s %d %d %d %d %lf\n", ignore, &t1, &t2, &t3, &t4, &J(l,q1,q2))!=6) printf("Error: 6 entries expected.\n");
          }
        }
      }
    }
    for(int i=0; i<N; i++){
      for(int q1=0; q1<q; q1++){
        if(fscanf(f, "%s %d %d %lf\n", ignore, &t1, &t2, &h(i,q1))!=4) printf("Error: 4 entries expected.\n");
      }
    }
  }
  fclose(f);
  f=NULL;
}

/**************************/
/*** Printing Functions ***/
/**************************/

void print_traj_scalar(double *s, int n, const char *dir, const char *name){

  FILE *f;
  char dest[200];
  strcpy(dest, dir);
  strcat(dest, name);
  f = fopen(dest, "w");
  for(int i=0; i<n; i++){
    fprintf(f, "%d %f\n", i, s[i]);
  }
  fclose(f);
  f=NULL;
}

void pconf(int *seq, int time, const char *dir, const char *name){

  FILE *f;
  char dest[200];
  strcpy(dest, dir);
  strcat(dest, name);
  f = fopen(dest, "w");
  fprintf(f, "%d\n", N);
  fprintf(f, "%d\n", time);

  for(int i=0; i<N; i++){
    fprintf(f, "%d %d\n", i, seq[i]);
  }
  fclose(f);
  f=NULL;
}
