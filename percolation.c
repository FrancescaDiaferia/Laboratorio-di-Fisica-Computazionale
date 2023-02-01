//PROGRAMMA SCRITTO SU WINDOWS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define N_args 4

//variabili globali
int length, n_sample, N, num;
FILE *fp;
int *parent;
int counter_t = 0; //serve a contare il #rietichettamenti 

typedef struct{
  int is;
  int S;
}percolans;


void mergeCluster(int index, int j);
int truelabel(int i);
void label(int *bucket);
void relabel(int index);
void IsPerc(percolans *perc, int* size);

int main(int argc, char *argv[]){
  clock_t begin = clock();
  int seed, sum_S2=0, M=0; //M = #buchi
  int n=0, s=0, k=0; //indici dei for
  int *size, *bucket;
  char nome[100];
  percolans perc;
  double t_spent;

  if(argc!=N_args){
    fprintf(stderr, "\nErrore nell'inserimento degli argomenti!\nUsage: nome, seed, L, #samples\nRiesegui.\n");
    exit(EXIT_FAILURE);
  }

  seed = atoi(argv[1]);
  srand(seed);

  length = atoi(argv[2]);
  n_sample = atoi(argv[3]);
  N = length * length;
  num = N;


  if( (size=(int*)calloc(N+1, sizeof(int) ) )==NULL){
    exit(EXIT_FAILURE);
  }

  if( (bucket=(int*)calloc((length * length), sizeof(int)) )==NULL){
    exit(EXIT_FAILURE);
  }

  if( (parent=(int*)calloc(N, sizeof(int) ) )==NULL){
    exit(EXIT_FAILURE);
  }


  sprintf(nome, "perc3_L%d.dat", length);
  fp=fopen(nome, "w");
  fprintf(fp, "#1:L\t2:numbuchi\t3:isperc\t4:S_sperc\t5:sumS2\t6:t_spent\n");


//Ad ogni sample ricomincio da capo e riazzero tutto
for(s=0; s<n_sample; s++){
  counter_t = 0;
  num = N;

  for(k=0; k<N; k++){
    parent[k] = N;
    bucket[k] = k;
  }

  perc.is = 0;
  sum_S2 = 0;
  perc.S = 0;

//Aggiungo ogni volta 1 buco
  for(M=1; M<=N; M++){

    for(k=0; k<N+1; k++){
      size[k]=0; //azzero l'array
    }
    sum_S2 = 0;

    label(bucket);


  //conto quante volte compare ogni etichetta
      for(n=0; n<N; n++){
        if(parent[n]!=N){
          size[truelabel(n)]++;
        }
      }


  //Calcolo la somma della taglia dei cluster al quadrato
      for(n=0; n<N; n++){
        sum_S2+= (size[n]*size[n]);
      }

      IsPerc(&perc, size);


      fprintf(fp, "%i\t%i\t%i\t\t%i\t%i\t%i\n", length, M, perc.is, perc.S, sum_S2, counter_t);

  }

}
clock_t end = clock();
t_spent = (double) (end-begin)/CLOCKS_PER_SEC;
fprintf(stderr, "L'algoritmo ha impiegato %g s.\n", t_spent);

  free(size);
  free(bucket);
  free(parent);

  fclose(fp);
  return EXIT_SUCCESS;
}


void label(int *bucket){
  int rnd, index;

    rnd = (int) ((rand() / (RAND_MAX + 1.) ) * num);
    index = bucket[rnd];

    parent[index] = index;

    relabel(index);

    bucket[rnd]=bucket[--num];

}

int truelabel(int i){
  if(parent[i]==N){
    return N;
  }else if(parent[i]==i){
    return i;
  }else{
    return truelabel(parent[i]);
  }

/* //Versione non ricorsiva
  int lab=i;
  while(lab!=parent[lab]){
    lab = parent[lab];
  }
  return lab;
  */
}

void mergeCluster(int index, int j){
  counter_t ++;
  if(truelabel(index)<truelabel(j)){
    parent[truelabel(j)] = truelabel(index);
  }else{
    parent[truelabel(index)] = truelabel(j);
  }

}


void relabel(int index){


    //vado a destra
    if(index<N-1 && parent[index+1]!=N && truelabel(index)!=truelabel(index+1) ){
      mergeCluster(index, index+1);
    }
    //vado a sinistra
    if(index>0 && parent[index-1]!=N && truelabel(index)!=truelabel(index-1) ){
      mergeCluster(index, index-1);
    }
    //vado in alto
    if(index<N-length && parent[index+length]!=N && truelabel(index)!=truelabel(index+length) ){
      mergeCluster(index, index+length);
    }
    //vado in basso
    if(index>=length && parent[index-length]!=N && truelabel(index)!=truelabel(index-length) ){
      mergeCluster(index, index-length);
    }

}

void IsPerc(percolans *perc, int*size){
  int i;

  for(i=N-length; i<N; i++){

    if(truelabel(i) < length){
      perc->is = 1; //se percola
      perc->S=size[truelabel(i)]; //taglia del percolante
    }

  }
}
