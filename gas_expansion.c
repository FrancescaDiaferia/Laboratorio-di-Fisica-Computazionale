//PROGRAMMA SCRITTO SU WINDOWS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N_args 4
#define N_iterations 10

typedef struct{
  int i_x, i_y, dist_x, dist_y;
}particle;

void usage(void);
void fill(particle *part, int **grid, int *bucket);
void moving(particle *part, int **grid, int *bucket);
void clearArray(int **grid);

FILE *fp, *fp2;
int length=0, N_particles=0, N_passi=0;
double R2 = 0.;

int main(int argc, char *argv[]){
  int seed;
  int i=0;
  particle *part;
  int **grid, *bucket, Vol;

  if(argc!=N_args){
    fprintf(stderr, "\nErrore nell'inserimento degli argomenti!\n");
    usage();
    fprintf(stderr, "Riesegui.\n");
    exit(EXIT_FAILURE);
  }

  seed = atoi(argv[1]);
  srand(seed);

  length = atoi(argv[2]);
  N_passi = atoi(argv[3]);

  Vol = length * length;



  if( (bucket=(int*)calloc((Vol), sizeof(int)) )==NULL){
    exit(EXIT_FAILURE);
  }

  if( (grid = (int **) calloc(length, sizeof(int *)) )==NULL){
    exit(EXIT_FAILURE);
  }
  for(i=0; i<length; i++){
    if ( (grid[i] = (int *) calloc(length, sizeof(int)) )==NULL) {
      exit(EXIT_FAILURE);
    }
  }


  fp=fopen("diffusione.dat", "w");
  fprintf(fp, "#1:N\t 2:L\t 3:t\t 4:R^2\n");

  fp2=fopen("alltimes.dat", "w");
  fprintf(fp2, "#1:N\t 2:L\t 3:t\t 4:R^2\n");

/*Itero il procedimento per stimare al meglio il coefficiente di diffusione D (la media verrà fatta su gnuplot con smooth unique)*/
  for(i=0; i<N_iterations; i++){

    /*Inizializzo tutto l'array con componenti -1 così che quando gli assegnerò di volta in volta le posizioni riempite posso distringuere le caselle vuote da quelle occupate*/
    clearArray(grid);

/*Faccio un for sul numero di particelle in modo tale da avere le densità da 0.1 a 0.9 con passo 0.1*/
      for(N_particles= (Vol *0.1); N_particles<=Vol; N_particles+=(Vol * 0.1)){

//Ogni volta lo alloco dinamicamente col nuovo numero di particelle
        if( (part=(particle*)calloc(N_particles, sizeof(particle)) )==NULL){
          exit(EXIT_FAILURE);
        }

          clearArray(grid); //riporto gli array alle condizioni iniziali

          fill(part, grid, bucket);

          moving(part, grid, bucket);

          fprintf(fp2, "\n\n");

          free(part);
      }
  }

  for(i=0; i< length; i++){
    free(grid[i]);
  }
  free(grid);

  free(bucket);


  fclose(fp);

  return EXIT_SUCCESS;
}

void usage(void){
  fprintf(stderr, "Questo programma riceve in input %i argomenti: nome, seed, L, numero di passi.\n", N_args);
}

void fill(particle *part, int **grid, int *bucket){
  int num=0, k=0, n=0, rnd=0;

  num = length * length;
  for(k=0; k<num; k++){
    bucket[k] = k;
  }

  for(n = 0; n < N_particles; n++){

    rnd = (int) ((rand() / (RAND_MAX + 1.) ) * num);

    /*Assegno alla casella dell'array scelta casualmente, tramite estrazione delle sue coordinate, la posizione estratta*/
    grid[(int) (bucket[rnd]/length)][bucket[rnd]%length]=n;
    part[n].i_x = bucket[rnd] / length;
    part[n].i_y = bucket[rnd] % length;

    part[n].dist_x = 0;
    part[n].dist_y = 0;

    bucket[rnd]=bucket[--num];
  }


}

void moving(particle *part, int **grid, int *bucket){
  int t=0, n=0;
  double rnd;

  for(t=1; t<N_passi; t++){
    R2 = 0.;

    for(n=0; n<N_particles; n++){
      rnd=rand() / (RAND_MAX + 1.);

          if( rnd < 0.25 ) {
            //vado a destra se è vuota la posizione (c'è -1)
            if(grid[(part[n].i_x+1) % length][part[n].i_y]==-1){
              //sposto la particella a destra e libero la posizione che lascio
              grid[(part[n].i_x+1) % length][part[n].i_y] = n;
              grid[(part[n].i_x)][part[n].i_y] = -1;

              //aggiorno part con la nuova mossa
              part[n].i_x = (part[n].i_x+1) % length;
              part[n].dist_x+=1;
            }
          }else if( rnd < 0.5 ){
            //vado in alto se è vuota la posizione (c'è -1)
            if(grid[part[n].i_x][(part[n].i_y+1) % length]==-1){
              //sposto la particella in alto e libero la posizione che lascio
              grid[part[n].i_x][(part[n].i_y+1) % length] = n;
              grid[part[n].i_x][part[n].i_y] = -1;

              //aggiorno part con la nuova mossa
              part[n].i_y = (part[n].i_y+1) % length;
              part[n].dist_y+=1;
            }
          }else if( rnd < 0.75 ){
            //vado a sinistra se è vuota la posizione (c'è -1)
            if(grid[(part[n].i_x+length-1) % length][part[n].i_y]==-1){
              //sposto la particella a sinistra e libero la posizione che lascio
              grid[(part[n].i_x+length-1) % length][part[n].i_y] = n;
              grid[(part[n].i_x)][part[n].i_y] = -1;

              //aggiorno part con la nuova mossa
              part[n].i_x = (part[n].i_x+length-1) % length;
              part[n].dist_x-=1;
            }
          }else{
            //vado in basso se è vuota la posizione (c'è -1)
            if(grid[part[n].i_x][(part[n].i_y+length-1) % length]==-1){
              //sposto la particella in basso e libero la posizione che lascio
              grid[part[n].i_x][(part[n].i_y+length-1) % length] = n;
              grid[(part[n].i_x)][part[n].i_y] = -1;

              //aggiorno part con la nuova mossa
              part[n].i_y = (part[n].i_y+length-1) % length;
              part[n].dist_y-=1;
            }
          }
          R2 += part[n].dist_x * part[n].dist_x + part[n].dist_y * part[n].dist_y;
    }
    R2 /= (double) N_particles;
    fprintf(fp2, "%i\t %i\t %i\t %lf\n", N_particles, length, t, R2);

  }
  fprintf(fp, "%i\t %i\t %i\t %g\n", N_particles, length, t, R2);

}

void clearArray(int **grid){
  int i, j;

  for(i=0; i<length; i++){
    for(j=0; j<length; j++){
      grid[i][j] = -1;
    }
  }

}
