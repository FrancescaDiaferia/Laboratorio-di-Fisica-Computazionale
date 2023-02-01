//ATTENZIONE
/*Questo programma è stato scritto su windows, pertanto ad esempio è stato necessario usare lo specificatore dei long long int adatto a questo sistema operativo e bisogna inserire il seed da linea di comando non esendo presente il file devran*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N_args 5
#define print_traj 10

FILE *fp, *fp_2;

int find_Nhisto(int b, int t_max);

int main(int argc, char *argv[]){
  long long int seed, x=0;
  int i=0, t=0, num_traj=0, lenght=0, N_histo, base, tmp_base, k, j;
  long long int *sum1, *sum2, *sum4;
  int **counters;
  double mean_1=0., mean_2=0., mean_4=0., sigma_mean1=0., sigma_mean2=0.;

  if(argc!=N_args){
    fprintf(stderr, "\nErrore nell'inserimento degli argomenti!\n");
    fprintf(stderr, "Questo programma riceve in input %i argomenti: <%s>, <seed>, <L> lunghezza traiettoria, <N> numero di traiettorie, <b> base.\n", N_args, argv[0]);
  }
    fprintf(stderr, "Riesegui.\n");
    exit(EXIT_FAILURE);
  }

  /*Su windows non funziona però
  FILE* devran;
  if((devran=fopen("/dev/urandom", "r")) == NULL)
  {
    fprintf(stderr, "Errore nell'apertura di un file.\n");
    return EXIT_FAILURE;
  }
  fread(&seed, sizeof(seed), 1, devran);
  fclose(devran);
  */

  seed = atoi(argv[1]);
  srand(seed);

  lenght = atoi(argv[2]);
  num_traj = atoi(argv[3]);
  base = atoi(argv[4]);

  N_histo=find_Nhisto(base, lenght);


  if( (sum1=(long long int *)calloc(lenght+1, sizeof(long long int)) )==NULL){
    exit(EXIT_FAILURE);
  }

  if( (sum2=(long long int *)calloc(lenght+1, sizeof(long long int)) )==NULL){
    exit(EXIT_FAILURE);
  }

  if( (sum4=(long long int *)calloc(lenght+1, sizeof(long long int)) )==NULL){
    exit(EXIT_FAILURE);
  }

  counters = (int **) calloc(N_histo, sizeof(int *));
  fp_2=fopen("histo.dat", "w");


  tmp_base = base;
  for(k=0; k<N_histo; k++){
    counters[k] = (int *) calloc( (tmp_base*2 + 1) , sizeof(int));
    tmp_base*=base;
  }

  fp=fopen("medie.dat", "w");
  fprintf(fp, "#t\tmean_x\tsigma_mean_x\tmean_x^2\tsigma_mean_x^2\n ");


  for(i=0; i<num_traj; i++){
    if(i<print_traj){
      printf("#indice\tx\n");
    }
    x = 0;
    k = 0;
    tmp_base = base;

    for(t=0; t<=lenght; t++){
      if(i<print_traj){
        printf("%i\t%I64d\n", t, x);
      }


      if(t==tmp_base){
        (counters[k][x + t])++;
        k++;
        tmp_base *= base;
      }

      sum1[t] += x;
      sum2[t] += x * x;
      sum4[t] += x * x * x * x;

      if( (rand() / (RAND_MAX + 1.) ) > 0.5){
        x+=1;
      }else{
        x-=1;
      }

    }
    if(i<print_traj){
      printf("\n\n");
    }

  }

  for(t=0; t<=lenght; t++){
    mean_1 = ((double)sum1[t] / num_traj);
    mean_2 = ((double)sum2[t] / num_traj);
    mean_4 = ((double)sum4[t] / num_traj);
    sigma_mean1 = sqrt( (mean_2 - mean_1*mean_1) / (num_traj - 1) );
    sigma_mean2 = sqrt( (mean_4 - mean_2*mean_2) / (num_traj - 1) );
    fprintf(fp, "%d\t%g\t%g\t%g\t%g\n", t, mean_1, sigma_mean1, mean_2, sigma_mean2);


      tmp_base = base;
      for(i=0; i<N_histo; i++){
        if(t==tmp_base){
          fprintf(fp_2, "#%i\n", tmp_base);
          for(j=0; j<tmp_base*2 +1; j++){
              fprintf(fp_2, "%i\t%lf\n", (j-t), (counters[i][j]/(double)num_traj) );
            }
          fprintf(fp_2, "\n\n");
        }
        tmp_base*=base;
      }

  }


  //I64d specificatore di long long int

  fclose(fp);
  fclose(fp_2);




  for(k=0; k<N_histo; k++){
    free(counters[k]);
  }
  free(counters);

  return EXIT_SUCCESS;
}


int find_Nhisto(int b, int t_max){
  int Nhisto, tmp;

  tmp = b;
  Nhisto = 0;
  while(tmp<=t_max){
    Nhisto++;
    tmp*=b;
  }
  return Nhisto;
}
