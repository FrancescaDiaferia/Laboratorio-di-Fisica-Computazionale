#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N_args 4


double k, m, T_max, omega2, E_0, x_0, v_0, choice; //variabili globali
FILE *fp_e, *fp_ec, *fp_p, *fp_lf, *fp_vv, *fp_vp, *fp_rk4, *fp_T;

typedef struct{
  double x;
  double v;
}vec;


void leggi_input(void);
void eulero(double dt, int N);
void eulero_cromer(double dt, int N);
void pnt_centrale(double dt, int N);
void leap_frog(double dt, int N);
void verlet_v(double dt, int N);
void verlet_p(double dt, int N);
double energy(double vel, double pos);
double a(double pos, double vel, double t);
void choose_algo(char *a);
void usage(void);
void trajectory(double traiettoria[], double dt, int N);
vec somma(vec a, vec b);
vec prodSc(vec a, double k);
vec f(vec y, double t);
void RK4(double dt, int N);

int main(int argc, char *argv[]){
  double dt, n_dt;
  int i, N;

  if(argc!=N_args){
    fprintf(stderr, "\nErrore!\n");
    usage();
    fprintf(stderr, "Riesegui.\n");
    exit(EXIT_FAILURE);
  }


  fprintf(stderr,"Benvenuto! Sto eseguendo %s.\n", argv[0]);
  usage();
  leggi_input();



  dt=atof(argv[2]); //questo sarà il tempo d'integrazione massimo
  n_dt=atoi(argv[3]); //numero di passi d'integrazione
  N=(int) (T_max/dt +0.5); //+0.5 affinché approssimi meglio

  choose_algo(argv[1]);

  for(i=0; i<n_dt; i++){
    if(choice==1){
      eulero(dt, N);
    }else if(choice==2){
      eulero_cromer(dt, N);
    }else if(choice==3){
      pnt_centrale(dt, N);
    }else if(choice==4){
      leap_frog(dt, N);
    }else if(choice==5){
      verlet_v(dt, N);
    }else if(choice==6){
      verlet_p(dt, N);
    }else if(choice==7){
      RK4(dt, N);
    }else{
      fprintf(stderr, "\nErrore1!\n");
      usage();
      fprintf(stderr, "Riesegui.\n");
      exit(EXIT_FAILURE);
    }

    N=N*2;
    dt=dt*0.5;
  }

  fclose(fp_e);
  fclose(fp_ec);
  fclose(fp_p);
  fclose(fp_lf);
  fclose(fp_vp);
  fclose(fp_vv);
  fclose(fp_rk4);
  fclose(fp_T);

  return EXIT_SUCCESS;
}

void leggi_input(void){


  fprintf(stderr,"Costante k?\n");
  scanf("k %lf\n", &k);
  fprintf(stderr,"Massa m?\n");
  scanf("m %lf\n", &m);
  omega2=k/m;
  fprintf(stderr,"La frequenza pertanto vale:%g\n", sqrt(omega2));
  fprintf(stderr,"Posizione iniziale x_0?\n");
  scanf("x_0 %lf\n", &x_0);
  fprintf(stderr,"Velocita' iniziale v_0?\n");
  scanf("v_0 %lf\n", &v_0);
  fprintf(stderr,"Tempo massimo di integrazione T_max?\n");
  scanf("T_max %lf\n", &T_max);
  E_0=energy(v_0,x_0);
  fprintf(stderr,"L'energia iniziale e': %g\n", E_0);
  fprintf(stderr,"Hai scelto: k=%g, m=%g, x_0=%g, v_0=%g, T_max=%g\n", k, m, x_0, v_0, T_max);


}

inline double energy(double vel, double pos){
    return 0.5*m*vel*vel + 0.5*k*pos*pos;
}

inline double a(double pos, double vel, double t){
  return -omega2 * pos;
  //return 0; //per controllare che funzioni bene l'integrazione
}



void eulero(double dt, int N){
	int i;
	double t, x_temp, E, delta_Er, x, v;
	char nome[100];
	FILE *fp;
  double *x_save;

	x=x_0;
	v=v_0;

	sprintf(nome, "eulero_dt%g.out", dt);
	fp=fopen(nome, "w");

	fprintf(fp, "#dt=%g\n", dt);
	fprintf(fp,"#tempo posizione velocità energia delta_E \n");
	fprintf(fp,"#0  %g %g %g \n", x, v, E_0);

  if( (x_save=(double *)calloc(N+1, sizeof(double)) )==NULL){
    exit(EXIT_FAILURE);
  }

	for(i=1; i<=N; i++){

		t= i * dt;
		x_temp=x;
		x= x + v*dt;
    x_save[i-1]=x;
		v= v + a(x_temp, v, t) * dt;
		E=energy(v,x);
		delta_Er=fabs(E-E_0) / E_0;

		fprintf(fp,"%g  %g %g %g %g\n", t, x, v, E, delta_Er);
	}

  trajectory(x_save, dt, N);

  fprintf(fp_e, "%g %g\n", dt, delta_Er);

	fclose(fp);
}

void eulero_cromer(double dt, int N){

	int i;
	double t, E, delta_Er, x, v;
	char nome[100];
	FILE *fp;
  double *x_save;

	x=x_0;
	v=v_0;

	sprintf(nome, "eulerocromer_dt%g.out", dt);
	fp=fopen(nome, "w");

	fprintf(fp, "#dt=%g\n", dt);
	fprintf(fp,"# tempo posizione velocità energia delta_E \n");
	fprintf(fp,"#0  %g %g %g  \n", x, v, E_0);

  if( (x_save=(double *)calloc(N+1, sizeof(double)) )==NULL){
    exit(EXIT_FAILURE);
  }


	for(i=1; i<=N; i++){

		t= i * dt;

		v= v + a(x,v,t) * dt;
		x= x + v * dt;
    x_save[i-1]=x;

		E=energy(v,x);
		delta_Er=fabs(E-E_0) / E_0;

		fprintf(fp,"%g  %g %g %g %g\n", t, x, v, E, delta_Er);

	}

  trajectory(x_save, dt, N);

  fprintf(fp_ec, "%g %g\n", dt, delta_Er);
	fclose(fp);
}

void choose_algo(char *a){
  choice= atoi(a);

  if(choice==1){
    fp_e=fopen("tail_e.dat", "w");
    fp_T=fopen("eu_semiperiodi.dat", "w");
    fprintf(stderr, "Hai scelto di integrare con Eulero.\n");

  }else if(choice==2){
    fp_ec=fopen("tail_ec.dat", "w");
    fp_T=fopen("ec_semiperiodi.dat", "w");
    fprintf(stderr, "Hai scelto di integrare con Eulero-Cromer.\n");

  }else if(choice==3){
    fp_p=fopen("tail_p.dat", "w");
    fp_T=fopen("pc_semiperiodi.dat", "w");
    fprintf(stderr, "Hai scelto di integrare con punto centrale.\n");

  }else if(choice==4){
    fp_lf=fopen("tail_lf.dat", "w");
      fp_T=fopen("lf_semiperiodi.dat", "w");
    fprintf(stderr, "Hai scelto di integrare con leap-frog.\n");

  }else if(choice==5){
    fp_vv=fopen("tail_vv.dat", "w");
    fp_T=fopen("vv_semiperiodi.dat", "w");
    fprintf(stderr, "Hai scelto di integrare con Verlet-velocity.\n");
  }else if(choice==6){
    fp_vp=fopen("tail_vp.dat", "w");
    fp_T=fopen("vp_semiperiodi.dat", "w");
    fprintf(stderr, "Hai scelto di integrare con Verlet-position.\n");
  }else if(choice==7){
    fp_rk4=fopen("tail_rk4.dat", "w");
    fp_T=fopen("rk4_semiperiodi.dat", "w");
    fprintf(stderr, "Hai scelto di integrare con RK4.\n");

  }else{
    fprintf(stderr, "\nErrore!\n");
    usage();
    fprintf(stderr, "Riesegui.\n");
    exit(EXIT_FAILURE);
  }

}

void usage(){
  fprintf(stderr, "Questo algoritmo deve ricevere in input %i argomenti:\n1)<nome.exe>\n2)algoritmo d'integrazione:<1>=eulero, <2>=eulero-cromer, <3>=punto centrale, <4>=leap frog, <5>=verlet velocity, <6>=verlet position, <7>=RK4\n3)<dt>=passo d'integrazione MASSIMO \n4)<n_dt>=numero di passi d'integrazione\n", N_args);
}

void pnt_centrale(double dt, int N){

	int i;
	double t, E, delta_Er, x, v, v_temp;
	char nome[100];
	FILE *fp;
  double *x_save;

	x=x_0;
	v=v_0;

	sprintf(nome, "pntcentrale_dt%g.out", dt);
	fp=fopen(nome, "w");

	fprintf(fp, "#dt=%g\n", dt);
	fprintf(fp,"# tempo posizione velocità energia delta_E \n");
  fprintf(fp,"#0  %g %g %g  \n", x, v, E_0);

  if( (x_save=(double *)calloc(N+1, sizeof(double)) )==NULL){
    exit(EXIT_FAILURE);
  }

	for(i=1; i<=N; i++){

	  t= i * dt;

	  v_temp=v;
	  v= v + a(x,v,t)*dt;
	  x= x + ((v_temp+v)*0.5) * dt;
    x_save[i-1]=x;

	  E=energy(x,v);
	  delta_Er=fabs(E-E_0) / E_0;

	  fprintf(fp,"%g  %g %g %g %g\n", t, x, v, E, delta_Er);

	}

  trajectory(x_save, dt, N);

  fprintf(fp_p, "%g %g\n", dt, delta_Er);
	fclose(fp);
}

void leap_frog(double dt, int N){

	int i;
	double t, E, delta_Er, x, v;
	char nome[100];
	FILE *fp;
  double *x_save;

	x=x_0;
	v=v_0;
  t=0.;

	sprintf(nome, "leapfrog_dt%g.out", dt);
	fp=fopen(nome, "w");

	fprintf(fp, "#dt=%g\n", dt);
	fprintf(fp,"# tempo posizione velocità energia delta_E \n");
	fprintf(fp,"0  %g %g %g  \n", x, v, E_0);

  if( (x_save=(double *)calloc(N+1, sizeof(double)) )==NULL){
    exit(EXIT_FAILURE);
  }

	v = v + a(x,v,t) * (0.5*dt);

	for(i=1; i<=N; i++){

	  t= i * dt;

	  x= x + v*dt;
    x_save[i-1]=x;
	  v= v + a(x,v,t)*dt;

	  E=energy(x, (v- a(x, v, t)* (dt*0.5))); //affinché siamo x e v dello stesso tempo
	  delta_Er=fabs(E-E_0)/E_0;

	  fprintf(fp,"%g  %g %g %g %g\n", t, x, (v-a(x,v,t)*(dt*0.5)), E, delta_Er);
	}

	v=v- a(x,v,t)*(dt*0.5);
  trajectory(x_save, dt, N);

  fprintf(fp_lf, "%g %g\n", dt, delta_Er);
	fclose(fp);
}


void verlet_v(double dt, int N){
  int i;
	double t, E, delta_Er, x, v, x_temp;
	char nome[100];
	FILE *fp;
  double *x_save;

	x=x_0;
	v=v_0;

	sprintf(nome, "verlet_v_dt%g.out", dt);
	fp=fopen(nome, "w");

  fprintf(fp, "#dt=%g\n", dt);
	fprintf(fp,"# tempo posizione velocità energia delta_E \n");
	fprintf(fp,"0  %g %g %g  \n", x, v, E_0);

  if( (x_save=(double *)calloc(N+1, sizeof(double)) )==NULL){
    exit(EXIT_FAILURE);
  }

  for(i=1; i<=N; i++){
    t= i * dt;
    x_temp=x;
    x=x + v*dt + 0.5*a(x,v,t)*dt*dt;
    x_save[i-1]=x;
    v=v + ( a(x_temp, v, t) + a(x, v, t) )*0.5*dt;
    E=energy(x, v);

    delta_Er=fabs(E-E_0)/E_0;

	  fprintf(fp,"%g  %g %g %g %g\n", t, x, v, E, delta_Er);
  }

  trajectory(x_save, dt, N);

  fprintf(fp_vv, "%g %g\n", dt, delta_Er);
  fclose(fp);

}

void verlet_p(double dt, int N){
  int i;
	double t, E, delta_Er, x, v, x_old, x_new;
	char nome[100];
	FILE *fp;
  double *x_save;

	sprintf(nome, "verlet_p_dt%g.out", dt);
	fp=fopen(nome, "w");

  fprintf(fp, "#dt=%g\n", dt);
	fprintf(fp,"# tempo posizione velocità energia delta_E \n");
	fprintf(fp,"0  %g %g %g  \n", x_0, v_0, E_0);

  if( (x_save=(double *)calloc(N+1, sizeof(double)) )==NULL){
    exit(EXIT_FAILURE);
  }
  t=0.;
  x_old=x_0;
  v=v_0;
  x= x_0 + v_0*dt + 0.5*a(x_0, v_0, t)*dt*dt;
  for(i=1; i<=N; i++){
    t= i * dt;

    x_save[i-1]=x;
    x_new=2 * x-x_old +a(x, v, t)*dt*dt; //posizione al tempo x(n+1)
    v=(x_new-x_old) / (2*dt); //velocità al tempo v(n)

    E=energy(x, v);
    delta_Er=fabs(E-E_0)/E_0;

    fprintf(fp,"%g  %g %g %g %g\n", t, x, v, E, delta_Er);

    x_old=x;
    x=x_new;
  }


  trajectory(x_save, dt, N);

  fprintf(fp_vp, "%g %g\n", dt, delta_Er);
  fclose(fp);

}



void trajectory(double traiettoria[], double dt, int N){
  int j;
  double T_half, t1, t2;
  double t_star, t_temp;

  fprintf(fp_T, "#dt=%g\n", dt);
  fprintf(fp_T, "#T/2(atteso): %g\n", (M_PI/sqrt(omega2)));
  fprintf(fp_T,"#T/2 \t T/2-T/2(atteso)\n ");

  t_temp=-1.0;
  for(j=0; j<N+1; j++){
        if(traiettoria[j]*traiettoria[j+1]<0){
          t1=(j+1)*dt; //perchè il ciclo inizia da 0
          t2=(j+2)*dt;
          t_star= t1 - (traiettoria[j]* (t2-t1) )/(traiettoria[j+1]-traiettoria[j]);
          if(t_temp!=-1.0){
            T_half=t_star-t_temp;
            fprintf(fp_T, "%g \t %lg\n", T_half, ((M_PI/sqrt(omega2))-T_half));
          }
          t_temp=t_star;
        }
      }

     fprintf(fp_T, "\n\n");

}

vec somma(vec a, vec b){

  vec res;

    res.x=a.x+b.x;
    res.v=a.v+b.v;

  return res;
}

vec prodSc(vec a, double k){
  vec res;

    res.x=a.x*k;
    res.v=a.v*k;

  return res;
}

vec f(vec y, double t){
  vec dy;

    dy.x=y.v;  //velocità
    dy.v=a(y.x, y.v, t); //accelerazione

  return dy;

}


void RK4(double dt, int N){
  int i;
  double t, E, delta_Er;
  vec y, k1, k2, k3, k4;
  FILE *fp;
  char nome[100];
  double *x_save;

  y.x=x_0;
  y.v=v_0;

  sprintf(nome, "rk4_dt%g.out", dt);
	fp=fopen(nome, "w");

  fprintf(fp, "#dt=%g\n", dt);
  fprintf(fp,"# tempo posizione velocità energia delta_E \n");
  fprintf(fp,"0  %g %g %g  \n", y.x, y.v, E_0);

  if( (x_save=(double *)calloc(N+1, sizeof(double)) )==NULL){
    exit(EXIT_FAILURE);
  }

  for(i=1; i<=N; i++){

       t=i*dt;
       k1=f(y, t);
       k2=f( (somma(y, prodSc(k1,0.5*dt) ) ), (t+dt*0.5) );
       k3=f( (somma(y, prodSc(k2,0.5*dt) ) ), (t+dt*0.5) );
       k4=f( (somma(y, prodSc(k3,dt) ) ), (t+dt) );

       y= somma(y, prodSc( (somma( somma( somma(k1, prodSc(k2,2.0)) , prodSc(k3,2.0)) , k4)), (dt/6.)) );
       x_save[i-1]=y.x;
       E=energy(y.x, y.v);
       delta_Er=fabs(E-E_0)/E_0;

      fprintf(fp,"%g  %g %g %g %g\n", t, y.x, y.v, E, delta_Er);
  }
  fprintf(fp_rk4, "%g %g\n", dt, delta_Er);

  trajectory(x_save, dt, N);

  fclose(fp);
}
