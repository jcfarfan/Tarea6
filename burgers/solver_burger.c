#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float *reserva_matriz(int n_points_x, int n_points_y);
void *inicializar(float *matrix, int n_points_y);
void *imprimir(float *matrix, int n_points_x, int n_points_y);
void *unos(float *matrix, int n_points_x, int n_points_y);
void *escribirArchivo(float *matrix, int n_points_x, int n_points_y, int n);
void *copy(float *origen, float *destino, int n_points);
void *asignar(float *matrix, int n_points, float suma);
void *paso(float *u, float *v, float *un, float *vn, int n_points_x, int n_points_y);

int main(){
  //variables
  float *matrixU, *matrixV,*un, *vn;
  int n_points_x, n_points_y,i,j;
  n_points_x = 41; 
  n_points_y = 41;
  //matrices llenas de ceros 41x41
  matrixU = reserva_matriz(n_points_x,n_points_y);
  matrixV = reserva_matriz(n_points_x,n_points_y);
  un = reserva_matriz(n_points_x,n_points_y);
  vn = reserva_matriz(n_points_x,n_points_y);
  //matrices llenas de unos
  unos(matrixU,n_points_y,n_points_x);
  unos(matrixV,n_points_y,n_points_x);

  //Matrices inicializadas con las condiciones iniciales
  inicializar(matrixU,n_points_y);
  inicializar(matrixV,n_points_y);
  
  //paso
  paso(matrixU,matrixV,un,vn,n_points_x, n_points_y);
  
  return 0;
}


float *reserva_matriz(int n_points_x, int n_points_y){
  float *array;
  int i,j;
  if(!(array = malloc(n_points_x * n_points_y * sizeof(float)))){
    printf("Problema en reserva_matriz\n");
    exit(1);
  }
  for(i=0;i<n_points_x;i++){
    for(j=0;j<n_points_y;j++){
      array[i+(n_points_y*j)] = 0.0;
    }
  }
  return array;
}

void *inicializar(float *matrix, int n_points_y){
  int i,j;
  for(i=10;i<21;i++){
    for(j=10;j<21;j++){
      matrix[i+(n_points_y*j)] = 2.0;
    }
  }
}

void *imprimir(float *matrix, int n_points_x, int n_points_y){
  int i,j,pos;
  for(i=0;i<n_points_x;i++){
    for(j=0;j<n_points_y;j++){
      pos = j + (n_points_y * i);/*position in the array*/
      fprintf(stdout, " %f ",matrix[pos]);
    }
    fprintf(stdout, "\n");
  }
}

void *unos(float *matrix, int n_points_x, int n_points_y){
  int i,j;
  for(i=0;i<n_points_x;i++){
    for(j=0;j<n_points_y;j++){
      matrix[j+(n_points_y*i)] = 1.0;
    }
  }

}

void *escribirArchivo(float *matrix, int n_points_x, int n_points_y, int n){
  int i,j;
  FILE* archivo;
  char filename[100];
  sprintf(filename,"matrizdata%d.dat",n);
  archivo = fopen(filename, "w");
  for(i=0;i<n_points_x;i++){
    for(j=0;j<n_points_y;j++){
      fprintf(archivo, "%f\n", matrix[i+(n_points_y*j)]);
    }
  }
  fclose(archivo);
}

void *copy(float *origen, float *copia, int n_points){
  int i,j;
  for(i=0;i<n_points;i++){
    for(j=0;j<n_points;j++){
      copia[i+(n_points*j)] = origen[i+(n_points*j)];
    }
  }
}

void *asignar(float *matrix, int n_points, float suma){
  int i,j;
  for(i=0;i<n_points;i++){
    for(j=0;j<n_points;j++){
      matrix[i+(n_points*j)] = suma;
      suma = suma + 1.0;
    }
  }
}

void *paso(float *u, float *v, float *un, float *vn, int n_points_x, int n_points_y){
  //matrix un is the copy of u
  
  int i,j, n_points,t;
  float nx,ny,nt,c,dx,dy,sigma,nu,dt;
  for(t=1;t<501;t++){
    nx = 41.0;
    ny = 41.0;
    nt = 120.0;
    c = 1.0;
    dx = 2.0/(nx-1.0);
    dy = 2.0/(ny-1.0);
    sigma = .0009;
    nu = 0.01;
    dt = sigma*dx*dy/nu;
    n_points = n_points_x;

    copy(u, un, n_points_x);
    copy(v, vn, n_points_x);

    for(j=1;j<n_points-1;j++){
      for(i=1;i<n_points-1;i++){
	u[j+(n_points*i)] = un[j+(n_points*i)] - dt/dx*un[j+(n_points*i)]*(un[j+(n_points*i)]-un[(j-1)+(n_points*i)])-dt/dy*vn[j+(n_points*i)]*(un[j+(n_points*i)]-un[j+(n_points*(i-1))])+nu*dt/(dx*dx)*((un[(j+1)+(n_points*i)]-2*un[j+(n_points*i)])+un[(j-1)+(n_points*i)])+nu*dt/(dy*dy)*(un[j+((i+1)*n_points)]-2*un[j+(n_points*i)]+un[j+((i+1)*n_points)]);

	v[j+(n_points*i)] = vn[j+(n_points*i)] - dt/dx*vn[j+(n_points*i)]*(vn[j+(n_points*i)]-vn[(j-1)+(n_points*i)])-dt/dy*un[j+(n_points*i)]*(vn[j+(n_points*i)]-vn[j+(n_points*(i-1))])+nu*dt/(dx*dx)*((vn[(j+1)+(n_points*i)]-2*vn[j+(n_points*i)])+vn[(j-1)+(n_points*i)])+nu*dt/(dy*dy)*(vn[j+((i+1)*n_points)]-2*vn[j+(n_points*i)]+vn[j+((i+1)*n_points)]);
      }
  
      u[0+(j*n_points)]= 1.0;
      u[(n_points-1)+(j*n_points)]= 1.0;
      u[i+0*n_points]= 1.0;
      u[i+((n_points-1)*n_points)]= 1.0;

      v[0+(j*n_points)]= 1.0;
      v[(n_points-1)+(j*n_points)]= 1.0;
      v[i+0*n_points]= 1.0;
      v[i+((n_points-1)*n_points)]= 1.0;
    }
    escribirArchivo(u,n_points_x,n_points_y,t);
  }
}					  
