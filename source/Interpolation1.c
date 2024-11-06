#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"global.h"
#include<time.h>
#include<string.h>
#define IMUL    314159269
#define IADD    453806245
#define MASK    2147483647
#define SCALE   0.4656612873e-9

#define MASTER 0

 void Interpolation1(){
 
 int i, j, k, n, m;

 double dx, dy, tx, ty;
 
 //del_x = region[1]/(Nx-1);
 //del_y = region[2]/(Ny-1);
 
 //h_area = 1.0/(del_x*del_y);


 printf("interpolation 1 start\n");

 for(n=1;n<=nAtom;n++){


  xs[n] = rx[n] + regionH[1];


 }

 for(n=1;n<=Nx;n++){
 xg[n] = (n-1)*del_x;
 }

/* for(n=1;n<=Ny;n++){
 yg[n] = (n-1)*del_y;
 }
*/ 
 
for(i=1;i<=Nx;i++){
 
   
   qg[i] = 0.0;   

   //jgy[i] = 0.0;   

   //printf("qg = %lf\n", qg[n]);

 } 


for(n=1;n<=nAtom;n++){

  i = (int)(xs[n]/del_x) + 1; 
  //j = (int)(ys[n]/del_y) + 1; 

  j = i+1;

 
  tx = xg[j] - xs[n];
  ty = xs[n] - xg[i];
 
  //For Charge Density

  qg[i]       =  qg[i]     +  Q[n]*tx*h_area; //(i, j)  

  qg[j]       =  qg[j]     +  Q[n]*ty*h_area; //(i, j)  


  }

  /*

 //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   //For Current Density

  jgy[i]      =  jgy[i]    +  jy[n]*tx*h_area; //(i, j)  

  jgy[j]      =  jgy[j]    +  jy[n]*ty*h_area; //(i, j)  



 for(i=2;i<=Nx;i++){
 
  jgy[i] = 0.5*( jgy[i] + jgy[i-1] );

 }
 
 jgy[1] = 0.5*( jgy[1] + jgy[Nx] );

 

 //Subtractig the mean from the current density

 double avjgy = 0.0;
  
 for(i=1;i<=Nx;i++){
 
  avjgy += jgy[i];

 }

 avjgy = avjgy/Nx;

 for(i=1;i<=Nx;i++){
 
  jgy[i] = jgy[i] - avjgy;

 }

*/

 double sum = 0;

 FILE *fp;
 fp = fopen("rho_info", "w");


 double fac = ((double)Nx)/((double) nAtom);

 for(i=1;i<=Nx;i++){
 
   
   //rho[n] = qg[n]/(del_x * del_y);   
   //rho[n] = qg[n]/(region[1] * region[2]);   
   //rho[i] = qg[i]*(512/40000.0);   
   rho[i] = qg[i]*fac;   

   sum += qg[i];
   //printf("qg = %lf\t rho = %lf %d\n", qg[i], rho[i], Nx);
   fprintf(fp, "%lf\n", rho[i]);

 }

 fclose(fp);

 for(i=1;i<=Nx;i++){
 
   rho1[i-1] = rho[i] -1;   

 }

 printf("sum = %lf\t fac = %lf\t Nx = %d\t nAtom = %d\n", sum, fac, Nx, nAtom);
 
 
}
