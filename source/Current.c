#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"global.h"
#include<time.h>
#include<string.h>


 void Current(){

 int i, j, k, n, m;

 double dx, dy, tx, ty;


 for(n=1;n<=nAtom;n++){

  xs[n] = rx[n] + regionH[1];

  jy[n] = Q[n]*vy[n];  //Calculation of the Current density jy

  jz[n] = Q[n]*vz[n];  //Calculation of the Current density jz

 }

 for(n=1;n<=Nx;n++){
 xg[n] = (n-1)*del_x;
 }


 for(i=1;i<=Nx;i++){
 
   jgy[i] = 0.0;   

   jgz[i] = 0.0;   

 } 


for(n=1;n<=nAtom;n++){

  i = (int)(xs[n]/del_x) + 1; 
  //j = (int)(ys[n]/del_y) + 1; 

  j = i+1;

 
  tx = xg[j] - xs[n];
  ty = xs[n] - xg[i];
 

   //For Current Density jgy

  jgy[i]      =  jgy[i]    +  jy[n]*tx*h_area; //(i, j)  

  jgy[j]      =  jgy[j]    +  jy[n]*ty*h_area; //(i, j)  



  //For Current Density jgz

  jgz[i]      =  jgz[i]    +  jz[n]*tx*h_area; //(i, j)  

  jgz[j]      =  jgz[j]    +  jz[n]*ty*h_area; //(i, j)  

 }


 for(i=2;i<=Nx;i++){
 
  jgy[i] = 0.5*( jgy[i] + jgy[i-1] );

 }
 
  jgy[1] = 0.5*( jgy[1] + jgy[Nx] );

 

 //Subtractig the mean from the current density

 double avjgy = 0.0, avjgz = 0.0;
  
 for(i=1;i<=Nx;i++){
 
  avjgy += jgy[i];

  avjgz += jgz[i];

 }

 avjgy = avjgy/Nx;

 avjgz = avjgz/Nx;

 for(i=1;i<=Nx;i++){
 
  jgy[i] = jgy[i] - avjgy;

  jgz[i] = jgz[i] - avjgz;

 }











 }
