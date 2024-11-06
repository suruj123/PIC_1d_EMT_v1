//#include "mpi.h"
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


void FS(){

int i;


//Calculate By and Ex


 for(i=1;i<=Nx-1;i++){

  Bzg[i]  += - 0.5*(deltaT/del_x)*(ef_yg[i+1] - ef_yg[i]);

 }

  Bzg[Nx] += - 0.5*(deltaT/del_x)*(ef_yg[1] - ef_yg[Nx]);   //Implementation of PBC on Magnetic Field





 for(i=2;i<=Nx;i++){

  Byg[i]  +=  0.5*(deltaT/del_x)*(ef_zg[i] - ef_zg[i-1]);

 }

  Byg[1] +=  0.5*(deltaT/del_x)*(ef_zg[1] - ef_zg[Nx]);   //Implementation of PBC on Magnetic Field




 //Byg is calculated at the half steps

 for(i=1;i<=Nx-1;i++){

  Byg[i] = 0.5*( Byg[i+1] + Byg[i] );

 }

  Byg[Nx] = 0.5*( Byg[1] + Byg[Nx] );


/*

 for(i=2;i<=Nx;i++){

 ef_yg[i] += -(deltaT)*( ( Bzg[i] - Bzg[i-1] )/del_x + jgy[i] );

 }

 ef_yg[1] += -(deltaT)*( ( Bzg[1] - Bzg[Nx] )/del_x + jgy[1] );  //Implementation of PBC on Electric Field

*/





}









