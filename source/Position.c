#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"global.h"
#include<time.h>
#include<string.h>


void Position(){


   int n;




   //Position is updated

   for(n=1;n<=nAtom;n++){

   //vx[n]  += Q[n] * 0.5 * deltaT * ax[n];
 
   rx[n]  += deltaT * vx[n]*0.5;

   }



   //Periodic Boundary condition

   for(n=1;n<=nAtom;n++){

   rx[n] -= region[1]*rint(rx[n]/region[1]);
   //ry[n] -= region[2]*rint(ry[n]/region[2]);

   }




}
