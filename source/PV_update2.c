//#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"global.h"
#include<time.h>
#include<string.h>


void PV_update2(){

int i, n;

double t_dash, s_dash;


 /*for(n=1;n<=nAtom;n++){

   ax[n] = -Q[n]*ef_x[n];
   ay[n] = -Q[n]*ef_y[n];
 }*/


  //Calculation of v_minus

  for(n=1;n<=nAtom;n++){

   ax[n] = -Q[n]*ef_x[n];

   ay[n] = -Q[n]*ef_y[n];

   az[n] = -Q[n]*ef_z[n];


   vx[n]  += Q[n] * 0.5 * deltaT * ax[n];

   vy[n]  += Q[n] * 0.5 * deltaT * ay[n];
   
   vz[n]  += Q[n] * 0.5 * deltaT * az[n];

  }

   //Rotation

   //t_dash = 0.5*wc*deltaT;

   //s_dash = 2*(t_dash)/(1 + (t_dash)*(t_dash));

  
   
   for(n=1;n<=nAtom;n++){


    //t_dash = 0.5*(Bzn[n] + B0)*deltaT;

    vx_dash[n] = vx[n] + vy[n]*(Bzn[n] + B0 )- vz[n]*Byn[n]; //Calculation of V dash using V minus

    vy_dash[n] = vy[n] - vx[n]*(Bzn[n] + B0);

    vz_dash[n] = vz[n] + vx[n]*Byn[n];


   }   
  
/*
   for(n=1;n<=nAtom;n++){                    

    t_dash = 0.5*(Bzn[n] + B0)*deltaT;

    s_dash = 2*(t_dash)/(1 + (t_dash)*(t_dash));

    vx[n] +=  (vy_dash[n])*s_dash; //Calculation of V plus form V minus, V dash and S dash

    vy[n] +=  (-vx_dash[n])*s_dash;

   }
*/


  double boris, cons = 1.0;

  for(n=1;n<=nAtom;n++){

    boris = 2.0/( (cons + Sqr(bx)) + Sqr(Byn[n]) + Sqr(Bzn[n] + B0) );

    vx[n] +=   boris*(vy_dash[n]*(Bzn[n] + B0) - vz_dash[n]*Byn[n]); 

    vy[n] +=  boris*( vz_dash[n]*(bx)  -   vx_dash[n]*(Bzn[n] + B0) );

    vz[n] +=   boris*(vx_dash[n]*Byn[n] -  vy_dash[n]*bx);

  }


   //Second Half Acceleration

   for(n=1;n<=nAtom;n++){

    vx[n]  += Q[n] * 0.5 * deltaT * ax[n]; // This calculates Vx (n + 1/2)

    vy[n]  += Q[n] * 0.5 * deltaT * ay[n]; // This calculates Vy (n + 1/2)

    vz[n]  += Q[n] * 0.5 * deltaT * az[n]; // This calculates Vy (n + 1/2)

   }

/*
   for(i=1;i<=Nx-1;i++){

    Bzg[i]  += - 0.5*(deltaT/del_x)*(ef_yg[i+1] - ef_yg[i]);

   }

    Bzg[Nx] += - 0.5*(deltaT/del_x)*(ef_yg[1] - ef_yg[Nx]);   //Implementation of PBC on Magnetic Field

 */


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



/*   
   for(n=1;n<=nAtom;n++){

   vx[n] += cos(wc * deltaT) * vx[n] + sin(wc * deltaT) * vy[n];

   vy[n] += -sin(wc * deltaT) * vx[n] + cos(wc * deltaT) * vy[n];

   }
*/


   double  vv = 0;

   for(n=1;n<=nAtom;n++){

   vx_dummy[n] = 0.0;
   vy_dummy[n] = 0.0;

   }


   for(n=1;n<=nAtom;n++){

   vx_dummy[n] = vx[n] + Q[n] * 0.5 * deltaT * ax[n];

   }



   for(n=1;n<=nAtom;n++){

   vx_dummy[n] += cos(wc * deltaT) * vx[n] + sin(wc * deltaT) * vy[n];

   vy_dummy[n] += -sin(wc * deltaT) * vx[n] + cos(wc * deltaT) * vy[n];

   }


   for(n=1;n<=nAtom;n++){

   vv += vx_dummy[n]*vx_dummy[n] + vy_dummy[n]*vy_dummy[n];

   }


   KE = 0.5*vv;

   KE = KE/nAtom;

   PE = 0;

   for(n=1;n<=nAtom;n++){

   PE += (ef_x[n])*(ef_x[n]);// + (ef_y[n])*(ef_y[n]);

   }

   PE = (0.5*PE)/nAtom;

   //printf("pe = %lf %lf\n", PE, timeNow);

/*
   for(n=1;n<=nAtom;n++){

   rx[n] -= region[1]*rint(rx[n]/region[1]);
   //ry[n] -= region[2]*rint(ry[n]/region[2]);

   }
*/

}
