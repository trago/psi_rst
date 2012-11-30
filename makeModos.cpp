/**
 * @file makemodos.cpp
 * @author Orlando M Medina Cazares
 * @date 09/11/2011 - 13:16:00
 * @version 2.0
 * 
 * @brief Este programa implementa una secuencia de modos de vibracion para
 *          para una placa. La version anterior (1.0) tiene un error a la
 *          indexacion, en esta version se corrigio este problema.
 */

#include "mex.h"
#include <stdio.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.1416
#endif
#define PI M_PI

void makeModos(double *Iout, int Rows, int Cols,int frames, int modosX, int modosY, double wt = PI/2 ){

    double *frame = NULL;
    
    for(int k=0; k<frames; k++){
        frame = Iout + Rows*Cols*k ;
        for(int i=0; i < Rows;i++)
    			for(int j=0; j < Cols;j++)
        		frame[i+j*Rows] = sin((PI*modosX*j)/Cols) * sin((PI*modosY*i)/Rows) * cos(wt*k);
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  int mI        = *mxGetPr(prhs[0]);
  int nI        = *mxGetPr(prhs[1]);
  double frames  = *mxGetPr(prhs[2]);
  int modosX    = *mxGetPr(prhs[3]);
  int modosY    = *mxGetPr(prhs[4]);
  double wt     = *mxGetPr(prhs[5]);
  
  int Dims [] = {mI,nI,frames};
  plhs[0] = mxCreateNumericArray(3,Dims,mxDOUBLE_CLASS,mxREAL);
  double* Iout = mxGetPr(plhs[0]);
  makeModos(Iout,mI,nI,frames,modosX,modosY,wt);
}

