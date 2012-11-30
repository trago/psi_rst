#include <mex.h>
#include <stdio.h>


void makeRampa(double *Iout,double wx,double wy,int Rows,int Cols){
    for(int i=0; i < Rows;i++)
        for(int j=0; j < Cols;j++)
            Iout[j*Rows+i] = wx*i + wy*j;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  if (nrhs!=4) {
    mexErrMsgTxt("Se requieren 4 parametros.");
  }
  double wx = *mxGetPr(prhs[0]);
  double wy = *mxGetPr(prhs[1]);
  int    mI = *mxGetPr(prhs[2]);
  int    nI = *mxGetPr(prhs[3]);

  
  plhs[0]         = mxCreateDoubleMatrix(mI,nI,mxREAL);
  double* gRampa  = mxGetPr(plhs[0]);
  
  makeRampa(gRampa,wx,wy,mI,nI);
 }