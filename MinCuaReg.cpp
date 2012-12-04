
#include "mex.h"
//#include <cmath>

void MinCuaReg(double *gDC,double *gPhi,double *gPsi,double *DC,double *Phi,double *Psi,
               double *SSalpha,double *CCalpha,double *I,double lambdaf,double lambdaDC,int mI,int nI,int kI){
    
    double sCk=0,sSk=0;
    double sCCk=0,sSSk=0;
    double sSCk=0,sIk=0;
    double sIS=0,sIC=0;
    int    xy=0,xyk =0;
    double reg_Phi = 0;
    double reg_Psi = 0;
    double reg_DC  = 0;
    
    for(int j=0; j < nI;j++){ // Recorre las Columnas.
        for(int i=0; i < mI;i++){ // Recorre los Renglones.
            sCk =0;sSk =0;
            sCCk=0;sSSk=0;
            sSCk=0;sIk =0;
            sIS =0;sIC =0;
            xy = j*mI+i;
            for(int k=0; k < kI; k++){ // Recorre los k interferogramas.
                xyk = (mI*nI*k)+xy;
                
                sCk     += CCalpha[k];
                sSk     += SSalpha[k];
                sCCk    += CCalpha[k]*CCalpha[k];
                sSSk    += SSalpha[k]*SSalpha[k];
                sSCk    += CCalpha[k]*SSalpha[k];
                
                sIk     += I[xyk];
                sIS     += I[xyk]*SSalpha[k];
                sIC     += I[xyk]*CCalpha[k];
            }
            reg_Phi = 0;
            reg_Psi = 0;
            reg_DC  = 0;
            int N  = 0;
            if(i-1 >= 0){
                reg_Phi += Phi[xy-1];
                reg_Psi += Psi[xy-1];
                reg_DC  += DC [xy-1];
                N+=1;
            }
            if(i+1 < mI){
                reg_Phi += Phi[xy+1];
                reg_Psi += Psi[xy+1];
                reg_DC  += DC [xy+1];
                N+=1;
            }
            if(j-1 >= 0){
                reg_Phi += Phi[xy-mI];
                reg_Psi += Psi[xy-mI];
                reg_DC  += DC [xy-mI];
                N+=1;
            }
            if(j+1 < nI){
                reg_Phi += Phi[xy+mI];
                reg_Psi += Psi[xy+mI];
                reg_DC  += DC [xy+mI];
                N+=1;
            }
            gDC[xy]  = (-sCk*Phi[xy] + sSk * Psi[xy] + sIk + lambdaDC*reg_DC )/(N*lambdaDC+kI);
            gPhi[xy] = (-sCk*gDC[xy] + sSCk* Psi[xy] + sIC + lambdaf *reg_Phi)/(N*lambdaf +sCCk);
            gPsi[xy] = ( sSk*gDC[xy] + sSCk*gPhi[xy] - sIS + lambdaf *reg_Psi)/(N*lambdaf +sSSk);
            
        }
    }
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  if (nrhs!=7) {
    mexErrMsgTxt("Se requieren 7 parametros");
  }
  double *I       = mxGetPr(prhs[0]);
  double *Phi     = mxGetPr(prhs[1]);
  double *Psi     = mxGetPi(prhs[1]);
  double *DC      = mxGetPr(prhs[2]);
  double *SSalpha = mxGetPr(prhs[3]);
  double *CCalpha = mxGetPr(prhs[4]);
  double lambdaf  = *mxGetPr(prhs[5]);
  double lambdaDC = *mxGetPr(prhs[6]);
  
  const int *Dims = mxGetDimensions(prhs[0]);

  int mI          = Dims[0];
  int nI          = Dims[1];
  int kI          = Dims[2];
  int Dims2[2]    = {mI,nI};
  
  plhs[0]         = mxCreateNumericArray(2,Dims2,mxDOUBLE_CLASS,mxCOMPLEX);
  double* gPhi    = mxGetPr(plhs[0]);
  double* gPsi    = mxGetPi(plhs[0]);
  plhs[1]         = mxCreateNumericArray(2,Dims2,mxDOUBLE_CLASS,mxREAL);
  double* gDC     = mxGetPr(plhs[1]);
  
  MinCuaReg(gDC,gPhi,gPsi,DC,Phi,Psi,SSalpha,CCalpha,I,lambdaf,lambdaDC,mI,nI,kI);
}