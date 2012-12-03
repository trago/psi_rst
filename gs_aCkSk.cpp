
#include "mex.h"
//#include <stdio.h>
#include <cmath>

void getDCxy(double *gDC,double *DC,double *I,double *Phi,double *Psi,double *SSalpha,
             double *CCalpha,double lambda,int mI,int nI,int kI,int i,int j);

void getSkCkxy(double *gSSalpha,double *gCCalpha,double *I,double *Phi,double *Psi,
               double *DC,double *SSalpha,double *CCalpha,double lambda,int mI,int nI,int kI,int i,int j);


void gs_DCCkSk_RST(double *gDC,double *gSSalpha,double *gCCalpha,double *I,double *Phi,double *Psi,
                   double *DC,double *SSalpha,double *CCalpha,double lambdaDC,double lambdaSC,int mI,int nI,int kI)
{
    //int w=1; // limite de la sumatoria para la ventana para el punto xy.
       
    for(int j=0; j < nI;j++){ // Recorre las Columnas donde entra la ventana
        for(int i=0; i < mI;i++){ // Recorre los Renglones donde entra la ventana
            
//             getDCxy(gDC,DC,I,Phi,Psi,SSalpha,CCalpha,lambda,mI,nI,kI,i,j); // encuentra el termino de DC en el punto xy=j*mI+i;
//             getSkCkxy(gSSalpha,gCCalpha,I,Phi,Psi,gDC,SSalpha,CCalpha,mI,nI,kI,i,j);// encuentra los k-esimos cosenos y senos en el punto xy=j*mI+i;
            
            
            getSkCkxy(gSSalpha,gCCalpha,I,Phi,Psi,DC,SSalpha,CCalpha,lambdaSC,mI,nI,kI,i,j);// encuentra los k-esimos cosenos y senos en el punto xy=j*mI+i;
            getDCxy(gDC,DC,I,Phi,Psi,gSSalpha,gCCalpha,lambdaDC,mI,nI,kI,i,j); // encuentra el termino de DC en el punto xy=j*mI+i;
            
            
        }
    }
}

void getSkCkxy(double *gSSalpha,double *gCCalpha,double *I,double *Phi,double *Psi,
               double *DC,double *SSalpha,double *CCalpha,double lambda,int mI,int nI,int kI,int i,int j){
    
    double sPhi=0,sPsi=0;
    double sPhi2=0,sPsi2=0;
    double sPhiPsi=0;
    double sIkPhi=0,sIkPsi=0;
    double reg_Ck=0,reg_Sk=0;
    
    int xyk=0;
    int m_xy=0,m_xyk=0;
    int xy = j*mI+i;
    
    int mi=-1,mf=1;
    int ni=-1,nf=1;
    
    // Verifica en que borde se encuentra la ventana
    if(i ==    0){mi=0 ;mf=1;}
    if(i == mI-1){mi=-1;mf=0;}
    if(j ==    0){ni=0 ;nf=1;}
    if(j == nI-1){ni=-1;nf=0;}

    
    for(int k=0;k < kI;k++){ // Recorre los kI Interferogramas
        xyk = (mI*nI*k)+xy;

        sIkPhi = 0;
        sIkPsi = 0;

        for(int m=mi;m <= mf;m++){ // Desplaza los Renglones
            for(int n=ni;n <= nf;n++){ // Desplaza las Columnas
                m_xy  = (j+n)*mI+(i+m);
                m_xyk = (mI*nI*k)+m_xy;

                // Esta condicion es para calcular solo una vez las sumatorias sobre Phi y Psi.
                if(k==0){ 
                    sPhi    += Phi[m_xy];
                    sPsi    += Psi[m_xy];
                    sPhi2   += Phi[m_xy]*Phi[m_xy];
                    sPsi2   += Psi[m_xy]*Psi[m_xy];
                    sPhiPsi += Phi[m_xy]*Psi[m_xy];
                }
                sIkPhi  += (I[m_xyk]*Phi[m_xy]);
                sIkPsi  += (I[m_xyk]*Psi[m_xy]);
            }
        }
        
      // Seccion de regularizacion.
      reg_Ck = 0;
      reg_Sk = 0;
      int N  = 0;
      if(i-1 >= 0){
        reg_Ck += CCalpha[xyk-1];
        reg_Sk += SSalpha[xyk-1];
        N+=1;
      }
      if(i+1 < mI){
        reg_Ck += CCalpha[xyk+1];
        reg_Sk += SSalpha[xyk+1];
        N+=1;
      }
      if(j-1 >= 0){
        reg_Ck += CCalpha[xyk-mI];
        reg_Sk += SSalpha[xyk-mI];
        N+=1;
      }
      if(j+1 < nI){
        reg_Ck += CCalpha[xyk+mI];
        reg_Sk += SSalpha[xyk+mI];
        N+=1;
      }
        
        gCCalpha[xyk] = (-DC[xy]*sPhi + sPhiPsi* SSalpha[xyk] + sIkPhi + lambda*reg_Ck) / (sPhi2 + N*lambda);
        gSSalpha[xyk] = ( DC[xy]*sPsi + sPhiPsi*gCCalpha[xyk] - sIkPsi + lambda*reg_Sk) / (sPsi2 + N*lambda);
    }
}


void getDCxy(double *gDC,double *DC,double *I,double *Phi,double *Psi,double *SSalpha,
             double *CCalpha,double lambda,int mI,int nI,int kI,int i,int j){
    
    int xyk=0;
    int xy = j*mI+i;
    int elementos=9; // Numero de elementos que tendra la ventana.
    int m_xy=0,m_xyk=0;
    double sPhi=0,sPsi=0;
    double sI=0;
    double sCC=0,sSS=0;
    double reg = 0;
    int N  = 0;
    
    int mi=-1,mf=1;
    int ni=-1,nf=1;
    
    // Verifica en que borde se encuentra la ventana
    if(i ==    0){mi=0 ;mf=1;}//elementos-=3;}
    if(i == mI-1){mi=-1;mf=0;}//elementos-=3;}
    if(j ==    0){ni=0 ;nf=1;}//elementos-=3;}
    if(j == nI-1){ni=-1;nf=0;}//elementos-=3;}
    
//     if(i==0    & j==0){    elementos=3;}
//     if(i==mI-1 & j==nI-1){ elementos=3;}
//     if(i==0    & j==nI-1){ elementos=3;}
//     if(i==mI-1 & j==0){    elementos=3;}
    
    for(int k=0;k < kI;k++){ // Recorre los kI Interferogramas
        xyk = (mI*nI*k)+xy;

        for(int m=mi;m <= mf;m++){ // Desplaza los Renglones
            for(int n=ni;n <= nf;n++){ // Desplaza las Columnas
                m_xy  = (j+n)*mI+(i+m);
                m_xyk = (mI*nI*k)+m_xy;

                // Esta condicion es para calcular solo una vez las sumatorias sobre Phi y Psi.
                if(k==0){ 
                    sPhi += Phi[m_xy];
                    sPsi += Psi[m_xy];
                }
                sI += I[m_xyk];
            }
        }
        sCC += CCalpha[xyk];
        sSS += SSalpha[xyk];
    }
    
      reg = 0;
      N  = 0;
      if(i-1 >= 0){
        reg += DC[j*mI+(i-1)];
        N+=1;
      }
      if(i+1 < mI){
        reg += DC[j*mI+(i+1)];
        N+=1;
      }
      if(j-1 >= 0){
        reg += DC[(j-1)*mI+i];
        N+=1;
      }
      if(j+1 < nI){
        reg += DC[(j+1)*mI+i];
        N+=1;
      }
    
    gDC[xy] = (-sPhi*sCC + sPsi*sSS + sI + lambda*reg)/(kI*elementos+N*lambda);

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  if (nrhs!=7) {
    mexErrMsgTxt("Se requieren 7 parametros (el segundo se espera que sea Complejo)");
  }
  double *I       = mxGetPr(prhs[0]);
  double *Phi     = mxGetPr(prhs[1]);
  double *Psi     = mxGetPi(prhs[1]);
  double *DC      = mxGetPr(prhs[2]);
  double *SSalpha = mxGetPr(prhs[3]);
  double *CCalpha = mxGetPr(prhs[4]);
  double lambdaDC   = *mxGetPr(prhs[5]);
  double lambdaSC   = *mxGetPr(prhs[6]);
  
  const int *Dims = mxGetDimensions(prhs[0]);

  int mI          = Dims[0];
  int nI          = Dims[1];
  int kI          = Dims[2];
  int Dims2[2]    = {mI,nI};
  
  plhs[0]         = mxCreateNumericArray(2,Dims2,mxDOUBLE_CLASS,mxREAL);
  double* gDC     = mxGetPr(plhs[0]);
  plhs[1]         = mxCreateNumericArray(3,Dims ,mxDOUBLE_CLASS,mxREAL);
  double* gSSalpha= mxGetPr(plhs[1]);
  plhs[2]         = mxCreateNumericArray(3,Dims ,mxDOUBLE_CLASS,mxREAL);
  double* gCCalpha= mxGetPr(plhs[2]);
  
  gs_DCCkSk_RST(gDC,gSSalpha,gCCalpha,I,Phi,Psi,DC,SSalpha,CCalpha,lambdaDC,lambdaSC,mI,nI,kI);
}