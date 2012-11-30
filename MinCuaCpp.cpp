
#include "mex.h"
#include <stdio.h>
#include <cmath>

/*Esta funcion obtiene la inversa de una matriz de 3x3.
 *gA: Inversa de la matriz A.
 *A: Matriz a invertir de 3x3.
 *Orlando Medina
 *26 Nov 2012.
 */
void invA(double *gA,double *A){
    
    double detA = A[0]*(A[4]*A[8]-A[7]*A[5]) - A[3]*(A[1]*A[8]-A[7]*A[2]) + A[6]*(A[1]*A[5]-A[4]*A[2]);
//     if(detA <= 0.001){
//         printf("Warning: Singular Matrix.\n");
//         return;
//     }
    
    gA[0] =  (A[4]*A[8]-A[7]*A[5])/detA;    
    gA[1] = -(A[1]*A[8]-A[7]*A[2])/detA;
    gA[2] =  (A[1]*A[5]-A[4]*A[2])/detA;
    
    gA[3] = -(A[3]*A[8]-A[6]*A[5])/detA;    
    gA[4] =  (A[0]*A[8]-A[6]*A[2])/detA;
    gA[5] = -(A[0]*A[5]-A[3]*A[2])/detA;
    
    gA[6] =  (A[3]*A[7]-A[6]*A[4])/detA;
    gA[7] = -(A[0]*A[7]-A[6]*A[1])/detA;
    gA[8] =  (A[0]*A[4]-A[3]*A[1])/detA;
}

/*Esta funcion crea la matriz A para el metodo de Minimos Cuadrados del articulo AIA.
 *gA: matriz A de salida.
 *k: numero de pasos.
 *S: seno de los pasos = sin(steps).
 *C: coseno de los pasos = cos(steps).
 *Autor: Orlando Medina.
 *Fecha: 26 Nov 2012.
 */
void make_A(double *gA,int k,double *S,double *C){
    
    double sS=0 ,sC=0;
    double sSS=0,sCC=0;
    double sSC=0;
    
    for(int n=0;n<k;n++){
        sS += S[n];
        sC += C[n];
        sSS += S[n]*S[n];
        sCC += C[n]*C[n];
        sSC += S[n]*C[n];
    }
    
    gA[0] = k;      gA[3] = sC;     gA[6] = sS;
    gA[1] = sC;     gA[4] = sCC;    gA[7] = sSC;
    gA[2] = sS;     gA[5] = sSC;    gA[8] = sSS;
}

/*Esta funcion crea vector b para el metodo de Minimos Cuadrados del articulo AIA.
 *b: Vector b de salida.
 *I: Interferogramas M*N*k.
 *S: seno de los pasos = sin(steps).
 *C: coseno de los pasos = cos(steps).
 *xy: punto sobre el que se obtiene el vector i
 *Autor: Orlando Medina.
 *Fecha: 26 Nov 2012.
 */
void make_b(double *b,double *I,double *S,double *C,int mI,int nI,int kI,int xy){
    double sI =0;
    double sIS=0;
    double sIC=0;
    
    for(int n=0; n<kI ;n++){
        sI  += I[mI*nI*n+xy];
        sIS += I[mI*nI*n+xy] * S[n];
        sIC += I[mI*nI*n+xy] * C[n];
    }
    b[0]=sI; b[1]=sIC; b[2]=sIS;
}

/* Esta funcion multiplica la Matriz A por el vector b para el metodo de 
 * Minimos Cuadrados del articulo AIA.
 *Ab: Vector Ab de salida.
 *A: Matriz A de 3x3.
 *b: Vector b de 3x1.
 *Autor: Orlando Medina.
 *Fecha: 26 Nov 2012.
 */
void mult_iAb(double *Ab, double *A, double *b){
    
    Ab[0] = A[0]*b[0] + A[3]*b[1] + A[6]*b[2]; 
    Ab[1] = A[1]*b[0] + A[4]*b[1] + A[7]*b[2];
    Ab[2] = A[2]*b[0] + A[5]*b[1] + A[8]*b[2];
}
/* Esta funcion se encarga de OBTENER FASE por Minimos Cuadrados como 
 * el articulo AIA lo propone.
 *DC: Termino de DC calculado.
 *Phi: Puntero a la parte real del campo complejo calculado (fase).
 *Psi:Puntero a la parte imaginaria del campo complejo calculado (fase).
 *I: Interferogramas a demodular.
 *S: seno de los pasos = sin(steps).
 *C: coseno de los pasos = cos(steps).
 *mI: Numero de renglones de I.
 *nI: Numero de columnas de I.
 *kI: Numero de interferogramas (pasos).
 *Autor: Orlando Medina.
 *Fecha: 26 Nov 2012.
 */
void MinimosCuadrados(double *DC,double *Phi,double *Psi,double *I,double *S,double *C,int mI,int nI,int kI){
    
    int xy = 0;
    double A[9];
    double iA[9];
    double b[3];
    double X[3];
    
    make_A(A,kI,S,C);
    invA(iA,A);
    
    for(int r=0;r<mI;r++){
        for(int c=0;c<nI;c++){
            xy = c*mI+r;
            make_b(b,I,S,C,mI,nI,kI,xy);
            mult_iAb(X,iA,b);
            DC[xy] =X[0];
            Phi[xy]=X[1];
            Psi[xy]=X[2];
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  double* I = mxGetPr(prhs[0]);
  double* S = mxGetPr(prhs[1]);
  double* C = mxGetPr(prhs[2]);
    
  const int *Dims = mxGetDimensions(prhs[0]);
  int mI          = Dims[0];
  int nI          = Dims[1];
  int kI          = Dims[2];
  
  int Dims2 [] = {mI,nI};
  
  plhs[0]     = mxCreateNumericArray(2,Dims2,mxDOUBLE_CLASS,mxREAL);
  double*DC   = mxGetPr(plhs[0]);
  plhs[1]     = mxCreateNumericArray(2,Dims2,mxDOUBLE_CLASS,mxCOMPLEX);
  double*Phi  = mxGetPr(plhs[1]);
  double*Psi  = mxGetPi(plhs[1]);
  
  MinimosCuadrados(DC,Phi,Psi,I,S,C,mI,nI,kI);
}