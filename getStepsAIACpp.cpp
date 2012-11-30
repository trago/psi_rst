#include "mex.h"
#include <stdio.h>
//#include <cmath>

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

/*Esta funcion crea la matriz A de coeficientes para obtener los pasos como
 *lo muestra el algoritmo AIA.
 *gA: matriz A de salida.
 *Phi: Puntero a la parte real del campo complejo calculado (fase).
 *Psi:Puntero a la parte imaginaria del campo complejo calculado (fase).
 *mI: Numero de renglones de I.
 *nI: Numero de columnas de I.
 *Autor: Orlando Medina.
 *Fecha: 27 Nov 2012.
 */
void make_A(double *gA,double *Phi,double *Psi,int mI,int nI){
    
    double sSf=0 ,sCf=0;
    double sSSf=0,sCCf=0;
    double sSCf=0;
    
    for(int n=0;n<mI*nI;n++){
        sSf += Psi[n];
        sCf += Phi[n];
        sSSf += Psi[n]*Psi[n];
        sCCf += Phi[n]*Phi[n];
        sSCf += Psi[n]*Phi[n];
    }
    
    gA[0] = mI*nI;   gA[3] = sCf;     gA[6] = sSf;
    gA[1] = sCf;     gA[4] = sCCf;    gA[7] = sSCf;
    gA[2] = sSf;     gA[5] = sSCf;    gA[8] = sSSf;
}

/*Esta funcion crea vector b para obtener los pasos como lo muestra el
 *algoritmo AIA.
 *b: Vector b de salida.
 *I: Interferogramas de M*N*k.
 *Phi: Puntero a la parte real del campo complejo calculado (fase).
 *Psi:Puntero a la parte imaginaria del campo complejo calculado (fase).
 *mI: Numero de renglones de I.
 *nI: Numero de columnas de I.
 *k: Numero de interferograma a tratar.
 *Autor: Orlando Medina.
 *Fecha: 27 Nov 2012.
 */
void make_b(double *b,double *I,double *Phi,double *Psi,int mI,int nI,int k){
    double sI =0;
    double sIS=0;
    double sIC=0;
    int    xyk=0;
    
    for(int n=0; n<mI*nI ;n++){
        xyk = mI*nI*k+n;
        sI  += I[xyk];
        sIS += I[xyk] * Psi[n];
        sIC += I[xyk] * Phi[n];
    }
    b[0]=sI; b[1]=sIC; b[2]=sIS;
}

/* Esta funcion multiplica la Matriz A por el vector b para el metodo de 
 * Minimos Cuadrados del articulo AIA.
 *Ab: Vector Ab de salida 3x1.
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
 *Fecha: 27 Nov 2012.
 */
void getStepsAIA(double *S,double *C,double *I,double *Phi,double *Psi,int mI,int nI,int kI){
    
    double A[9];
    double iA[9];
    double b[3];
    double X[3];
    
    make_A(A,Phi,Psi,mI,nI);
    invA(iA,A);
    
    for(int k=0;k<kI;k++){
            make_b(b,I,Phi,Psi,mI,nI,k);
            mult_iAb(X,iA,b);
            //DC[k] = X[0];
            C [k] = X[1];
            S [k] = X[2];
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  double* I   = mxGetPr(prhs[0]);
  double* Phi = mxGetPr(prhs[1]);
  double* Psi = mxGetPi(prhs[1]);
    
  const int *Dims = mxGetDimensions(prhs[0]);
  int mI          = Dims[0];
  int nI          = Dims[1];
  int kI          = Dims[2];
  
  int Dims2 [] = {mI,nI};
  int Dims3 [] = {1,kI};
  
  //plhs[0]     = mxCreateNumericArray(2,Dims3,mxDOUBLE_CLASS,mxREAL);
  //double*DC   = mxGetPr(plhs[0]);
  plhs[0]     = mxCreateNumericArray(2,Dims3,mxDOUBLE_CLASS,mxREAL);
  double*S    = mxGetPr(plhs[0]);
  plhs[1]     = mxCreateNumericArray(2,Dims3,mxDOUBLE_CLASS,mxREAL);
  double*C    = mxGetPr(plhs[1]);
  
  getStepsAIA(S,C,I,Phi,Psi,mI,nI,kI);
}