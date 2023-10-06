
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <math.h>
#include <complex.h>
#include <algorithm>
////testing the fortran interface
extern "C" void Eprim(int nsource,int ntarget, double *rs,double* js,double* robs, double* Eprimary, double iprec);
extern "C" void Hprim(int nsource,int ntarget, double *rs,double* js,double* robs, double* Hprimary, double iprec);
extern "C" void Aprim(int nsource,int ntarget, double *rs,double* js,double* robs, double* Aprim,double*curlAprim, double iprec);
extern "C" void Ephiprim(int nsource,int ntarget, double *rs,double* rho,double* robs, double* Ephiprimary,double iprec);
extern "C" void Eprimrsjsksrho(int nsource,int ntarget, double *rs,double* js,double* ks,double* rho,double* robs, double* Eprimary,double iprec);
extern "C" void  computeeprimary_(double *rs,double *js,double* robs,double* Eprimary,int *ntarget,int *npoints,double *iprec);
extern "C" void  computehprimary_(double *rs,double *js,double* robs,double* Eprimary,int *ntarget,int *npoints,double *iprec);
extern "C" void  computeaprimary_(double *rs,double *js,double* robs,double* Aprimary,double* curlAprimary,int *ntarget,int *npoints,double *iprec);
extern "C" void  computeephiprimary_(double *rs,double *rho,double* robs,double* Ephiprimary,int *ntarget,int *npoints,double *iprec);


void Eprim(int nsource,int ntarget, double *rs,double* js,double* robs, double* Eprimary,double iprec) {


computeeprimary_(rs,js,robs,Eprimary,&ntarget,&nsource,&iprec);
}

void Eprimrsjsksrho(int nsource,int ntarget, double *rs,double* js,double* ks,double* rho,double* robs, double* Eprimary,double iprec) {
const double mu0=1.256637061435917e-06;
                     double *Eprim2 = (double *)malloc(sizeof(double)*ntarget*3);//unique face 2 tetra array
int j;
computeeprimary_(rs,js,robs,Eprimary,&ntarget,&nsource,&iprec);
computehprimary_(rs,ks,robs,Eprim2,&ntarget,&nsource,&iprec);
		for (j=0; j<ntarget*3; j++) {			 
					 Eprimary[j]=Eprimary[j]-Eprim2[j]/mu0;
		}
computeephiprimary_(rs,rho,robs,Eprim2,&ntarget,&nsource,&iprec);
		for (j=0; j<ntarget*3; j++) {			 
					 Eprimary[j]=Eprimary[j]+Eprim2[j];
		}

}

void Hprim(int nsource,int ntarget, double *rs,double* js,double* robs, double* Hprimary,double iprec) {


computehprimary_(rs,js,robs,Hprimary,&ntarget,&nsource,&iprec);
}

void Aprim(int nsource,int ntarget, double *rs,double* js,double* robs, double* Aprimary,double* curlAprimary,double iprec) {


computeaprimary_(rs,js,robs,Aprimary,curlAprimary,&ntarget,&nsource,&iprec);
}


void Ephiprim(int nsource,int ntarget, double *rs,double* rho,double* robs, double* Ephiprimary,double iprec) {


computeephiprimary_(rs,rho,robs,Ephiprimary,&ntarget,&nsource,&iprec);
}

