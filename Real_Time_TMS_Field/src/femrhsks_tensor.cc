#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <math.h>
#include <complex.h>
#include <algorithm>
extern "C" void  computehprimary_(double *rs,double *js,double* robs,double* Eprimary,int *ntarget,int *npoints,double *iprec);
extern "C" void rhsfunckstensor(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *rs,double* js,double* rhs,double iprec);

// these will be done later
//extern "C" void rhsfunc2ndks_tensor(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *rs,double* js,double* rhs,int iprec);
//extern "C" void rhsfunc3rdks_tensor(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *rs,double* js,double* rhs,int iprec);



void rhsfunckstensor(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *rs,double* js, double* rhs,double iprec) {
  const double mu0=1.25663706e-6;
int i,j,k,iquad;
double Gtemp[12];
double Gtemp2[12];
double jv[3];
double volume;
int te_stride=4;
int ntarget;
double *rv = (double *)malloc(sizeof(double)*nte*3);//unique face 2 tetra array
double *Eprimary = (double *)malloc(sizeof(double)*nte*3);//unique face 2 tetra array



//double *Eprimarypt = (double *)malloc(sizeof(double)*np*3);//unique face 2 tetra array
//for (i=0;i<3*np;i++) {
//Eprimarypt[i]=0.0;
//}
//for (i=0;i<3*np;i++) {
//rv[i]=p[i];
//}
//ntarget=np;
//computehprimary_(rs,js,rv,Eprimarypt,&ntarget,&nsource,&iprec);


//int nquad=5;
//double qwt []={-8.000000000000000e-01,4.500000000000000e-01,4.500000000000000e-01,4.500000000000000e-01,4.500000000000000e-01};
//double qpt []={2.500000000000000e-01,2.500000000000000e-01,2.500000000000000e-01,2.500000000000000e-01,
//               5.000000000000000e-01,1.666666666666667e-01,1.666666666666667e-01,1.666666666666667e-01,
//               1.666666666666667e-01,5.000000000000000e-01,1.666666666666667e-01,1.666666666666667e-01,
//               1.666666666666667e-01,1.666666666666667e-01,5.000000000000000e-01,1.666666666666667e-01,
//               1.666666666666667e-01,1.666666666666667e-01,1.666666666666667e-01,5.000000000000000e-01};
int nquad=1;
double qwt []={1};
double qpt []={2.500000000000000e-01,2.500000000000000e-01,2.500000000000000e-01,2.500000000000000e-01};

  for (iquad=0;iquad<nquad;iquad++) {
for (i=0;i<nte;i++) {
rv[3*i+0]=(p[te2p[i*4+0]*3+0]*qpt[4*iquad]+  p[te2p[i*4+1]*3+0]*qpt[4*iquad+1]
          +p[te2p[i*4+2]*3+0]*qpt[4*iquad+2]+p[te2p[i*4+3]*3+0]*qpt[4*iquad+3]);
rv[3*i+1]=(p[te2p[i*4+0]*3+1]*qpt[4*iquad]+  p[te2p[i*4+1]*3+1]*qpt[4*iquad+1]
          +p[te2p[i*4+2]*3+1]*qpt[4*iquad+2]+p[te2p[i*4+3]*3+1]*qpt[4*iquad+3]);
rv[3*i+2]=(p[te2p[i*4+0]*3+2]*qpt[4*iquad]+  p[te2p[i*4+1]*3+2]*qpt[4*iquad+1]
          +p[te2p[i*4+2]*3+2]*qpt[4*iquad+2]+p[te2p[i*4+3]*3+2]*qpt[4*iquad+3]);
}
ntarget=nte;
for (i=0;i<3*nte;i++) {
Eprimary[i]=0.0;
}
      computehprimary_(rs,js,rv,Eprimary,&ntarget,&nsource,&iprec);


//for (i=0;i<nte;i++) {
//Eprimary[3*i+0]=(Eprimarypt[te2p[i*4+0]*3+0]*qpt[4*iquad]+Eprimarypt[te2p[i*4+1]*3+0]*qpt[4*iquad+1]
//        +Eprimarypt[te2p[i*4+2]*3+0]*qpt[4*iquad+2]+Eprimarypt[te2p[i*4+3]*3+0]*qpt[4*iquad+3]);
//Eprimary[3*i+1]=(Eprimarypt[te2p[i*4+0]*3+1]*qpt[4*iquad]+Eprimarypt[te2p[i*4+1]*3+1]*qpt[4*iquad+1]
//        +Eprimarypt[te2p[i*4+2]*3+1]*qpt[4*iquad+2]+Eprimarypt[te2p[i*4+3]*3+1]*qpt[4*iquad+3]);
//Eprimary[3*i+2]=(Eprimarypt[te2p[i*4+0]*3+2]*qpt[4*iquad]+Eprimarypt[te2p[i*4+1]*3+2]*qpt[4*iquad+1]
//        +Eprimarypt[te2p[i*4+2]*3+2]*qpt[4*iquad+2]+Eprimarypt[te2p[i*4+3]*3+2]*qpt[4*iquad+3]);
//}

for (i=0;i<3*nte;i++) {
Eprimary[i]=Eprimary[i]*0.1666666666666666*qwt[iquad];
}
  for (i=0;i<nte;i++) {


      Gtemp[0]=p[te2p[i*4+1]*3+0]-p[te2p[i*4]*3+0];
      Gtemp[1]=p[te2p[i*4+2]*3+0]-p[te2p[i*4]*3+0];
      Gtemp[2]=p[te2p[i*4+3]*3+0]-p[te2p[i*4]*3+0];
      Gtemp[3]=p[te2p[i*4+1]*3+1]-p[te2p[i*4]*3+1];
      Gtemp[4]=p[te2p[i*4+2]*3+1]-p[te2p[i*4]*3+1];
      Gtemp[5]=p[te2p[i*4+3]*3+1]-p[te2p[i*4]*3+1];
      Gtemp[6]=p[te2p[i*4+1]*3+2]-p[te2p[i*4]*3+2];
      Gtemp[7]=p[te2p[i*4+2]*3+2]-p[te2p[i*4]*3+2];
      Gtemp[8]=p[te2p[i*4+3]*3+2]-p[te2p[i*4]*3+2];
      Gtemp2[3 ]=Gtemp[4]*Gtemp[8]-Gtemp[7]*Gtemp[5];
      Gtemp2[4 ]=Gtemp[7]*Gtemp[2]-Gtemp[1]*Gtemp[8];
      Gtemp2[5 ]=Gtemp[1]*Gtemp[5]-Gtemp[4]*Gtemp[2];
      Gtemp2[6 ]=Gtemp[6]*Gtemp[5]-Gtemp[3]*Gtemp[8];
      Gtemp2[7 ]=Gtemp[0]*Gtemp[8]-Gtemp[6]*Gtemp[2];
      Gtemp2[8 ]=Gtemp[3]*Gtemp[2]-Gtemp[0]*Gtemp[5];
      Gtemp2[9 ]=Gtemp[3]*Gtemp[7]-Gtemp[6]*Gtemp[4];
      Gtemp2[10]=Gtemp[6]*Gtemp[1]-Gtemp[0]*Gtemp[7];
      Gtemp2[11]=Gtemp[0]*Gtemp[4]-Gtemp[3]*Gtemp[1];
     // volume=1/(Gtemp[0]*Gtemp2[3]+Gtemp[3]*Gtemp2[4]+Gtemp[6]*Gtemp2[5]);
      Gtemp2[0]=-Gtemp2[3]-Gtemp2[6]-Gtemp2[9];
      Gtemp2[1]=-Gtemp2[4]-Gtemp2[7]-Gtemp2[10];
      Gtemp2[2]=-Gtemp2[5]-Gtemp2[8]-Gtemp2[11];


jv[0]=reg[0+9*i]*Eprimary[3*i+0]+reg[3+9*i]*Eprimary[3*i+1]+reg[6+9*i]*Eprimary[3*i+2];
jv[1]=reg[1+9*i]*Eprimary[3*i+0]+reg[4+9*i]*Eprimary[3*i+1]+reg[7+9*i]*Eprimary[3*i+2];
jv[2]=reg[2+9*i]*Eprimary[3*i+0]+reg[5+9*i]*Eprimary[3*i+1]+reg[8+9*i]*Eprimary[3*i+2];



  for (k=0;k<4;k++)
    {
      rhs[te2p[k+4*i]]=rhs[te2p[k+4*i]]+(jv[0]*Gtemp2[3*k]+jv[1]*Gtemp2[1+3*k]+jv[2]*Gtemp2[2+3*k]);
    }

}
}

}


// these will be done later!

//void rhsfunc2ndks(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *rs,double* js, double* rhs,int iprec) {
//  const double mu0=1.25663706e-6;

//}





//void rhsfunc3rdks(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *rs,double* js, double* rhs,int iprec) {
//  const double mu0=1.25663706e-6;

//}
