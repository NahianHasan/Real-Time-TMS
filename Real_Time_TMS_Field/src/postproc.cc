#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <cmath>
#include <complex.h>
#include <algorithm>

extern "C" void  computeeprimary_(const double *rs,const double *js,const double* robs,double* Eprimary,const int *ntarget,const int *npoints,const double *iprec);

extern "C" void  computehprimary_(const double *rs,const double *js,const double* robs,double* Eprimary,const int *ntarget,const int *npoints,const double *iprec);

extern "C" void pointLocation_f(const  double* robs,const  int nobs,const  int* te2p,const  int nte,const double* p
  ,const int te_stride,const int* te2te,  int* teid, double* bari);
extern "C" void secondaryEfield(const  double* robs,const  int nobs,const int* te2te,const  int* te2p,const  int nte,const double* p
  ,const double* u,  double* soln,const double iprec);
extern "C" void secondaryEfield2nd(const  double* robs,const  int nobs,const int* te2te,const  int* te2p,const  int nte,const double* p
  ,const double* u,  double* soln,const double iprec);
extern "C" void secondaryEfield3rd(const  double* robs,const  int nobs,const int* te2te,const  int* te2p,const  int nte,const double* p
  ,const double* u,  double* soln,const double iprec);




extern "C" void secondaryEfieldc(const  double* robs,const  int nobs,const int* te2te,const  int* te2p,const  int nte,const double* p
  ,const double* u,  double* soln,const double iprec);
extern "C" void secondaryEfield2ndc(const  double* robs,const  int nobs,const int* te2te,const  int* te2p,const  int nte,const double* p
  ,const double* u,  double* soln,const double iprec);
extern "C" void secondaryEfield3rdc(const  double* robs,const  int nobs,const int* te2te,const  int* te2p,const  int nte,const double* p
  ,const double* u,  double* soln,const double iprec);


inline double fdet(const double x1,const double x2,const double x3,const double x4,
           const double y1,const double y2,const double y3,const double y4,
           const double z1,const double z2,const double z3,const double z4) {
return (x2*(y3*z4-y4*z3)-x3*(y2*z4-y4*z2)+x4*(y2*z3-y3*z2))
      -(x1*(y3*z4-y4*z3)-x3*(y1*z4-y4*z1)+x4*(y1*z3-y3*z1))
      +(x1*(y2*z4-y4*z2)-x2*(y1*z4-y4*z1)+x4*(y1*z2-y2*z1))
      -(x1*(y2*z3-y3*z2)-x2*(y1*z3-y3*z1)+x3*(y1*z2-y2*z1));
}


void pointLocation_f(const  double* robs,const  int nobs,const  int* te2p,const  int nte,const double* p
  ,const int te_stride,const int* te2te,  int* teid, double* bari)
{

double *tecen = (double *)malloc(sizeof(double)*3*nte);




int i;
int j;
double volume;


for (i=0;i<nte;i++){
tecen[3*i+0]=(p[te2p[i*te_stride+0]*3+0]+p[te2p[i*te_stride+1]*3+0]+p[te2p[i*te_stride+2]*3+0]+p[te2p[i*te_stride+3]*3+0])*0.25;
tecen[3*i+1]=(p[te2p[i*te_stride+0]*3+1]+p[te2p[i*te_stride+1]*3+1]+p[te2p[i*te_stride+2]*3+1]+p[te2p[i*te_stride+3]*3+1])*0.25;
tecen[3*i+2]=(p[te2p[i*te_stride+0]*3+2]+p[te2p[i*te_stride+1]*3+2]+p[te2p[i*te_stride+2]*3+2]+p[te2p[i*te_stride+3]*3+2])*0.25;

}

for (i=0;i<4*nobs;i++){
bari[i]=-1;
}

for (i=0;i<nobs;i++){
teid[i]=0;
}

  for (i=0;i<nobs;i++) {

  if ((i!=0)&&(teid[i-1]!=-1)) teid[i]=teid[i-1];
  if (nobs==nte){
teid[i]=i;
}
while (((bari[4*i]<-1e-16) || (bari[4*i+1]<-1e-16) || (bari[4*i+2]<-1e-16) || (bari[4*i+3]<-1e-16))&&(teid[i]!=-1))  {
volume=fdet(
p[te2p[teid[i]*te_stride+0]*3+0],p[te2p[teid[i]*te_stride+1]*3+0],p[te2p[teid[i]*te_stride+2]*3+0],p[te2p[teid[i]*te_stride+3]*3+0],
p[te2p[teid[i]*te_stride+0]*3+1],p[te2p[teid[i]*te_stride+1]*3+1],p[te2p[teid[i]*te_stride+2]*3+1],p[te2p[teid[i]*te_stride+3]*3+1],
p[te2p[teid[i]*te_stride+0]*3+2],p[te2p[teid[i]*te_stride+1]*3+2],p[te2p[teid[i]*te_stride+2]*3+2],p[te2p[teid[i]*te_stride+3]*3+2]);
volume=1/volume;

bari[4*i]=fdet(
robs[3*i+0],p[te2p[teid[i]*te_stride+1]*3+0],p[te2p[teid[i]*te_stride+2]*3+0],p[te2p[teid[i]*te_stride+3]*3+0],
robs[3*i+1],p[te2p[teid[i]*te_stride+1]*3+1],p[te2p[teid[i]*te_stride+2]*3+1],p[te2p[teid[i]*te_stride+3]*3+1],
robs[3*i+2],p[te2p[teid[i]*te_stride+1]*3+2],p[te2p[teid[i]*te_stride+2]*3+2],p[te2p[teid[i]*te_stride+3]*3+2]);
bari[4*i]=bari[4*i]*volume;

bari[4*i+1]=fdet(
p[te2p[teid[i]*te_stride+0]*3+0],robs[3*i+0],p[te2p[teid[i]*te_stride+2]*3+0],p[te2p[teid[i]*te_stride+3]*3+0],
p[te2p[teid[i]*te_stride+0]*3+1],robs[3*i+1],p[te2p[teid[i]*te_stride+2]*3+1],p[te2p[teid[i]*te_stride+3]*3+1],
p[te2p[teid[i]*te_stride+0]*3+2],robs[3*i+2],p[te2p[teid[i]*te_stride+2]*3+2],p[te2p[teid[i]*te_stride+3]*3+2]);
bari[4*i+1]=bari[4*i+1]*volume;

bari[4*i+2]=fdet(
p[te2p[teid[i]*te_stride+0]*3+0],p[te2p[teid[i]*te_stride+1]*3+0],robs[3*i+0],p[te2p[teid[i]*te_stride+3]*3+0],
p[te2p[teid[i]*te_stride+0]*3+1],p[te2p[teid[i]*te_stride+1]*3+1],robs[3*i+1],p[te2p[teid[i]*te_stride+3]*3+1],
p[te2p[teid[i]*te_stride+0]*3+2],p[te2p[teid[i]*te_stride+1]*3+2],robs[3*i+2],p[te2p[teid[i]*te_stride+3]*3+2]);
bari[4*i+2]=bari[4*i+2]*volume;

bari[4*i+3]=fdet(
p[te2p[teid[i]*te_stride+0]*3+0],p[te2p[teid[i]*te_stride+1]*3+0],p[te2p[teid[i]*te_stride+2]*3+0],robs[3*i+0],
p[te2p[teid[i]*te_stride+0]*3+1],p[te2p[teid[i]*te_stride+1]*3+1],p[te2p[teid[i]*te_stride+2]*3+1],robs[3*i+1],
p[te2p[teid[i]*te_stride+0]*3+2],p[te2p[teid[i]*te_stride+1]*3+2],p[te2p[teid[i]*te_stride+2]*3+2],robs[3*i+2]);
bari[4*i+3]=bari[4*i+3]*volume;

  if (bari[4*i]<-1e-12) {
teid[i]=te2te[4*teid[i]];
}
  else if (bari[4*i+1]<-1e-12) {
teid[i]=te2te[4*teid[i]+1];
}
 else if (bari[4*i+2]<-1e-12) {
teid[i]=te2te[4*teid[i]+2];
}
 else if (bari[4*i+3]<-1e-12) {
teid[i]=te2te[4*teid[i]+3];
}

}
if (teid[i]==-1){

for (j=0;j<nte;j++){
if ((fabs(robs[3*i+0]-tecen[3*j])+fabs(robs[3*i+1]-tecen[3*i+1])+fabs(robs[3*i+2]-tecen[3*j+2]))<.005) {
teid[i]=j;
volume=fdet(
p[te2p[teid[i]*te_stride+0]*3+0],p[te2p[teid[i]*te_stride+1]*3+0],p[te2p[teid[i]*te_stride+2]*3+0],p[te2p[teid[i]*te_stride+3]*3+0],
p[te2p[teid[i]*te_stride+0]*3+1],p[te2p[teid[i]*te_stride+1]*3+1],p[te2p[teid[i]*te_stride+2]*3+1],p[te2p[teid[i]*te_stride+3]*3+1],
p[te2p[teid[i]*te_stride+0]*3+2],p[te2p[teid[i]*te_stride+1]*3+2],p[te2p[teid[i]*te_stride+2]*3+2],p[te2p[teid[i]*te_stride+3]*3+2]);
volume=1/volume;

bari[4*i]=fdet(
robs[3*i+0],p[te2p[teid[i]*te_stride+1]*3+0],p[te2p[teid[i]*te_stride+2]*3+0],p[te2p[teid[i]*te_stride+3]*3+0],
robs[3*i+1],p[te2p[teid[i]*te_stride+1]*3+1],p[te2p[teid[i]*te_stride+2]*3+1],p[te2p[teid[i]*te_stride+3]*3+1],
robs[3*i+2],p[te2p[teid[i]*te_stride+1]*3+2],p[te2p[teid[i]*te_stride+2]*3+2],p[te2p[teid[i]*te_stride+3]*3+2]);
bari[4*i]=bari[4*i]*volume;

bari[4*i+1]=fdet(
p[te2p[teid[i]*te_stride+0]*3+0],robs[3*i+0],p[te2p[teid[i]*te_stride+2]*3+0],p[te2p[teid[i]*te_stride+3]*3+0],
p[te2p[teid[i]*te_stride+0]*3+1],robs[3*i+1],p[te2p[teid[i]*te_stride+2]*3+1],p[te2p[teid[i]*te_stride+3]*3+1],
p[te2p[teid[i]*te_stride+0]*3+2],robs[3*i+2],p[te2p[teid[i]*te_stride+2]*3+2],p[te2p[teid[i]*te_stride+3]*3+2]);
bari[4*i+1]=bari[4*i+1]*volume;

bari[4*i+2]=fdet(
p[te2p[teid[i]*te_stride+0]*3+0],p[te2p[teid[i]*te_stride+1]*3+0],robs[3*i+0],p[te2p[teid[i]*te_stride+3]*3+0],
p[te2p[teid[i]*te_stride+0]*3+1],p[te2p[teid[i]*te_stride+1]*3+1],robs[3*i+1],p[te2p[teid[i]*te_stride+3]*3+1],
p[te2p[teid[i]*te_stride+0]*3+2],p[te2p[teid[i]*te_stride+1]*3+2],robs[3*i+2],p[te2p[teid[i]*te_stride+3]*3+2]);
bari[4*i+2]=bari[4*i+2]*volume;

bari[4*i+3]=fdet(
p[te2p[teid[i]*te_stride+0]*3+0],p[te2p[teid[i]*te_stride+1]*3+0],p[te2p[teid[i]*te_stride+2]*3+0],robs[3*i+0],
p[te2p[teid[i]*te_stride+0]*3+1],p[te2p[teid[i]*te_stride+1]*3+1],p[te2p[teid[i]*te_stride+2]*3+1],robs[3*i+1],
p[te2p[teid[i]*te_stride+0]*3+2],p[te2p[teid[i]*te_stride+1]*3+2],p[te2p[teid[i]*te_stride+2]*3+2],robs[3*i+2]);
bari[4*i+3]=bari[4*i+3]*volume;

if ((bari[4*i]>-1e-16) && (bari[4*i+1]>-1e-16) && (bari[4*i+2]>-1e-16) && (bari[4*i+3]>-1e-16)) {

break;
}

}

}

}

}

}


void secondaryEfield(const  double* robs,const  int nobs,const int* te2te,const  int* te2p,const  int nte,const double* p
  ,const double* u, double* soln,const double iprec)
{
int i,j,k;
double Gtemp[12];
double Gtemp2[12];
double volume;
  for (i=0;i<3*nobs;i++) {
soln[i]=0;
}
int te_stride=4;
int *teid = (int *)malloc(sizeof(int)*nobs);
double *bari = (double *)malloc(sizeof(double)*4*nobs);
pointLocation_f(robs,nobs,te2p,nte,p,te_stride,te2te,teid,bari);


  for (i=0;i<nobs;i++) {

if (teid[i]!=-1) {
      Gtemp[0]=p[te2p[teid[i]*4+1]*3+0]-p[te2p[teid[i]*4]*3+0];
      Gtemp[1]=p[te2p[teid[i]*4+2]*3+0]-p[te2p[teid[i]*4]*3+0];
      Gtemp[2]=p[te2p[teid[i]*4+3]*3+0]-p[te2p[teid[i]*4]*3+0];
      Gtemp[3]=p[te2p[teid[i]*4+1]*3+1]-p[te2p[teid[i]*4]*3+1];
      Gtemp[4]=p[te2p[teid[i]*4+2]*3+1]-p[te2p[teid[i]*4]*3+1];
      Gtemp[5]=p[te2p[teid[i]*4+3]*3+1]-p[te2p[teid[i]*4]*3+1];
      Gtemp[6]=p[te2p[teid[i]*4+1]*3+2]-p[te2p[teid[i]*4]*3+2];
      Gtemp[7]=p[te2p[teid[i]*4+2]*3+2]-p[te2p[teid[i]*4]*3+2];
      Gtemp[8]=p[te2p[teid[i]*4+3]*3+2]-p[te2p[teid[i]*4]*3+2];
      Gtemp2[3 ]=Gtemp[4]*Gtemp[8]-Gtemp[7]*Gtemp[5];
      Gtemp2[4 ]=Gtemp[7]*Gtemp[2]-Gtemp[1]*Gtemp[8];
      Gtemp2[5 ]=Gtemp[1]*Gtemp[5]-Gtemp[4]*Gtemp[2];
      Gtemp2[6 ]=Gtemp[6]*Gtemp[5]-Gtemp[3]*Gtemp[8];
      Gtemp2[7 ]=Gtemp[0]*Gtemp[8]-Gtemp[6]*Gtemp[2];
      Gtemp2[8 ]=Gtemp[3]*Gtemp[2]-Gtemp[0]*Gtemp[5];
      Gtemp2[9 ]=Gtemp[3]*Gtemp[7]-Gtemp[6]*Gtemp[4];
      Gtemp2[10]=Gtemp[6]*Gtemp[1]-Gtemp[0]*Gtemp[7];
      Gtemp2[11]=Gtemp[0]*Gtemp[4]-Gtemp[3]*Gtemp[1];
      volume=1/(Gtemp[0]*Gtemp2[3]+Gtemp[3]*Gtemp2[4]+Gtemp[6]*Gtemp2[5]);
      Gtemp2[0]=-Gtemp2[3]-Gtemp2[6]-Gtemp2[9];
      Gtemp2[1]=-Gtemp2[4]-Gtemp2[7]-Gtemp2[10];
      Gtemp2[2]=-Gtemp2[5]-Gtemp2[8]-Gtemp2[11];
soln[3*i  ]=soln[3*i  ]-(Gtemp2[0]*u[te2p[teid[i]*4]]
                        +Gtemp2[3]*u[te2p[teid[i]*4+1]]
                        +Gtemp2[6]*u[te2p[teid[i]*4+2]]
                        +Gtemp2[9]*u[te2p[teid[i]*4+3]])*volume;
soln[3*i+1]=soln[3*i+1]-(Gtemp2[1]*u[te2p[teid[i]*4]]
                        +Gtemp2[4]*u[te2p[teid[i]*4+1]]
                        +Gtemp2[7]*u[te2p[teid[i]*4+2]]
                       +Gtemp2[10]*u[te2p[teid[i]*4+3]])*volume;
soln[3*i+2]=soln[3*i+2]-(Gtemp2[2]*u[te2p[teid[i]*4]]
                        +Gtemp2[5]*u[te2p[teid[i]*4+1]]
                        +Gtemp2[8]*u[te2p[teid[i]*4+2]]
                       +Gtemp2[11]*u[te2p[teid[i]*4+3]])*volume;
}
}

}




void secondaryEfield2nd(const  double* robs,const  int nobs,const int* te2te,const  int* te2p,const  int nte,const double* p
  ,const double* u, double* soln,const double iprec)
{

int i,j,k;
double Gtemp[12];
double Gtemp2[12];
double volume;
  for (i=0;i<3*nobs;i++) {
soln[i]=0;
}
int te_stride=10;
int *teid = (int *)malloc(sizeof(int)*nobs);
double *bari = (double *)malloc(sizeof(double)*4*nobs);
pointLocation_f(robs,nobs,te2p,nte,p,te_stride,te2te,teid,bari);



  for (i=0;i<nobs;i++) {

if (teid[i]!=-1) {
      Gtemp[0]=p[te2p[teid[i]*10+1]*3+0]-p[te2p[teid[i]*10]*3+0];
      Gtemp[1]=p[te2p[teid[i]*10+2]*3+0]-p[te2p[teid[i]*10]*3+0];
      Gtemp[2]=p[te2p[teid[i]*10+3]*3+0]-p[te2p[teid[i]*10]*3+0];
      Gtemp[3]=p[te2p[teid[i]*10+1]*3+1]-p[te2p[teid[i]*10]*3+1];
      Gtemp[4]=p[te2p[teid[i]*10+2]*3+1]-p[te2p[teid[i]*10]*3+1];
      Gtemp[5]=p[te2p[teid[i]*10+3]*3+1]-p[te2p[teid[i]*10]*3+1];
      Gtemp[6]=p[te2p[teid[i]*10+1]*3+2]-p[te2p[teid[i]*10]*3+2];
      Gtemp[7]=p[te2p[teid[i]*10+2]*3+2]-p[te2p[teid[i]*10]*3+2];
      Gtemp[8]=p[te2p[teid[i]*10+3]*3+2]-p[te2p[teid[i]*10]*3+2];
      Gtemp2[3 ]=Gtemp[4]*Gtemp[8]-Gtemp[7]*Gtemp[5];
      Gtemp2[4 ]=Gtemp[7]*Gtemp[2]-Gtemp[1]*Gtemp[8];
      Gtemp2[5 ]=Gtemp[1]*Gtemp[5]-Gtemp[4]*Gtemp[2];
      Gtemp2[6 ]=Gtemp[6]*Gtemp[5]-Gtemp[3]*Gtemp[8];
      Gtemp2[7 ]=Gtemp[0]*Gtemp[8]-Gtemp[6]*Gtemp[2];
      Gtemp2[8 ]=Gtemp[3]*Gtemp[2]-Gtemp[0]*Gtemp[5];
      Gtemp2[9 ]=Gtemp[3]*Gtemp[7]-Gtemp[6]*Gtemp[4];
      Gtemp2[10]=Gtemp[6]*Gtemp[1]-Gtemp[0]*Gtemp[7];
      Gtemp2[11]=Gtemp[0]*Gtemp[4]-Gtemp[3]*Gtemp[1];
      volume=1/(Gtemp[0]*Gtemp2[3]+Gtemp[3]*Gtemp2[4]+Gtemp[6]*Gtemp2[5]);
      Gtemp2[0]=-Gtemp2[3]-Gtemp2[6]-Gtemp2[9];
      Gtemp2[1]=-Gtemp2[4]-Gtemp2[7]-Gtemp2[10];
      Gtemp2[2]=-Gtemp2[5]-Gtemp2[8]-Gtemp2[11];
for (j=0;j<3;j++){
soln[3*i+j]=soln[3*i+j]-(
(4*bari[4*i+0]-1)*Gtemp2[j+0]*u[te2p[teid[i]*10+0]]+
(4*bari[4*i+1]-1)*Gtemp2[j+3]*u[te2p[teid[i]*10+1]]+
(4*bari[4*i+2]-1)*Gtemp2[j+6]*u[te2p[teid[i]*10+2]]+
(4*bari[4*i+3]-1)*Gtemp2[j+9]*u[te2p[teid[i]*10+3]]+
 4*(bari[4*i+0]*Gtemp2[j+3]+bari[4*i+1]*Gtemp2[j+0])*u[te2p[teid[i]*10+4]]+
 4*(bari[4*i+0]*Gtemp2[j+6]+bari[4*i+2]*Gtemp2[j+0])*u[te2p[teid[i]*10+5]]+
 4*(bari[4*i+0]*Gtemp2[j+9]+bari[4*i+3]*Gtemp2[j+0])*u[te2p[teid[i]*10+6]]+
 4*(bari[4*i+1]*Gtemp2[j+6]+bari[4*i+2]*Gtemp2[j+3])*u[te2p[teid[i]*10+7]]+
 4*(bari[4*i+1]*Gtemp2[j+9]+bari[4*i+3]*Gtemp2[j+3])*u[te2p[teid[i]*10+8]]+
 4*(bari[4*i+2]*Gtemp2[j+9]+bari[4*i+3]*Gtemp2[j+6])*u[te2p[teid[i]*10+9]]
)*volume;
}


}
}

}





void secondaryEfield3rd(const  double* robs,const  int nobs,const int* te2te,const  int* te2p,const  int nte,const double* p
  ,const double* u, double* soln,const double iprec)
{
int i,j,k;
double Gtemp[12];
double Gtemp2[12];
double volume;
  for (i=0;i<3*nobs;i++) {
soln[i]=0;
}
int te_stride=20;
int *teid = (int *)malloc(sizeof(int)*nobs);
double *bari = (double *)malloc(sizeof(double)*4*nobs);
pointLocation_f(robs,nobs,te2p,nte,p,te_stride,te2te,teid,bari);



  for (i=0;i<nobs;i++) {

if (teid[i]!=-1) {
      Gtemp[0]=p[te2p[teid[i]*20+1]*3+0]-p[te2p[teid[i]*20]*3+0];
      Gtemp[1]=p[te2p[teid[i]*20+2]*3+0]-p[te2p[teid[i]*20]*3+0];
      Gtemp[2]=p[te2p[teid[i]*20+3]*3+0]-p[te2p[teid[i]*20]*3+0];
      Gtemp[3]=p[te2p[teid[i]*20+1]*3+1]-p[te2p[teid[i]*20]*3+1];
      Gtemp[4]=p[te2p[teid[i]*20+2]*3+1]-p[te2p[teid[i]*20]*3+1];
      Gtemp[5]=p[te2p[teid[i]*20+3]*3+1]-p[te2p[teid[i]*20]*3+1];
      Gtemp[6]=p[te2p[teid[i]*20+1]*3+2]-p[te2p[teid[i]*20]*3+2];
      Gtemp[7]=p[te2p[teid[i]*20+2]*3+2]-p[te2p[teid[i]*20]*3+2];
      Gtemp[8]=p[te2p[teid[i]*20+3]*3+2]-p[te2p[teid[i]*20]*3+2];
      Gtemp2[3 ]=Gtemp[4]*Gtemp[8]-Gtemp[7]*Gtemp[5];
      Gtemp2[4 ]=Gtemp[7]*Gtemp[2]-Gtemp[1]*Gtemp[8];
      Gtemp2[5 ]=Gtemp[1]*Gtemp[5]-Gtemp[4]*Gtemp[2];
      Gtemp2[6 ]=Gtemp[6]*Gtemp[5]-Gtemp[3]*Gtemp[8];
      Gtemp2[7 ]=Gtemp[0]*Gtemp[8]-Gtemp[6]*Gtemp[2];
      Gtemp2[8 ]=Gtemp[3]*Gtemp[2]-Gtemp[0]*Gtemp[5];
      Gtemp2[9 ]=Gtemp[3]*Gtemp[7]-Gtemp[6]*Gtemp[4];
      Gtemp2[10]=Gtemp[6]*Gtemp[1]-Gtemp[0]*Gtemp[7];
      Gtemp2[11]=Gtemp[0]*Gtemp[4]-Gtemp[3]*Gtemp[1];
      volume=1/(Gtemp[0]*Gtemp2[3]+Gtemp[3]*Gtemp2[4]+Gtemp[6]*Gtemp2[5]);
      Gtemp2[0]=-Gtemp2[3]-Gtemp2[6]-Gtemp2[9];
      Gtemp2[1]=-Gtemp2[4]-Gtemp2[7]-Gtemp2[10];
      Gtemp2[2]=-Gtemp2[5]-Gtemp2[8]-Gtemp2[11];
for (j=0;j<3;j++){
soln[3*i+j]=soln[3*i+j]-(
(27*bari[4*i+0]*bari[4*i+0]-18*bari[4*i+0]+2)*0.5*Gtemp2[j+0]*u[te2p[teid[i]*20+0]]+
(27*bari[4*i+1]*bari[4*i+1]-18*bari[4*i+1]+2)*0.5*Gtemp2[j+3]*u[te2p[teid[i]*20+1]]+
(27*bari[4*i+2]*bari[4*i+2]-18*bari[4*i+2]+2)*0.5*Gtemp2[j+6]*u[te2p[teid[i]*20+2]]+
(27*bari[4*i+3]*bari[4*i+3]-18*bari[4*i+3]+2)*0.5*Gtemp2[j+9]*u[te2p[teid[i]*20+3]]+
(9*(6*bari[4*i+1]*bari[4*i+0]-bari[4*i+1])*0.5*Gtemp2[j+0]+
 9*(3*bari[4*i+0]*bari[4*i+0]-bari[4*i+0])*0.5*Gtemp2[j+3])*u[te2p[teid[i]*20+4]]+
(9*(6*bari[4*i+0]*bari[4*i+1]-bari[4*i+0])*0.5*Gtemp2[j+3]+
 9*(3*bari[4*i+1]*bari[4*i+1]-bari[4*i+1])*0.5*Gtemp2[j+0])*u[te2p[teid[i]*20+5]]+
(9*(6*bari[4*i+2]*bari[4*i+0]-bari[4*i+2])*0.5*Gtemp2[j+0]+
 9*(3*bari[4*i+0]*bari[4*i+0]-bari[4*i+0])*0.5*Gtemp2[j+6])*u[te2p[teid[i]*20+6]]+
(9*(6*bari[4*i+0]*bari[4*i+2]-bari[4*i+0])*0.5*Gtemp2[j+6]+
 9*(3*bari[4*i+2]*bari[4*i+2]-bari[4*i+2])*0.5*Gtemp2[j+0])*u[te2p[teid[i]*20+7]]+
(9*(6*bari[4*i+3]*bari[4*i+0]-bari[4*i+3])*0.5*Gtemp2[j+0]+
 9*(3*bari[4*i+0]*bari[4*i+0]-bari[4*i+0])*0.5*Gtemp2[j+9])*u[te2p[teid[i]*20+8]]+
(9*(6*bari[4*i+0]*bari[4*i+3]-bari[4*i+0])*0.5*Gtemp2[j+9]+
 9*(3*bari[4*i+3]*bari[4*i+3]-bari[4*i+3])*0.5*Gtemp2[j+0])*u[te2p[teid[i]*20+9]]+
(9*(6*bari[4*i+2]*bari[4*i+1]-bari[4*i+2])*0.5*Gtemp2[j+3]+
 9*(3*bari[4*i+1]*bari[4*i+1]-bari[4*i+1])*0.5*Gtemp2[j+6])*u[te2p[teid[i]*20+10]]+
(9*(6*bari[4*i+1]*bari[4*i+2]-bari[4*i+1])*0.5*Gtemp2[j+6]+
 9*(3*bari[4*i+2]*bari[4*i+2]-bari[4*i+2])*0.5*Gtemp2[j+3])*u[te2p[teid[i]*20+11]]+
(9*(6*bari[4*i+3]*bari[4*i+1]-bari[4*i+3])*0.5*Gtemp2[j+3]+
 9*(3*bari[4*i+1]*bari[4*i+1]-bari[4*i+1])*0.5*Gtemp2[j+9])*u[te2p[teid[i]*20+12]]+
(9*(6*bari[4*i+1]*bari[4*i+3]-bari[4*i+1])*0.5*Gtemp2[j+9]+
 9*(3*bari[4*i+3]*bari[4*i+3]-bari[4*i+3])*0.5*Gtemp2[j+3])*u[te2p[teid[i]*20+13]]+
(9*(6*bari[4*i+3]*bari[4*i+2]-bari[4*i+3])*0.5*Gtemp2[j+6]+
 9*(3*bari[4*i+2]*bari[4*i+2]-bari[4*i+2])*0.5*Gtemp2[j+9])*u[te2p[teid[i]*20+14]]+
(9*(6*bari[4*i+2]*bari[4*i+3]-bari[4*i+2])*0.5*Gtemp2[j+9]+
 9*(3*bari[4*i+3]*bari[4*i+3]-bari[4*i+3])*0.5*Gtemp2[j+6])*u[te2p[teid[i]*20+15]]+
27*(bari[4*i+0]*bari[4*i+1]*Gtemp2[j+6]+
    bari[4*i+1]*bari[4*i+2]*Gtemp2[j+0]+
    bari[4*i+2]*bari[4*i+0]*Gtemp2[j+3])*u[te2p[teid[i]*20+16]]+
27*(bari[4*i+0]*bari[4*i+1]*Gtemp2[j+9]+
    bari[4*i+1]*bari[4*i+3]*Gtemp2[j+0]+
    bari[4*i+3]*bari[4*i+0]*Gtemp2[j+3])*u[te2p[teid[i]*20+17]]+
27*(bari[4*i+0]*bari[4*i+3]*Gtemp2[j+6]+
    bari[4*i+3]*bari[4*i+2]*Gtemp2[j+0]+
    bari[4*i+2]*bari[4*i+0]*Gtemp2[j+9])*u[te2p[teid[i]*20+18]]+
27*(bari[4*i+3]*bari[4*i+1]*Gtemp2[j+6]+
    bari[4*i+1]*bari[4*i+2]*Gtemp2[j+9]+
    bari[4*i+2]*bari[4*i+3]*Gtemp2[j+3])*u[te2p[teid[i]*20+19]]
)*volume;
}


}
}

}










void secondaryEfieldc(const  double* robs,const  int nobs,const int* te2te,const  int* te2p,const  int nte,const double* p
  ,const double* u, double* soln,const double iprec)
{
int i,j,k;
double Gtemp[12];
double Gtemp2[12];
double volume;
  for (i=0;i<3*nobs;i++) {
soln[i]=0;
}
int te_stride=4;
int *teid = (int *)malloc(sizeof(int)*nobs);
double *bari = (double *)malloc(sizeof(double)*4*nobs);

for (i=0;i<nte;i++) {
teid[i]=i;
}

for (i=0;i<4*nte;i++) {
bari[i]=0.25;
}

  for (i=0;i<nobs;i++) {

if (teid[i]!=-1) {
      Gtemp[0]=p[te2p[teid[i]*4+1]*3+0]-p[te2p[teid[i]*4]*3+0];
      Gtemp[1]=p[te2p[teid[i]*4+2]*3+0]-p[te2p[teid[i]*4]*3+0];
      Gtemp[2]=p[te2p[teid[i]*4+3]*3+0]-p[te2p[teid[i]*4]*3+0];
      Gtemp[3]=p[te2p[teid[i]*4+1]*3+1]-p[te2p[teid[i]*4]*3+1];
      Gtemp[4]=p[te2p[teid[i]*4+2]*3+1]-p[te2p[teid[i]*4]*3+1];
      Gtemp[5]=p[te2p[teid[i]*4+3]*3+1]-p[te2p[teid[i]*4]*3+1];
      Gtemp[6]=p[te2p[teid[i]*4+1]*3+2]-p[te2p[teid[i]*4]*3+2];
      Gtemp[7]=p[te2p[teid[i]*4+2]*3+2]-p[te2p[teid[i]*4]*3+2];
      Gtemp[8]=p[te2p[teid[i]*4+3]*3+2]-p[te2p[teid[i]*4]*3+2];
      Gtemp2[3 ]=Gtemp[4]*Gtemp[8]-Gtemp[7]*Gtemp[5];
      Gtemp2[4 ]=Gtemp[7]*Gtemp[2]-Gtemp[1]*Gtemp[8];
      Gtemp2[5 ]=Gtemp[1]*Gtemp[5]-Gtemp[4]*Gtemp[2];
      Gtemp2[6 ]=Gtemp[6]*Gtemp[5]-Gtemp[3]*Gtemp[8];
      Gtemp2[7 ]=Gtemp[0]*Gtemp[8]-Gtemp[6]*Gtemp[2];
      Gtemp2[8 ]=Gtemp[3]*Gtemp[2]-Gtemp[0]*Gtemp[5];
      Gtemp2[9 ]=Gtemp[3]*Gtemp[7]-Gtemp[6]*Gtemp[4];
      Gtemp2[10]=Gtemp[6]*Gtemp[1]-Gtemp[0]*Gtemp[7];
      Gtemp2[11]=Gtemp[0]*Gtemp[4]-Gtemp[3]*Gtemp[1];
      volume=1/(Gtemp[0]*Gtemp2[3]+Gtemp[3]*Gtemp2[4]+Gtemp[6]*Gtemp2[5]);
      Gtemp2[0]=-Gtemp2[3]-Gtemp2[6]-Gtemp2[9];
      Gtemp2[1]=-Gtemp2[4]-Gtemp2[7]-Gtemp2[10];
      Gtemp2[2]=-Gtemp2[5]-Gtemp2[8]-Gtemp2[11];
soln[3*i  ]=soln[3*i  ]-(Gtemp2[0]*u[te2p[teid[i]*4]]
                        +Gtemp2[3]*u[te2p[teid[i]*4+1]]
                        +Gtemp2[6]*u[te2p[teid[i]*4+2]]
                        +Gtemp2[9]*u[te2p[teid[i]*4+3]])*volume;
soln[3*i+1]=soln[3*i+1]-(Gtemp2[1]*u[te2p[teid[i]*4]]
                        +Gtemp2[4]*u[te2p[teid[i]*4+1]]
                        +Gtemp2[7]*u[te2p[teid[i]*4+2]]
                       +Gtemp2[10]*u[te2p[teid[i]*4+3]])*volume;
soln[3*i+2]=soln[3*i+2]-(Gtemp2[2]*u[te2p[teid[i]*4]]
                        +Gtemp2[5]*u[te2p[teid[i]*4+1]]
                        +Gtemp2[8]*u[te2p[teid[i]*4+2]]
                       +Gtemp2[11]*u[te2p[teid[i]*4+3]])*volume;
}
}

}




void secondaryEfield2ndc(const  double* robs,const  int nobs,const int* te2te,const  int* te2p,const  int nte,const double* p
  ,const double* u, double* soln,const double iprec)
{

int i,j,k;
double Gtemp[12];
double Gtemp2[12];
double volume;
  for (i=0;i<3*nobs;i++) {
soln[i]=0;
}
int te_stride=10;
int *teid = (int *)malloc(sizeof(int)*nobs);
double *bari = (double *)malloc(sizeof(double)*4*nobs);

for (i=0;i<nte;i++) {
teid[i]=i;
}

for (i=0;i<4*nte;i++) {
bari[i]=0.25;
}


  for (i=0;i<nobs;i++) {

if (teid[i]!=-1) {
      Gtemp[0]=p[te2p[teid[i]*10+1]*3+0]-p[te2p[teid[i]*10]*3+0];
      Gtemp[1]=p[te2p[teid[i]*10+2]*3+0]-p[te2p[teid[i]*10]*3+0];
      Gtemp[2]=p[te2p[teid[i]*10+3]*3+0]-p[te2p[teid[i]*10]*3+0];
      Gtemp[3]=p[te2p[teid[i]*10+1]*3+1]-p[te2p[teid[i]*10]*3+1];
      Gtemp[4]=p[te2p[teid[i]*10+2]*3+1]-p[te2p[teid[i]*10]*3+1];
      Gtemp[5]=p[te2p[teid[i]*10+3]*3+1]-p[te2p[teid[i]*10]*3+1];
      Gtemp[6]=p[te2p[teid[i]*10+1]*3+2]-p[te2p[teid[i]*10]*3+2];
      Gtemp[7]=p[te2p[teid[i]*10+2]*3+2]-p[te2p[teid[i]*10]*3+2];
      Gtemp[8]=p[te2p[teid[i]*10+3]*3+2]-p[te2p[teid[i]*10]*3+2];
      Gtemp2[3 ]=Gtemp[4]*Gtemp[8]-Gtemp[7]*Gtemp[5];
      Gtemp2[4 ]=Gtemp[7]*Gtemp[2]-Gtemp[1]*Gtemp[8];
      Gtemp2[5 ]=Gtemp[1]*Gtemp[5]-Gtemp[4]*Gtemp[2];
      Gtemp2[6 ]=Gtemp[6]*Gtemp[5]-Gtemp[3]*Gtemp[8];
      Gtemp2[7 ]=Gtemp[0]*Gtemp[8]-Gtemp[6]*Gtemp[2];
      Gtemp2[8 ]=Gtemp[3]*Gtemp[2]-Gtemp[0]*Gtemp[5];
      Gtemp2[9 ]=Gtemp[3]*Gtemp[7]-Gtemp[6]*Gtemp[4];
      Gtemp2[10]=Gtemp[6]*Gtemp[1]-Gtemp[0]*Gtemp[7];
      Gtemp2[11]=Gtemp[0]*Gtemp[4]-Gtemp[3]*Gtemp[1];
      volume=1/(Gtemp[0]*Gtemp2[3]+Gtemp[3]*Gtemp2[4]+Gtemp[6]*Gtemp2[5]);
      Gtemp2[0]=-Gtemp2[3]-Gtemp2[6]-Gtemp2[9];
      Gtemp2[1]=-Gtemp2[4]-Gtemp2[7]-Gtemp2[10];
      Gtemp2[2]=-Gtemp2[5]-Gtemp2[8]-Gtemp2[11];
for (j=0;j<3;j++){
soln[3*i+j]=soln[3*i+j]-(
(4*bari[4*i+0]-1)*Gtemp2[j+0]*u[te2p[teid[i]*10+0]]+
(4*bari[4*i+1]-1)*Gtemp2[j+3]*u[te2p[teid[i]*10+1]]+
(4*bari[4*i+2]-1)*Gtemp2[j+6]*u[te2p[teid[i]*10+2]]+
(4*bari[4*i+3]-1)*Gtemp2[j+9]*u[te2p[teid[i]*10+3]]+
 4*(bari[4*i+0]*Gtemp2[j+3]+bari[4*i+1]*Gtemp2[j+0])*u[te2p[teid[i]*10+4]]+
 4*(bari[4*i+0]*Gtemp2[j+6]+bari[4*i+2]*Gtemp2[j+0])*u[te2p[teid[i]*10+5]]+
 4*(bari[4*i+0]*Gtemp2[j+9]+bari[4*i+3]*Gtemp2[j+0])*u[te2p[teid[i]*10+6]]+
 4*(bari[4*i+1]*Gtemp2[j+6]+bari[4*i+2]*Gtemp2[j+3])*u[te2p[teid[i]*10+7]]+
 4*(bari[4*i+1]*Gtemp2[j+9]+bari[4*i+3]*Gtemp2[j+3])*u[te2p[teid[i]*10+8]]+
 4*(bari[4*i+2]*Gtemp2[j+9]+bari[4*i+3]*Gtemp2[j+6])*u[te2p[teid[i]*10+9]]
)*volume;
}


}
}

}





void secondaryEfield3rdc(const  double* robs,const  int nobs,const int* te2te,const  int* te2p,const  int nte,const double* p
  ,const double* u, double* soln,const double iprec)
{
int i,j,k;
double Gtemp[12];
double Gtemp2[12];
double volume;
  for (i=0;i<3*nobs;i++) {
soln[i]=0;
}
int te_stride=20;
int *teid = (int *)malloc(sizeof(int)*nobs);
double *bari = (double *)malloc(sizeof(double)*4*nobs);

for (i=0;i<nte;i++) {
teid[i]=i;
}

for (i=0;i<4*nte;i++) {
bari[i]=0.25;
}

  for (i=0;i<nobs;i++) {

if (teid[i]!=-1) {
      Gtemp[0]=p[te2p[teid[i]*20+1]*3+0]-p[te2p[teid[i]*20]*3+0];
      Gtemp[1]=p[te2p[teid[i]*20+2]*3+0]-p[te2p[teid[i]*20]*3+0];
      Gtemp[2]=p[te2p[teid[i]*20+3]*3+0]-p[te2p[teid[i]*20]*3+0];
      Gtemp[3]=p[te2p[teid[i]*20+1]*3+1]-p[te2p[teid[i]*20]*3+1];
      Gtemp[4]=p[te2p[teid[i]*20+2]*3+1]-p[te2p[teid[i]*20]*3+1];
      Gtemp[5]=p[te2p[teid[i]*20+3]*3+1]-p[te2p[teid[i]*20]*3+1];
      Gtemp[6]=p[te2p[teid[i]*20+1]*3+2]-p[te2p[teid[i]*20]*3+2];
      Gtemp[7]=p[te2p[teid[i]*20+2]*3+2]-p[te2p[teid[i]*20]*3+2];
      Gtemp[8]=p[te2p[teid[i]*20+3]*3+2]-p[te2p[teid[i]*20]*3+2];
      Gtemp2[3 ]=Gtemp[4]*Gtemp[8]-Gtemp[7]*Gtemp[5];
      Gtemp2[4 ]=Gtemp[7]*Gtemp[2]-Gtemp[1]*Gtemp[8];
      Gtemp2[5 ]=Gtemp[1]*Gtemp[5]-Gtemp[4]*Gtemp[2];
      Gtemp2[6 ]=Gtemp[6]*Gtemp[5]-Gtemp[3]*Gtemp[8];
      Gtemp2[7 ]=Gtemp[0]*Gtemp[8]-Gtemp[6]*Gtemp[2];
      Gtemp2[8 ]=Gtemp[3]*Gtemp[2]-Gtemp[0]*Gtemp[5];
      Gtemp2[9 ]=Gtemp[3]*Gtemp[7]-Gtemp[6]*Gtemp[4];
      Gtemp2[10]=Gtemp[6]*Gtemp[1]-Gtemp[0]*Gtemp[7];
      Gtemp2[11]=Gtemp[0]*Gtemp[4]-Gtemp[3]*Gtemp[1];
      volume=1/(Gtemp[0]*Gtemp2[3]+Gtemp[3]*Gtemp2[4]+Gtemp[6]*Gtemp2[5]);
      Gtemp2[0]=-Gtemp2[3]-Gtemp2[6]-Gtemp2[9];
      Gtemp2[1]=-Gtemp2[4]-Gtemp2[7]-Gtemp2[10];
      Gtemp2[2]=-Gtemp2[5]-Gtemp2[8]-Gtemp2[11];
for (j=0;j<3;j++){
soln[3*i+j]=soln[3*i+j]-(
(27*bari[4*i+0]*bari[4*i+0]-18*bari[4*i+0]+2)*0.5*Gtemp2[j+0]*u[te2p[teid[i]*20+0]]+
(27*bari[4*i+1]*bari[4*i+1]-18*bari[4*i+1]+2)*0.5*Gtemp2[j+3]*u[te2p[teid[i]*20+1]]+
(27*bari[4*i+2]*bari[4*i+2]-18*bari[4*i+2]+2)*0.5*Gtemp2[j+6]*u[te2p[teid[i]*20+2]]+
(27*bari[4*i+3]*bari[4*i+3]-18*bari[4*i+3]+2)*0.5*Gtemp2[j+9]*u[te2p[teid[i]*20+3]]+
(9*(6*bari[4*i+1]*bari[4*i+0]-bari[4*i+1])*0.5*Gtemp2[j+0]+
 9*(3*bari[4*i+0]*bari[4*i+0]-bari[4*i+0])*0.5*Gtemp2[j+3])*u[te2p[teid[i]*20+4]]+
(9*(6*bari[4*i+0]*bari[4*i+1]-bari[4*i+0])*0.5*Gtemp2[j+3]+
 9*(3*bari[4*i+1]*bari[4*i+1]-bari[4*i+1])*0.5*Gtemp2[j+0])*u[te2p[teid[i]*20+5]]+
(9*(6*bari[4*i+2]*bari[4*i+0]-bari[4*i+2])*0.5*Gtemp2[j+0]+
 9*(3*bari[4*i+0]*bari[4*i+0]-bari[4*i+0])*0.5*Gtemp2[j+6])*u[te2p[teid[i]*20+6]]+
(9*(6*bari[4*i+0]*bari[4*i+2]-bari[4*i+0])*0.5*Gtemp2[j+6]+
 9*(3*bari[4*i+2]*bari[4*i+2]-bari[4*i+2])*0.5*Gtemp2[j+0])*u[te2p[teid[i]*20+7]]+
(9*(6*bari[4*i+3]*bari[4*i+0]-bari[4*i+3])*0.5*Gtemp2[j+0]+
 9*(3*bari[4*i+0]*bari[4*i+0]-bari[4*i+0])*0.5*Gtemp2[j+9])*u[te2p[teid[i]*20+8]]+
(9*(6*bari[4*i+0]*bari[4*i+3]-bari[4*i+0])*0.5*Gtemp2[j+9]+
 9*(3*bari[4*i+3]*bari[4*i+3]-bari[4*i+3])*0.5*Gtemp2[j+0])*u[te2p[teid[i]*20+9]]+
(9*(6*bari[4*i+2]*bari[4*i+1]-bari[4*i+2])*0.5*Gtemp2[j+3]+
 9*(3*bari[4*i+1]*bari[4*i+1]-bari[4*i+1])*0.5*Gtemp2[j+6])*u[te2p[teid[i]*20+10]]+
(9*(6*bari[4*i+1]*bari[4*i+2]-bari[4*i+1])*0.5*Gtemp2[j+6]+
 9*(3*bari[4*i+2]*bari[4*i+2]-bari[4*i+2])*0.5*Gtemp2[j+3])*u[te2p[teid[i]*20+11]]+
(9*(6*bari[4*i+3]*bari[4*i+1]-bari[4*i+3])*0.5*Gtemp2[j+3]+
 9*(3*bari[4*i+1]*bari[4*i+1]-bari[4*i+1])*0.5*Gtemp2[j+9])*u[te2p[teid[i]*20+12]]+
(9*(6*bari[4*i+1]*bari[4*i+3]-bari[4*i+1])*0.5*Gtemp2[j+9]+
 9*(3*bari[4*i+3]*bari[4*i+3]-bari[4*i+3])*0.5*Gtemp2[j+3])*u[te2p[teid[i]*20+13]]+
(9*(6*bari[4*i+3]*bari[4*i+2]-bari[4*i+3])*0.5*Gtemp2[j+6]+
 9*(3*bari[4*i+2]*bari[4*i+2]-bari[4*i+2])*0.5*Gtemp2[j+9])*u[te2p[teid[i]*20+14]]+
(9*(6*bari[4*i+2]*bari[4*i+3]-bari[4*i+2])*0.5*Gtemp2[j+9]+
 9*(3*bari[4*i+3]*bari[4*i+3]-bari[4*i+3])*0.5*Gtemp2[j+6])*u[te2p[teid[i]*20+15]]+
27*(bari[4*i+0]*bari[4*i+1]*Gtemp2[j+6]+
    bari[4*i+1]*bari[4*i+2]*Gtemp2[j+0]+
    bari[4*i+2]*bari[4*i+0]*Gtemp2[j+3])*u[te2p[teid[i]*20+16]]+
27*(bari[4*i+0]*bari[4*i+1]*Gtemp2[j+9]+
    bari[4*i+1]*bari[4*i+3]*Gtemp2[j+0]+
    bari[4*i+3]*bari[4*i+0]*Gtemp2[j+3])*u[te2p[teid[i]*20+17]]+
27*(bari[4*i+0]*bari[4*i+3]*Gtemp2[j+6]+
    bari[4*i+3]*bari[4*i+2]*Gtemp2[j+0]+
    bari[4*i+2]*bari[4*i+0]*Gtemp2[j+9])*u[te2p[teid[i]*20+18]]+
27*(bari[4*i+3]*bari[4*i+1]*Gtemp2[j+6]+
    bari[4*i+1]*bari[4*i+2]*Gtemp2[j+9]+
    bari[4*i+2]*bari[4*i+3]*Gtemp2[j+3])*u[te2p[teid[i]*20+19]]
)*volume;
}


}
}

}
