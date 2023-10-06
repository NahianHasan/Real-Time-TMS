#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <math.h>
#include <complex.h>
#include <algorithm>

extern "C" void rhsrecip2ndpmd(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double* that,int* teid,double* rhs);
extern "C" void rhsrecip3rdpmd(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double* that,int* teid,double* rhs);

extern "C" void  rhsrecipmd(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double* that,int* teid,double* rhs);




void rhsrecipmd(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *that,int* teid,double* rhs) {

int i,j,k;
double Gtemp[12];
double Gtemp2[12];
double volume;
double totalvol=0;
totalvol=0;
volume=0;
for (i=0;i<np;i++) {
rhs[i]=0;
}
  for (i=0;i<nsource;i++) {

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
      volume=(Gtemp[0]*Gtemp2[3]+Gtemp[3]*Gtemp2[4]+Gtemp[6]*Gtemp2[5]);
      Gtemp2[0]=-Gtemp2[3]-Gtemp2[6]-Gtemp2[9];
      Gtemp2[1]=-Gtemp2[4]-Gtemp2[7]-Gtemp2[10];
      Gtemp2[2]=-Gtemp2[5]-Gtemp2[8]-Gtemp2[11];
totalvol=1/volume;
that[3*i+0]=that[3*i+0]*totalvol;
that[3*i+1]=that[3*i+1]*totalvol;
that[3*i+2]=that[3*i+2]*totalvol;
for (j=0;j<4;j++){
rhs[te2p[teid[i]*4+j]]=rhs[te2p[teid[i]*4+j]]+
                                  (Gtemp2[3*j+0]*that[3*i+0]
                                  +Gtemp2[3*j+1]*that[3*i+1]
                                  +Gtemp2[3*j+2]*that[3*i+2]);
}
that[3*i+0]=that[3*i+0]*6;
that[3*i+1]=that[3*i+1]*6;
that[3*i+2]=that[3*i+2]*6;


}




}



void rhsrecip2ndpmd(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *that,int* teid,double* rhs) {
int i,j,k;
double Gtemp[12];
double Gtemp2[12];
double volume;
int nquad=5;
double totalvol;
double qwt []={-0.8,0.45,0.45,0.45,0.45,0.45};
double bari []={
 0.2500000000000000,0.2500000000000000,0.2500000000000000,0.2500000000000000
,0.5000000000000000,0.1666666666666667,0.1666666666666667,0.1666666666666667
,0.1666666666666667,0.5000000000000000,0.1666666666666667,0.1666666666666667
,0.1666666666666667,0.1666666666666667,0.5000000000000000,0.1666666666666667
,0.1666666666666667,0.1666666666666667,0.1666666666666667,0.5000000000000000};

volume=0;
totalvol=0;
for (i=0;i<np;i++) {
rhs[i]=0;
}
  for (i=0;i<nsource;i++) {
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
      volume=(Gtemp[0]*Gtemp2[3]+Gtemp[3]*Gtemp2[4]+Gtemp[6]*Gtemp2[5]);
      Gtemp2[0]=-Gtemp2[3]-Gtemp2[6]-Gtemp2[9];
      Gtemp2[1]=-Gtemp2[4]-Gtemp2[7]-Gtemp2[10];
      Gtemp2[2]=-Gtemp2[5]-Gtemp2[8]-Gtemp2[11];

totalvol=1/volume;
that[3*i+0]=that[3*i+0]*totalvol;
that[3*i+1]=that[3*i+1]*totalvol;
that[3*i+2]=that[3*i+2]*totalvol;
for (j=0;j<3;j++){
for (k=0;k<nquad;k++){
rhs[te2p[teid[i]*10+0]]=rhs[te2p[teid[i]*10+0]]+qwt[k]*(4*bari[4*k+0]-1)*Gtemp2[j+0]*that[3*i+j];
rhs[te2p[teid[i]*10+1]]=rhs[te2p[teid[i]*10+1]]+qwt[k]*(4*bari[4*k+1]-1)*Gtemp2[j+3]*that[3*i+j];
rhs[te2p[teid[i]*10+2]]=rhs[te2p[teid[i]*10+2]]+qwt[k]*(4*bari[4*k+2]-1)*Gtemp2[j+6]*that[3*i+j];
rhs[te2p[teid[i]*10+3]]=rhs[te2p[teid[i]*10+3]]+qwt[k]*(4*bari[4*k+3]-1)*Gtemp2[j+9]*that[3*i+j];
rhs[te2p[teid[i]*10+4]]=rhs[te2p[teid[i]*10+4]]+qwt[k]*4*(bari[4*k+0]*Gtemp2[j+3]+bari[4*k+1]*Gtemp2[j+0])*that[3*i+j];
rhs[te2p[teid[i]*10+5]]=rhs[te2p[teid[i]*10+5]]+qwt[k]*4*(bari[4*k+0]*Gtemp2[j+6]+bari[4*k+2]*Gtemp2[j+0])*that[3*i+j];
rhs[te2p[teid[i]*10+6]]=rhs[te2p[teid[i]*10+6]]+qwt[k]*4*(bari[4*k+0]*Gtemp2[j+9]+bari[4*k+3]*Gtemp2[j+0])*that[3*i+j];
rhs[te2p[teid[i]*10+7]]=rhs[te2p[teid[i]*10+7]]+qwt[k]*4*(bari[4*k+1]*Gtemp2[j+6]+bari[4*k+2]*Gtemp2[j+3])*that[3*i+j];
rhs[te2p[teid[i]*10+8]]=rhs[te2p[teid[i]*10+8]]+qwt[k]*4*(bari[4*k+1]*Gtemp2[j+9]+bari[4*k+3]*Gtemp2[j+3])*that[3*i+j];
rhs[te2p[teid[i]*10+9]]=rhs[te2p[teid[i]*10+9]]+qwt[k]*4*(bari[4*k+2]*Gtemp2[j+9]+bari[4*k+3]*Gtemp2[j+6])*that[3*i+j];
}

}


that[3*i+0]=that[3*i+0]*6;
that[3*i+1]=that[3*i+1]*6;
that[3*i+2]=that[3*i+2]*6;

}

}



void rhsrecip3rdpmd(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *that,int* teid,double* rhs) {

int i,j,k;
double Gtemp[12];
double Gtemp2[12];
double volume;
int nquad=5;
double totalvol;
double qwt []={-0.8,0.45,0.45,0.45,0.45,0.45};
double bari []={
 0.2500000000000000,0.2500000000000000,0.2500000000000000,0.2500000000000000
,0.5000000000000000,0.1666666666666667,0.1666666666666667,0.1666666666666667
,0.1666666666666667,0.5000000000000000,0.1666666666666667,0.1666666666666667
,0.1666666666666667,0.1666666666666667,0.5000000000000000,0.1666666666666667
,0.1666666666666667,0.1666666666666667,0.1666666666666667,0.5000000000000000};

volume=0;
totalvol=0;
for (i=0;i<np;i++) {
rhs[i]=0;
}
  for (i=0;i<nsource;i++) {
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
      volume=(Gtemp[0]*Gtemp2[3]+Gtemp[3]*Gtemp2[4]+Gtemp[6]*Gtemp2[5]);
      Gtemp2[0]=-Gtemp2[3]-Gtemp2[6]-Gtemp2[9];
      Gtemp2[1]=-Gtemp2[4]-Gtemp2[7]-Gtemp2[10];
      Gtemp2[2]=-Gtemp2[5]-Gtemp2[8]-Gtemp2[11];
totalvol=1/volume;
that[3*i+0]=that[3*i+0]*totalvol;
that[3*i+1]=that[3*i+1]*totalvol;
that[3*i+2]=that[3*i+2]*totalvol;

for (j=0;j<3;j++){
for (k=0;k<nquad;k++){





rhs[te2p[teid[i]*20+0]]=rhs[te2p[teid[i]*20+0]]+qwt[k]*(27*bari[4*k+0]*bari[4*k+0]-18*bari[4*k+0]+2)*0.5*Gtemp2[j+0]*that[3*i+j];
rhs[te2p[teid[i]*20+1]]=rhs[te2p[teid[i]*20+1]]+qwt[k]*(27*bari[4*k+1]*bari[4*k+1]-18*bari[4*k+1]+2)*0.5*Gtemp2[j+3]*that[3*i+j];
rhs[te2p[teid[i]*20+2]]=rhs[te2p[teid[i]*20+2]]+qwt[k]*(27*bari[4*k+2]*bari[4*k+2]-18*bari[4*k+2]+2)*0.5*Gtemp2[j+6]*that[3*i+j];
rhs[te2p[teid[i]*20+3]]=rhs[te2p[teid[i]*20+3]]+qwt[k]*(27*bari[4*k+3]*bari[4*k+3]-18*bari[4*k+3]+2)*0.5*Gtemp2[j+9]*that[3*i+j];
rhs[te2p[teid[i]*20+4]]=rhs[te2p[teid[i]*20+4]]+qwt[k]*(9*(6*bari[4*k+1]*bari[4*k+0]-bari[4*k+1])*0.5*Gtemp2[j+0]+
 9*(3*bari[4*k+0]*bari[4*k+0]-bari[4*k+0])*0.5*Gtemp2[j+3])*that[3*i+j];

rhs[te2p[teid[i]*20+5]]=rhs[te2p[teid[i]*20+5]]+qwt[k]*(9*(6*bari[4*k+0]*bari[4*k+1]-bari[4*k+0])*0.5*Gtemp2[j+3]+
 9*(3*bari[4*k+1]*bari[4*k+1]-bari[4*k+1])*0.5*Gtemp2[j+0])*that[3*i+j];

rhs[te2p[teid[i]*20+6]]=rhs[te2p[teid[i]*20+6]]+qwt[k]*(9*(6*bari[4*k+2]*bari[4*k+0]-bari[4*k+2])*0.5*Gtemp2[j+0]+
 9*(3*bari[4*k+0]*bari[4*k+0]-bari[4*k+0])*0.5*Gtemp2[j+6])*that[3*i+j];

rhs[te2p[teid[i]*20+7]]=rhs[te2p[teid[i]*20+7]]+qwt[k]*(9*(6*bari[4*k+0]*bari[4*k+2]-bari[4*k+0])*0.5*Gtemp2[j+6]+
 9*(3*bari[4*k+2]*bari[4*k+2]-bari[4*k+2])*0.5*Gtemp2[j+0])*that[3*i+j];

rhs[te2p[teid[i]*20+8]]=rhs[te2p[teid[i]*20+8]]+qwt[k]*(9*(6*bari[4*k+3]*bari[4*k+0]-bari[4*k+3])*0.5*Gtemp2[j+0]+
 9*(3*bari[4*k+0]*bari[4*k+0]-bari[4*k+0])*0.5*Gtemp2[j+9])*that[3*i+j];

rhs[te2p[teid[i]*20+9]]=rhs[te2p[teid[i]*20+9]]+qwt[k]*(9*(6*bari[4*k+0]*bari[4*k+3]-bari[4*k+0])*0.5*Gtemp2[j+9]+
 9*(3*bari[4*k+3]*bari[4*k+3]-bari[4*k+3])*0.5*Gtemp2[j+0])*that[3*i+j];

rhs[te2p[teid[i]*20+10]]=rhs[te2p[teid[i]*20+10]]+qwt[k]*(9*(6*bari[4*k+2]*bari[4*k+1]-bari[4*k+2])*0.5*Gtemp2[j+3]+
 9*(3*bari[4*k+1]*bari[4*k+1]-bari[4*k+1])*0.5*Gtemp2[j+6])*that[3*i+j];

rhs[te2p[teid[i]*20+11]]=rhs[te2p[teid[i]*20+11]]+qwt[k]*(9*(6*bari[4*k+1]*bari[4*k+2]-bari[4*k+1])*0.5*Gtemp2[j+6]+
 9*(3*bari[4*k+2]*bari[4*k+2]-bari[4*k+2])*0.5*Gtemp2[j+3])*that[3*i+j];

rhs[te2p[teid[i]*20+12]]=rhs[te2p[teid[i]*20+12]]+qwt[k]*(9*(6*bari[4*k+3]*bari[4*k+1]-bari[4*k+3])*0.5*Gtemp2[j+3]+
 9*(3*bari[4*k+1]*bari[4*k+1]-bari[4*k+1])*0.5*Gtemp2[j+9])*that[3*i+j];

rhs[te2p[teid[i]*20+13]]=rhs[te2p[teid[i]*20+13]]+qwt[k]*(9*(6*bari[4*k+1]*bari[4*k+3]-bari[4*k+1])*0.5*Gtemp2[j+9]+
 9*(3*bari[4*k+3]*bari[4*k+3]-bari[4*k+3])*0.5*Gtemp2[j+3])*that[3*i+j];

rhs[te2p[teid[i]*20+14]]=rhs[te2p[teid[i]*20+14]]+qwt[k]*(9*(6*bari[4*k+3]*bari[4*k+2]-bari[4*k+3])*0.5*Gtemp2[j+6]+
 9*(3*bari[4*k+2]*bari[4*k+2]-bari[4*k+2])*0.5*Gtemp2[j+9])*that[3*i+j];

rhs[te2p[teid[i]*20+15]]=rhs[te2p[teid[i]*20+15]]+qwt[k]*(9*(6*bari[4*k+2]*bari[4*k+3]-bari[4*k+2])*0.5*Gtemp2[j+9]+
 9*(3*bari[4*k+3]*bari[4*k+3]-bari[4*k+3])*0.5*Gtemp2[j+6])*that[3*i+j];

rhs[te2p[teid[i]*20+16]]=rhs[te2p[teid[i]*20+16]]+qwt[k]*27*(bari[4*k+0]*bari[4*k+1]*Gtemp2[j+6]+
    bari[4*k+1]*bari[4*k+2]*Gtemp2[j+0]+
    bari[4*k+2]*bari[4*k+0]*Gtemp2[j+3])*that[3*i+j];

rhs[te2p[teid[i]*20+17]]=rhs[te2p[teid[i]*20+17]]+qwt[k]*27*(bari[4*k+0]*bari[4*k+1]*Gtemp2[j+9]+
    bari[4*k+1]*bari[4*k+3]*Gtemp2[j+0]+
    bari[4*k+3]*bari[4*k+0]*Gtemp2[j+3])*that[3*i+j];

rhs[te2p[teid[i]*20+18]]=rhs[te2p[teid[i]*20+18]]+qwt[k]*27*(bari[4*k+0]*bari[4*k+3]*Gtemp2[j+6]+
    bari[4*k+3]*bari[4*k+2]*Gtemp2[j+0]+
    bari[4*k+2]*bari[4*k+0]*Gtemp2[j+9])*that[3*i+j];

rhs[te2p[teid[i]*20+19]]=rhs[te2p[teid[i]*20+19]]+qwt[k]*27*(bari[4*k+3]*bari[4*k+1]*Gtemp2[j+6]+
    bari[4*k+1]*bari[4*k+2]*Gtemp2[j+9]+
    bari[4*k+2]*bari[4*k+3]*Gtemp2[j+3])*that[3*i+j];




}
}
that[3*i+0]=that[3*i+0]*6;
that[3*i+1]=that[3*i+1]*6;
that[3*i+2]=that[3*i+2]*6;
}

}
