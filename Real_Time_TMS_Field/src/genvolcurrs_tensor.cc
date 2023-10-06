#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <math.h>
#include <complex.h>
#include <algorithm>


extern "C" void conductioncurrTensor(const  int* te2p,const  int nte,const double* p
  ,const double* sigma, const double* u,  double* rv,  double* jv);
extern "C" void conductioncurr2ndTensor(const  int* te2p,const  int nte,const double* p
  ,const double* sigma, const double* u,  double* rv,  double* jv);
extern "C" void conductioncurr3rdTensor(const  int* te2p,const  int nte,const double* p
  ,const double* sigma, const double* u,  double* rv,  double* jv);



void conductioncurrTensor(const  int* te2p,const  int nte,const double* p
  ,const double* sigma, const double* u,  double* rv,  double* jv)
{

int i,j,k;
double Gtemp[12];
double Gtemp2[12];
double volume;
int te_stride=4;


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



Gtemp[0  ]=-(Gtemp2[0]*u[te2p[i*4]]
                    +Gtemp2[3]*u[te2p[i*4+1]]
                    +Gtemp2[6]*u[te2p[i*4+2]]
                    +Gtemp2[9]*u[te2p[i*4+3]])*0.16666666666666666;
Gtemp[1]=-(Gtemp2[1]*u[te2p[i*4]]
                    +Gtemp2[4]*u[te2p[i*4+1]]
                    +Gtemp2[7]*u[te2p[i*4+2]]
                    +Gtemp2[10]*u[te2p[i*4+3]])*0.16666666666666666;
Gtemp[2]=-(Gtemp2[2]*u[te2p[i*4]]
                    +Gtemp2[5]*u[te2p[i*4+1]]
                    +Gtemp2[8]*u[te2p[i*4+2]]
                    +Gtemp2[11]*u[te2p[i*4+3]])*0.16666666666666666;

jv[3*i  ]=jv[3*i  ]+sigma[0+9*i]*Gtemp[0]+sigma[3+9*i]*Gtemp[1]+sigma[6+9*i]*Gtemp[2];
jv[3*i+1]=jv[3*i+1]+sigma[1+9*i]*Gtemp[0]+sigma[4+9*i]*Gtemp[1]+sigma[7+9*i]*Gtemp[2];
jv[3*i+2]=jv[3*i+2]+sigma[2+9*i]*Gtemp[0]+sigma[5+9*i]*Gtemp[1]+sigma[8+9*i]*Gtemp[2];
rv[3*i+0]=(p[te2p[i*4+0]*3+0]+p[te2p[i*4+1]*3+0]+p[te2p[i*4+2]*3+0]+p[te2p[i*4+3]*3+0])*0.25;
rv[3*i+1]=(p[te2p[i*4+0]*3+1]+p[te2p[i*4+1]*3+1]+p[te2p[i*4+2]*3+1]+p[te2p[i*4+3]*3+1])*0.25;
rv[3*i+2]=(p[te2p[i*4+0]*3+2]+p[te2p[i*4+1]*3+2]+p[te2p[i*4+2]*3+2]+p[te2p[i*4+3]*3+2])*0.25;

}

}




void conductioncurr2ndTensor(const  int* te2p,const  int nte,const double* p
  ,const double* sigma, const double* u,  double* rv,  double* jv)
{

int i,j,k;
double Gtemp[12];
double Gtemp2[12];
double volume;

int te_stride=10;



  for (i=0;i<nte;i++) {


      Gtemp[0]=p[te2p[i*10+1]*3+0]-p[te2p[i*10]*3+0];
      Gtemp[1]=p[te2p[i*10+2]*3+0]-p[te2p[i*10]*3+0];
      Gtemp[2]=p[te2p[i*10+3]*3+0]-p[te2p[i*10]*3+0];
      Gtemp[3]=p[te2p[i*10+1]*3+1]-p[te2p[i*10]*3+1];
      Gtemp[4]=p[te2p[i*10+2]*3+1]-p[te2p[i*10]*3+1];
      Gtemp[5]=p[te2p[i*10+3]*3+1]-p[te2p[i*10]*3+1];
      Gtemp[6]=p[te2p[i*10+1]*3+2]-p[te2p[i*10]*3+2];
      Gtemp[7]=p[te2p[i*10+2]*3+2]-p[te2p[i*10]*3+2];
      Gtemp[8]=p[te2p[i*10+3]*3+2]-p[te2p[i*10]*3+2];
      Gtemp2[3 ]=Gtemp[4]*Gtemp[8]-Gtemp[7]*Gtemp[5];
      Gtemp2[4 ]=Gtemp[7]*Gtemp[2]-Gtemp[1]*Gtemp[8];
      Gtemp2[5 ]=Gtemp[1]*Gtemp[5]-Gtemp[4]*Gtemp[2];
      Gtemp2[6 ]=Gtemp[6]*Gtemp[5]-Gtemp[3]*Gtemp[8];
      Gtemp2[7 ]=Gtemp[0]*Gtemp[8]-Gtemp[6]*Gtemp[2];
      Gtemp2[8 ]=Gtemp[3]*Gtemp[2]-Gtemp[0]*Gtemp[5];
      Gtemp2[9 ]=Gtemp[3]*Gtemp[7]-Gtemp[6]*Gtemp[4];
      Gtemp2[10]=Gtemp[6]*Gtemp[1]-Gtemp[0]*Gtemp[7];
      Gtemp2[11]=Gtemp[0]*Gtemp[4]-Gtemp[3]*Gtemp[1];

      Gtemp2[0]=-Gtemp2[3]-Gtemp2[6]-Gtemp2[9];
      Gtemp2[1]=-Gtemp2[4]-Gtemp2[7]-Gtemp2[10];
      Gtemp2[2]=-Gtemp2[5]-Gtemp2[8]-Gtemp2[11];
volume=0.1666666666666666;
for (j=0;j<3;j++){
Gtemp[j]=-(
 4*(0.25*Gtemp2[j+3]+0.25*Gtemp2[j+0])*u[te2p[i*10+4]]+
 4*(0.25*Gtemp2[j+6]+0.25*Gtemp2[j+0])*u[te2p[i*10+5]]+
 4*(0.25*Gtemp2[j+9]+0.25*Gtemp2[j+0])*u[te2p[i*10+6]]+
 4*(0.25*Gtemp2[j+6]+0.25*Gtemp2[j+3])*u[te2p[i*10+7]]+
 4*(0.25*Gtemp2[j+9]+0.25*Gtemp2[j+3])*u[te2p[i*10+8]]+
 4*(0.25*Gtemp2[j+9]+0.25*Gtemp2[j+6])*u[te2p[i*10+9]]
)*volume;
}


jv[3*i  ]=jv[3*i  ]+sigma[0+9*i]*Gtemp[0]+sigma[3+9*i]*Gtemp[1]+sigma[6+9*i]*Gtemp[2];
jv[3*i+1]=jv[3*i+1]+sigma[1+9*i]*Gtemp[0]+sigma[4+9*i]*Gtemp[1]+sigma[7+9*i]*Gtemp[2];
jv[3*i+2]=jv[3*i+2]+sigma[2+9*i]*Gtemp[0]+sigma[5+9*i]*Gtemp[1]+sigma[8+9*i]*Gtemp[2];

rv[3*i+0]=(p[te2p[i*10+0]*3+0]+p[te2p[i*10+1]*3+0]+p[te2p[i*10+2]*3+0]+p[te2p[i*10+3]*3+0])*0.25;
rv[3*i+1]=(p[te2p[i*10+0]*3+1]+p[te2p[i*10+1]*3+1]+p[te2p[i*10+2]*3+1]+p[te2p[i*10+3]*3+1])*0.25;
rv[3*i+2]=(p[te2p[i*10+0]*3+2]+p[te2p[i*10+1]*3+2]+p[te2p[i*10+2]*3+2]+p[te2p[i*10+3]*3+2])*0.25;


}

}





void conductioncurr3rdTensor(const  int* te2p,const  int nte,const double* p
  ,const double* sigma, const double* u,  double* rv,  double* jv)
{
int i,j,k;
double Gtemp[12];
double Gtemp2[12];
double volume;

int te_stride=20;


  for (i=0;i<nte;i++) {

      Gtemp[0]=p[te2p[i*20+1]*3+0]-p[te2p[i*20]*3+0];
      Gtemp[1]=p[te2p[i*20+2]*3+0]-p[te2p[i*20]*3+0];
      Gtemp[2]=p[te2p[i*20+3]*3+0]-p[te2p[i*20]*3+0];
      Gtemp[3]=p[te2p[i*20+1]*3+1]-p[te2p[i*20]*3+1];
      Gtemp[4]=p[te2p[i*20+2]*3+1]-p[te2p[i*20]*3+1];
      Gtemp[5]=p[te2p[i*20+3]*3+1]-p[te2p[i*20]*3+1];
      Gtemp[6]=p[te2p[i*20+1]*3+2]-p[te2p[i*20]*3+2];
      Gtemp[7]=p[te2p[i*20+2]*3+2]-p[te2p[i*20]*3+2];
      Gtemp[8]=p[te2p[i*20+3]*3+2]-p[te2p[i*20]*3+2];
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
volume=0.1666666666666666;
for (j=0;j<3;j++){
Gtemp[j]=-(
(27*0.25*0.25-18*0.25+2)*0.5*Gtemp2[j+0]*u[te2p[i*20+0]]+
(27*0.25*0.25-18*0.25+2)*0.5*Gtemp2[j+3]*u[te2p[i*20+1]]+
(27*0.25*0.25-18*0.25+2)*0.5*Gtemp2[j+6]*u[te2p[i*20+2]]+
(27*0.25*0.25-18*0.25+2)*0.5*Gtemp2[j+9]*u[te2p[i*20+3]]+
(9*(6*0.25*0.25-0.25)*0.5*Gtemp2[j+0]+
 9*(3*0.25*0.25-0.25)*0.5*Gtemp2[j+3])*u[te2p[i*20+4]]+
(9*(6*0.25*0.25-0.25)*0.5*Gtemp2[j+3]+
 9*(3*0.25*0.25-0.25)*0.5*Gtemp2[j+0])*u[te2p[i*20+5]]+
(9*(6*0.25*0.25-0.25)*0.5*Gtemp2[j+0]+
 9*(3*0.25*0.25-0.25)*0.5*Gtemp2[j+6])*u[te2p[i*20+6]]+
(9*(6*0.25*0.25-0.25)*0.5*Gtemp2[j+6]+
 9*(3*0.25*0.25-0.25)*0.5*Gtemp2[j+0])*u[te2p[i*20+7]]+
(9*(6*0.25*0.25-0.25)*0.5*Gtemp2[j+0]+
 9*(3*0.25*0.25-0.25)*0.5*Gtemp2[j+9])*u[te2p[i*20+8]]+
(9*(6*0.25*0.25-0.25)*0.5*Gtemp2[j+9]+
 9*(3*0.25*0.25-0.25)*0.5*Gtemp2[j+0])*u[te2p[i*20+9]]+
(9*(6*0.25*0.25-0.25)*0.5*Gtemp2[j+3]+
 9*(3*0.25*0.25-0.25)*0.5*Gtemp2[j+6])*u[te2p[i*20+10]]+
(9*(6*0.25*0.25-0.25)*0.5*Gtemp2[j+6]+
 9*(3*0.25*0.25-0.25)*0.5*Gtemp2[j+3])*u[te2p[i*20+11]]+
(9*(6*0.25*0.25-0.25)*0.5*Gtemp2[j+3]+
 9*(3*0.25*0.25-0.25)*0.5*Gtemp2[j+9])*u[te2p[i*20+12]]+
(9*(6*0.25*0.25-0.25)*0.5*Gtemp2[j+9]+
 9*(3*0.25*0.25-0.25)*0.5*Gtemp2[j+3])*u[te2p[i*20+13]]+
(9*(6*0.25*0.25-0.25)*0.5*Gtemp2[j+6]+
 9*(3*0.25*0.25-0.25)*0.5*Gtemp2[j+9])*u[te2p[i*20+14]]+
(9*(6*0.25*0.25-0.25)*0.5*Gtemp2[j+9]+
 9*(3*0.25*0.25-0.25)*0.5*Gtemp2[j+6])*u[te2p[i*20+15]]+
27*(0.25*0.25*Gtemp2[j+6]+
    0.25*0.25*Gtemp2[j+0]+
    0.25*0.25*Gtemp2[j+3])*u[te2p[i*20+16]]+
27*(0.25*0.25*Gtemp2[j+9]+
    0.25*0.25*Gtemp2[j+0]+
    0.25*0.25*Gtemp2[j+3])*u[te2p[i*20+17]]+
27*(0.25*0.25*Gtemp2[j+6]+
    0.25*0.25*Gtemp2[j+0]+
    0.25*0.25*Gtemp2[j+9])*u[te2p[i*20+18]]+
27*(0.25*0.25*Gtemp2[j+6]+
    0.25*0.25*Gtemp2[j+9]+
    0.25*0.25*Gtemp2[j+3])*u[te2p[i*20+19]]
)*volume;
}



jv[3*i  ]=jv[3*i  ]+sigma[0+9*i]*Gtemp[0]+sigma[3+9*i]*Gtemp[1]+sigma[6+9*i]*Gtemp[2];
jv[3*i+1]=jv[3*i+1]+sigma[1+9*i]*Gtemp[0]+sigma[4+9*i]*Gtemp[1]+sigma[7+9*i]*Gtemp[2];
jv[3*i+2]=jv[3*i+2]+sigma[2+9*i]*Gtemp[0]+sigma[5+9*i]*Gtemp[1]+sigma[8+9*i]*Gtemp[2];

rv[3*i+0]=(p[te2p[i*20+0]*3+0]+p[te2p[i*20+1]*3+0]+p[te2p[i*20+2]*3+0]+p[te2p[i*20+3]*3+0])*0.25;
rv[3*i+1]=(p[te2p[i*20+0]*3+1]+p[te2p[i*20+1]*3+1]+p[te2p[i*20+2]*3+1]+p[te2p[i*20+3]*3+1])*0.25;
rv[3*i+2]=(p[te2p[i*20+0]*3+2]+p[te2p[i*20+1]*3+2]+p[te2p[i*20+2]*3+2]+p[te2p[i*20+3]*3+2])*0.25;
}




}
