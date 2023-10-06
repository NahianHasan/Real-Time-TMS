#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <math.h>
#include <complex.h>
#include <algorithm>
extern "C" void  computehprimary_(double *rs,double *js,double* robs,double* Eprimary,int *ntarget,int *npoints,double *iprec);
extern "C" void rhsfuncks(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *rs,double* js,double* rhs,double iprec);
extern "C" void rhsfunc2ndks(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *rs,double* js,double* rhs,double iprec);
extern "C" void rhsfunc3rdks(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *rs,double* js,double* rhs,double iprec);
extern "C" int compare_faces(const void * a, const void * b);


template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void rhsfuncks(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *rs,double* js, double* rhs,double iprec) {
  const double mu0=1.25663706e-6;

  //all this just makes sorted triangles with lists of tetra pointers
  int faces[3];
  int *faces2 = (int *)malloc(sizeof(int)*16*nte);
  const int sh[12]={1,2,3,2,3,0,3,0,1,0,1,2};
  int i,j,k,ii;
  for (i=0;i<np;i++) {
rhs[i]=0;
}
  //sort matrix entries
  for (k=0;k<nte;k++)  {
    for (j=0;j<4;j++)  {
      i=4*k+j;
      faces[0]=te2p[sh[3*j  ]+k*4];
      faces[1]=te2p[sh[3*j+1]+k*4];
      faces[2]=te2p[sh[3*j+2]+k*4];

      if (faces[0]>faces[1]){
	if (faces[1]>faces[2]){
	  faces2[4*i  ]=faces[2];  //third is smallest
	  faces2[4*i+1]=faces[1];
	  faces2[4*i+2]=faces[0];
	}
	else {
	    faces2[4*i  ]=faces[1];//second is smallest
	  if (faces[0]>faces[2]){
	    faces2[4*i+1]=faces[2];
	    faces2[4*i+2]=faces[0];
	  }
	  else{
	    faces2[4*i+1]=faces[0];
	    faces2[4*i+2]=faces[2];
	  }
	}
      }
      else {
//face 1 is bigger
	if (faces[0]>faces[2]){
	  faces2[4*i  ]=faces[2];  //third is smallest
	  faces2[4*i+1]=faces[0];
	  faces2[4*i+2]=faces[1];
	}
	else {
	    faces2[4*i  ]=faces[0];//first is smallest
	  if (faces[1]>faces[2]){
	    faces2[4*i+1]=faces[2];
	    faces2[4*i+2]=faces[1];
	  }
	  else{
	    faces2[4*i+1]=faces[1];
	    faces2[4*i+2]=faces[2];
	  }
	}
      }
      faces2[4*i+3]=i; //global id of triangle
    }
  }
  qsort((void*)faces2, 4*nte, 4*sizeof(int), compare_faces); //sorting triangles in ascending order

   int *ff = (int *)malloc(sizeof(int)*4*nte);
   int *ff2 = (int *)malloc(sizeof(int)*4*nte);
  int nf,ct2;

  nf=0;
  ff[faces2[3]]=nf;//ordered by 4*tetra pointing to face id
  ff2[nf]=0;//face to face2points
  for (i=1; i<4*nte; i++) {

    if (faces2[4*i  ]==faces2[4*(i-1)  ] &&
        faces2[4*i+1]==faces2[4*(i-1)+1] &&
        faces2[4*i+2]==faces2[4*(i-1)+2]){
    ff[faces2[4*i+3]]=nf;
    }
    else {//write a unique face
      nf++;
      ff[faces2[4*i+3]]=nf;
      ff2[nf]=i;//points to first occurrence of a face
    }
  }
  nf++;
  int *f2te = (int *)malloc(sizeof(int)*2*nf);//unique face 2 tetra array
  int *freenode = (int *)malloc(sizeof(int)*nf);//unique face 2 tetra array
  for (i=0; i<2*nf; i++) f2te[i]=-1;

  for (i=0; i<nte; i++) {
    for (j=0; j<4; j++) {
      if (f2te[2*ff[j+4*i]]==-1){
	f2te[2*ff[j+4*i]]=i;
	freenode[ff[j+4*i]]=te2p[4*i+j];
      }
      else {
	f2te[2*ff[j+4*i]+1]=i;
      }
    }
  }

 //get ready to actually compute E-primary

 int ninterface;
ct2=0;
 for (i=0; i<nf; i++) {
if (f2te[2*i+1]!=-1){
if(reg[f2te[2*i+1]]!=reg[f2te[2*i]]) ct2++;
}
else{
ct2++;}
}

 ninterface=ct2;
//choose favorite quadrature
int nquad=3;
double qwt []={0.3333333333333333,0.3333333333333333,0.3333333333333333};
double qpt []={0.1666666666666666,0.1666666666666666,0.6666666666666666,
               0.1666666666666666,0.6666666666666666,0.1666666666666666,
               0.6666666666666666,0.1666666666666666,0.1666666666666666};

                     int ntarget=nquad*ninterface;
                     double *robs = (double *)malloc(sizeof(double)*ntarget*3);//unique face 2 tetra array
                     int *faceindex = (int *)malloc(sizeof(int)*ninterface);//unique face 2 tetra array

                     ct2=0;
                     for (i=0;i<nf;i++) {
                     if (f2te[2*i+1]==-1){
                      faceindex[ct2]=i;
                       for (j=0;j<nquad;j++) {
                           robs[  3*(j+nquad*ct2)]=qpt[3*j  ]*p[3*faces2[4*ff2[i]  ]  ]+
                                                   qpt[3*j+1]*p[3*faces2[4*ff2[i]+1]  ]+
                                                   qpt[3*j+2]*p[3*faces2[4*ff2[i]+2]  ];
                           robs[1+3*(j+nquad*ct2)]=qpt[3*j  ]*p[3*faces2[4*ff2[i]  ]+1]+
                                                   qpt[3*j+1]*p[3*faces2[4*ff2[i]+1]+1]+
                                                   qpt[3*j+2]*p[3*faces2[4*ff2[i]+2]+1];
                           robs[2+3*(j+nquad*ct2)]=qpt[3*j  ]*p[3*faces2[4*ff2[i]  ]+2]+
                                                   qpt[3*j+1]*p[3*faces2[4*ff2[i]+1]+2]+
                                                   qpt[3*j+2]*p[3*faces2[4*ff2[i]+2]+2];
                       }

                       ct2++;
                     }
                     else{
                     if(reg[f2te[2*i+1]]!=reg[f2te[2*i]]){
                      faceindex[ct2]=i;
                       for (j=0;j<nquad;j++) {
                           robs[  3*(j+nquad*ct2)]=qpt[3*j  ]*p[3*faces2[4*ff2[i]  ]  ]+
                                                   qpt[3*j+1]*p[3*faces2[4*ff2[i]+1]  ]+
                                                   qpt[3*j+2]*p[3*faces2[4*ff2[i]+2]  ];
                           robs[1+3*(j+nquad*ct2)]=qpt[3*j  ]*p[3*faces2[4*ff2[i]  ]+1]+
                                                   qpt[3*j+1]*p[3*faces2[4*ff2[i]+1]+1]+
                                                   qpt[3*j+2]*p[3*faces2[4*ff2[i]+2]+1];
                           robs[2+3*(j+nquad*ct2)]=qpt[3*j  ]*p[3*faces2[4*ff2[i]  ]+2]+
                                                   qpt[3*j+1]*p[3*faces2[4*ff2[i]+1]+2]+
                                                   qpt[3*j+2]*p[3*faces2[4*ff2[i]+2]+2];
                       }

                       ct2++;
                     }
                     }
                     }



                     double *Eprimary = (double *)malloc(sizeof(double)*ntarget*3);//unique face 2 tetra array
                     computehprimary_(rs,js,robs,Eprimary,&ntarget,&nsource,&iprec);

                     double v1[3],v2[3],outn[3],nhat[3];
                     double dadtnhat,del;


for (j=0; j<ninterface; j++) {
    i=faceindex[j];//i now points to correct face
    v1[0]=p[3*faces2[4*ff2[i]  ]  ]-p[3*faces2[4*ff2[i]+2]  ];
    v1[1]=p[3*faces2[4*ff2[i]  ]+1]-p[3*faces2[4*ff2[i]+2]+1];
    v1[2]=p[3*faces2[4*ff2[i]  ]+2]-p[3*faces2[4*ff2[i]+2]+2];

    v2[0]=p[3*faces2[4*ff2[i]+1]  ]-p[3*faces2[4*ff2[i]+2]  ];
    v2[1]=p[3*faces2[4*ff2[i]+1]+1]-p[3*faces2[4*ff2[i]+2]+1];
    v2[2]=p[3*faces2[4*ff2[i]+1]+2]-p[3*faces2[4*ff2[i]+2]+2];

    outn[0]=p[3*faces2[4*ff2[i]+1]  ]-p[3*freenode[i]  ];
    outn[1]=p[3*faces2[4*ff2[i]+1]+1]-p[3*freenode[i]+1];
    outn[2]=p[3*faces2[4*ff2[i]+1]+2]-p[3*freenode[i]+2];
    nhat[0]=v1[1]*v2[2]-v1[2]*v2[1];
    nhat[1]=v1[2]*v2[0]-v1[0]*v2[2];
    nhat[2]=v1[0]*v2[1]-v1[1]*v2[0];
if (f2te[2*i+1]==-1){
del=reg[f2te[2*i]];
}
else {
del=(-reg[f2te[2*i+1]]+reg[f2te[2*i]]);
}

del=del*sgn(nhat[0]*outn[0]+nhat[1]*outn[1]+nhat[2]*outn[2]);
nhat[0]=nhat[0]*del;
nhat[1]=nhat[1]*del;
nhat[2]=nhat[2]*del;

for (k=0; k<nquad; k++) {
dadtnhat=(Eprimary[  3*(k+nquad*j)]*nhat[0]+
          Eprimary[1+3*(k+nquad*j)]*nhat[1]+
          Eprimary[2+3*(k+nquad*j)]*nhat[2]);

//rhs[3*j]=Eprimary[3*(k+nquad*j)];
//rhs[3*j+1]=Eprimary[3*(k+nquad*j)+1];
//rhs[3*j+2]=Eprimary[3*(k+nquad*j)+2];
    for (ii=0; ii<3; ii++) rhs[faces2[4*ff2[i]+ii]]=rhs[faces2[4*ff2[i]+ii]]
        + qpt[3*k+ii]*dadtnhat*qwt[k]*0.5;

  }


  }

}



void rhsfunc2ndks(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *rs,double* js, double* rhs,double iprec) {
  const double mu0=1.25663706e-6;

  //all this just makes sorted triangles with lists of tetra pointers
  int faces[6];
  int *faces2 = (int *)malloc(sizeof(int)*28*nte);
  const int sh[24]={1,2,3,7,8,9
                   ,2,3,0,9,5,6
                   ,3,0,1,6,8,4
                   ,0,1,2,4,5,7};//get triangles in correct order
  int i,j,k,ii;
  for (i=0;i<np;i++) {
rhs[i]=0;
}
  //sort matrix entries
  for (k=0;k<nte;k++)  {
    for (j=0;j<4;j++)  {
      i=4*k+j;
      faces[0]=te2p[sh[6*j  ]+k*10];
      faces[1]=te2p[sh[6*j+1]+k*10];
      faces[2]=te2p[sh[6*j+2]+k*10];
      faces[3]=te2p[sh[6*j+3]+k*10];
      faces[4]=te2p[sh[6*j+4]+k*10];
      faces[5]=te2p[sh[6*j+5]+k*10];

      if (faces[0]>faces[1]){
	if (faces[1]>faces[2]){
	  faces2[7*i  ]=faces[2];  //third is smallest
	  faces2[7*i+1]=faces[1];
	  faces2[7*i+2]=faces[0];
    faces2[7*i+3]=faces[5];
	  faces2[7*i+4]=faces[4];
    faces2[7*i+5]=faces[3];
	}
	else {
	  if (faces[0]>faces[2]){
	  faces2[7*i  ]=faces[1];//second is smallest
	  faces2[7*i+1]=faces[2];
	  faces2[7*i+2]=faces[0];
          faces2[7*i+3]=faces[5];
	  faces2[7*i+4]=faces[3];
          faces2[7*i+5]=faces[4];
	  }
	  else{
	  faces2[7*i  ]=faces[1];//second is smallest
	  faces2[7*i+1]=faces[0];
	  faces2[7*i+2]=faces[2];
          faces2[7*i+3]=faces[3];
	  faces2[7*i+4]=faces[5];
          faces2[7*i+5]=faces[4];
	  }
	}
      }
      else {
//face 1 is bigger
	if (faces[0]>faces[2]){
	  faces2[7*i  ]=faces[2];  //third is smallest
	  faces2[7*i+1]=faces[0];
	  faces2[7*i+2]=faces[1];
          faces2[7*i+3]=faces[4];
	  faces2[7*i+4]=faces[5];
          faces2[7*i+5]=faces[3];
	}
	else {
	  if (faces[1]>faces[2]){
            faces2[7*i  ]=faces[0];//first is smallest
	    faces2[7*i+1]=faces[2];
	    faces2[7*i+2]=faces[1];
            faces2[7*i+3]=faces[4];
	    faces2[7*i+4]=faces[3];
            faces2[7*i+5]=faces[5];
	  }
	  else{
            faces2[7*i  ]=faces[0];//first is smallest
	    faces2[7*i+1]=faces[1];
	    faces2[7*i+2]=faces[2];
            faces2[7*i+3]=faces[3];
	    faces2[7*i+4]=faces[4];
            faces2[7*i+5]=faces[5];
	  }
	}
      }
      faces2[7*i+6]=i; //global id of triangle
    }
  }
  qsort((void*)faces2, 4*nte, 7*sizeof(int), compare_faces); //sorting triangles in ascending order


   int *ff = (int *)malloc(sizeof(int)*4*nte);
   int *ff2 = (int *)malloc(sizeof(int)*4*nte);
  int nf,ct2;

  nf=0;
  ff[faces2[6]]=nf;//ordered by 4*tetra pointing to face id
  ff2[nf]=0;//face to face2points
  for (i=1; i<4*nte; i++) {

    if (faces2[7*i  ]==faces2[7*(i-1)  ] &&
        faces2[7*i+1]==faces2[7*(i-1)+1] &&
        faces2[7*i+2]==faces2[7*(i-1)+2]){
    ff[faces2[7*i+6]]=nf;
    }
    else {//write a unique face
      nf++;
      ff[faces2[7*i+6]]=nf;
      ff2[nf]=i;//points to first occurrence of a face
    }
  }
  nf++;
  int *f2te = (int *)malloc(sizeof(int)*2*nf);//unique face 2 tetra array
  int *freenode = (int *)malloc(sizeof(int)*nf);//unique face 2 tetra array
  for (i=0; i<2*nf; i++) f2te[i]=-1;

  for (i=0; i<nte; i++) {
    for (j=0; j<4; j++) {
      if (f2te[2*ff[j+4*i]]==-1){
	f2te[2*ff[j+4*i]]=i;
	freenode[ff[j+4*i]]=te2p[10*i+j];
      }
      else {
	f2te[2*ff[j+4*i]+1]=i;
      }
    }
  }

 //get ready to actually compute E-primary

 int ninterface;
ct2=0;
 for (i=0; i<nf; i++) {
if (f2te[2*i+1]!=-1){
if(reg[f2te[2*i+1]]!=reg[f2te[2*i]]) ct2++;
}
else{
ct2++;}
}

 ninterface=ct2;
//choose favorite quadrature

 int nquad=4;
 double qwt []={-0.5625000000000000,0.5208333333333300,0.5208333333333300,0.5208333333333300};
 double qpt []={0.3333333333333333,0.3333333333333333,0.3333333333333333,
                      0.2,0.2,0.6,
                      0.2,0.6,0.2,
                      0.6,0.2,0.2};

                      //

 int ntarget=nquad*ninterface;
 double *robs = (double *)malloc(sizeof(double)*ntarget*3);//unique face 2 tetra array
 int *faceindex = (int *)malloc(sizeof(int)*ninterface);//unique face 2 tetra array

ct2=0;
 for (i=0;i<nf;i++) {
if (f2te[2*i+1]==-1){
	faceindex[ct2]=i;
   for (j=0;j<nquad;j++) {
       robs[  3*(j+nquad*ct2)]=qpt[3*j  ]*p[3*faces2[7*ff2[i]  ]  ]+
                               qpt[3*j+1]*p[3*faces2[7*ff2[i]+1]  ]+
			                         qpt[3*j+2]*p[3*faces2[7*ff2[i]+2]  ];
       robs[1+3*(j+nquad*ct2)]=qpt[3*j  ]*p[3*faces2[7*ff2[i]  ]+1]+
			                         qpt[3*j+1]*p[3*faces2[7*ff2[i]+1]+1]+
			                         qpt[3*j+2]*p[3*faces2[7*ff2[i]+2]+1];
       robs[2+3*(j+nquad*ct2)]=qpt[3*j  ]*p[3*faces2[7*ff2[i]  ]+2]+
			                         qpt[3*j+1]*p[3*faces2[7*ff2[i]+1]+2]+
			                         qpt[3*j+2]*p[3*faces2[7*ff2[i]+2]+2];
	 }
   ct2++;
}
else{
if(reg[f2te[2*i+1]]!=reg[f2te[2*i]]){
	faceindex[ct2]=i;
   for (j=0;j<nquad;j++) {
       robs[  3*(j+nquad*ct2)]=qpt[3*j  ]*p[3*faces2[7*ff2[i]  ]  ]+
                               qpt[3*j+1]*p[3*faces2[7*ff2[i]+1]  ]+
			                         qpt[3*j+2]*p[3*faces2[7*ff2[i]+2]  ];
       robs[1+3*(j+nquad*ct2)]=qpt[3*j  ]*p[3*faces2[7*ff2[i]  ]+1]+
			                         qpt[3*j+1]*p[3*faces2[7*ff2[i]+1]+1]+
			                         qpt[3*j+2]*p[3*faces2[7*ff2[i]+2]+1];
       robs[2+3*(j+nquad*ct2)]=qpt[3*j  ]*p[3*faces2[7*ff2[i]  ]+2]+
			                         qpt[3*j+1]*p[3*faces2[7*ff2[i]+1]+2]+
			                         qpt[3*j+2]*p[3*faces2[7*ff2[i]+2]+2];
	 }
   ct2++;
}
}
}


 double *Eprimary = (double *)malloc(sizeof(double)*ntarget*3);//unique face 2 tetra array
computehprimary_(rs,js,robs,Eprimary,&ntarget,&nsource,&iprec);

double v1[3],v2[3],outn[3],nhat[3];
double dadtnhat,del;
double bak[6*nquad];

for (j=0; j<nquad; j++) {
bak[6*j+0]=qpt[3*j  ]*(2*qpt[3*j  ]-1);
bak[6*j+1]=qpt[3*j+1]*(2*qpt[3*j+1]-1);
bak[6*j+2]=qpt[3*j+2]*(2*qpt[3*j+2]-1);
bak[6*j+3]=4*qpt[3*j  ]*qpt[3*j+1];
bak[6*j+4]=4*qpt[3*j  ]*qpt[3*j+2];
bak[6*j+5]=4*qpt[3*j+1]*qpt[3*j+2];
}

for (j=0; j<ninterface; j++) {
    i=faceindex[j];//i now points to correct face
    v1[0]=p[3*faces2[7*ff2[i]  ]  ]-p[3*faces2[7*ff2[i]+2]  ];
    v1[1]=p[3*faces2[7*ff2[i]  ]+1]-p[3*faces2[7*ff2[i]+2]+1];
    v1[2]=p[3*faces2[7*ff2[i]  ]+2]-p[3*faces2[7*ff2[i]+2]+2];

    v2[0]=p[3*faces2[7*ff2[i]+1]  ]-p[3*faces2[7*ff2[i]+2]  ];
    v2[1]=p[3*faces2[7*ff2[i]+1]+1]-p[3*faces2[7*ff2[i]+2]+1];
    v2[2]=p[3*faces2[7*ff2[i]+1]+2]-p[3*faces2[7*ff2[i]+2]+2];

    outn[0]=p[3*faces2[7*ff2[i]+1]  ]-p[3*freenode[i]  ];
    outn[1]=p[3*faces2[7*ff2[i]+1]+1]-p[3*freenode[i]+1];
    outn[2]=p[3*faces2[7*ff2[i]+1]+2]-p[3*freenode[i]+2];
    nhat[0]=v1[1]*v2[2]-v1[2]*v2[1];
    nhat[1]=v1[2]*v2[0]-v1[0]*v2[2];
    nhat[2]=v1[0]*v2[1]-v1[1]*v2[0];
if (f2te[2*i+1]==-1){
del=reg[f2te[2*i]];
}
else {
del=(-reg[f2te[2*i+1]]+reg[f2te[2*i]]);
}

del=del*sgn(nhat[0]*outn[0]+nhat[1]*outn[1]+nhat[2]*outn[2]);
nhat[0]=nhat[0]*del;
nhat[1]=nhat[1]*del;
nhat[2]=nhat[2]*del;

for (k=0; k<nquad; k++) {
dadtnhat=(Eprimary[  3*(k+nquad*j)]*nhat[0]+
          Eprimary[1+3*(k+nquad*j)]*nhat[1]+
          Eprimary[2+3*(k+nquad*j)]*nhat[2]);

    for (ii=0; ii<6; ii++) rhs[faces2[7*ff2[i]+ii]]=rhs[faces2[7*ff2[i]+ii]]
        + bak[6*k+ii]*dadtnhat*qwt[k]*0.5;

  }


  }

}





void rhsfunc3rdks(const int nte,const int np,const  int* te2p,const  double* p,const  double* reg,int nsource, double *rs,double* js, double* rhs,double iprec) {
  const double mu0=1.25663706e-6;
  //all this just makes sorted triangles with lists of tetra pointers

  const int nelem=20;
  const int nfelem=10;
  const int nfelemp=11;
  int faces[nfelem];
  int *faces2 = (int *)malloc(sizeof(int)*4*nfelemp*nte);

  const int sh[4*nfelem]={1,2,3,10,11,12,13,14,15,19
                         ,2,3,0,14,15, 7, 6, 9, 8,18
                         ,3,0,1, 9, 8,13,12, 4, 5,17
                         ,0,1,2, 4, 5, 6, 7,10,11,16};//get triangles in correct order
  int i,j,k,ii,jj;
  for (i=0;i<np;i++) {
rhs[i]=0;
}
  //sort matrix entries
  for (k=0;k<nte;k++)  {
    for (j=0;j<4;j++)  {
      i=4*k+j;
  faces[0]=te2p[sh[nfelem*j+0]+k*nelem];
  faces[1]=te2p[sh[nfelem*j+1]+k*nelem];
  faces[2]=te2p[sh[nfelem*j+2]+k*nelem];
  faces[3]=te2p[sh[nfelem*j+3]+k*nelem];
  faces[4]=te2p[sh[nfelem*j+4]+k*nelem];
  faces[5]=te2p[sh[nfelem*j+5]+k*nelem];
  faces[6]=te2p[sh[nfelem*j+6]+k*nelem];
  faces[7]=te2p[sh[nfelem*j+7]+k*nelem];
  faces[8]=te2p[sh[nfelem*j+8]+k*nelem];
  faces[9]=te2p[sh[nfelem*j+9]+k*nelem];

      if (faces[0]>faces[1]){
	if (faces[1]>faces[2]){
	  faces2[nfelemp*i  ]=faces[2];  //third is smallest
	  faces2[nfelemp*i+1]=faces[1];
	  faces2[nfelemp*i+2]=faces[0];
    faces2[nfelemp*i+3]=faces[8];
    faces2[nfelemp*i+4]=faces[7];
	  faces2[nfelemp*i+5]=faces[6];
    faces2[nfelemp*i+6]=faces[5];
    faces2[nfelemp*i+7]=faces[4];
    faces2[nfelemp*i+8]=faces[3];
	}
	else {
	  if (faces[0]>faces[2]){
	  faces2[nfelemp*i  ]=faces[1];//second is smallest
	  faces2[nfelemp*i+1]=faces[2];
	  faces2[nfelemp*i+2]=faces[0];
    faces2[nfelemp*i+3]=faces[7];
    faces2[nfelemp*i+4]=faces[8];
	  faces2[nfelemp*i+5]=faces[4];
	  faces2[nfelemp*i+6]=faces[3];
    faces2[nfelemp*i+7]=faces[6];
    faces2[nfelemp*i+8]=faces[5];
	  }
	  else{
	  faces2[nfelemp*i  ]=faces[1];//second is smallest
	  faces2[nfelemp*i+1]=faces[0];
	  faces2[nfelemp*i+2]=faces[2];
    faces2[nfelemp*i+3]=faces[4];
    faces2[nfelemp*i+4]=faces[3];
	  faces2[nfelemp*i+5]=faces[7];
	  faces2[nfelemp*i+6]=faces[8];
    faces2[nfelemp*i+7]=faces[5];
    faces2[nfelemp*i+8]=faces[6];
	  }
	}
      }
      else {
//face 1 is bigger
	if (faces[0]>faces[2]){
	  faces2[nfelemp*i  ]=faces[2];  //third is smallest
	  faces2[nfelemp*i+1]=faces[0];
	  faces2[nfelemp*i+2]=faces[1];
    faces2[nfelemp*i+3]=faces[6];
    faces2[nfelemp*i+4]=faces[5];
	  faces2[nfelemp*i+5]=faces[8];
	  faces2[nfelemp*i+6]=faces[7];
    faces2[nfelemp*i+7]=faces[3];
    faces2[nfelemp*i+8]=faces[4];
	}
	else {
	  if (faces[1]>faces[2]){
      faces2[nfelemp*i  ]=faces[0];//first is smallest
	    faces2[nfelemp*i+1]=faces[2];
	    faces2[nfelemp*i+2]=faces[1];
      faces2[nfelemp*i+3]=faces[5];
      faces2[nfelemp*i+4]=faces[6];
	    faces2[nfelemp*i+5]=faces[3];
	    faces2[nfelemp*i+6]=faces[4];
      faces2[nfelemp*i+7]=faces[8];
      faces2[nfelemp*i+8]=faces[7];
	  }
	  else{
      faces2[nfelemp*i  ]=faces[0];//first is smallest
	    faces2[nfelemp*i+1]=faces[1];
	    faces2[nfelemp*i+2]=faces[2];
      faces2[nfelemp*i+3]=faces[3];
      faces2[nfelemp*i+4]=faces[4];
	    faces2[nfelemp*i+5]=faces[5];
	    faces2[nfelemp*i+6]=faces[6];
      faces2[nfelemp*i+7]=faces[7];
      faces2[nfelemp*i+8]=faces[8];
	  }
	}
}
      faces2[nfelemp*i+9]=faces[9]; //global id of triangle
      faces2[nfelemp*i+10]=i; //global id of triangle
    }
  }
  qsort((void*)faces2, 4*nte, nfelemp*sizeof(int), compare_faces); //sorting triangles in ascending order


   int *ff = (int *)malloc(sizeof(int)*4*nte);
   int *ff2 = (int *)malloc(sizeof(int)*4*nte);
  int nf,ct2;

  nf=0;
  ff[faces2[nfelem]]=nf;//ordered by 4*tetra pointing to face id
  ff2[nf]=0;//face to face2points
  for (i=1; i<4*nte; i++) {

    if (faces2[nfelemp*i  ]==faces2[nfelemp*(i-1)  ] &&
        faces2[nfelemp*i+1]==faces2[nfelemp*(i-1)+1] &&
        faces2[nfelemp*i+2]==faces2[nfelemp*(i-1)+2]){
    ff[faces2[nfelemp*i+nfelem]]=nf;
    }
    else {//write a unique face
      nf++;
      ff[faces2[nfelemp*i+nfelem]]=nf;
      ff2[nf]=i;//points to first occurrence of a face
    }
  }
  nf++;
  int *f2te = (int *)malloc(sizeof(int)*2*nf);//unique face 2 tetra array
  int *freenode = (int *)malloc(sizeof(int)*nf);//unique face 2 tetra array
  for (i=0; i<2*nf; i++) f2te[i]=-1;

  for (i=0; i<nte; i++) {
    for (j=0; j<4; j++) {
      if (f2te[2*ff[j+4*i]]==-1){
	f2te[2*ff[j+4*i]]=i;
	freenode[ff[j+4*i]]=te2p[nelem*i+j];
      }
      else {
	f2te[2*ff[j+4*i]+1]=i;
      }
    }
  }
 //get ready to actually compute E-primary

 int ninterface;
ct2=0;
 for (i=0; i<nf; i++) {
if (f2te[2*i+1]!=-1){
if(reg[f2te[2*i+1]]!=reg[f2te[2*i]]) ct2++;
}
else{
ct2++;}
}

 ninterface=ct2;
//choose favorite quadrature

 const int nquad=6;
 const double qwt []={0.2233815896780100,0.2233815896780100,0.2233815896780100,
                          0.1099517436553200,0.1099517436553200,0.1099517436553200};
 const double qpt []={0.1081030181680600,0.4459484909159700,0.4459484909159700,
                      0.4459484909159700,0.1081030181680600,0.4459484909159700,
                      0.4459484909159700,0.4459484909159700,0.1081030181680600,
                      0.8168475729804601,0.09157621350976999,0.09157621350976999,
                      0.09157621350976999,0.8168475729804601,0.09157621350976999,
                      0.09157621350976999,0.09157621350976999,0.8168475729804601};

//

 int ntarget=nquad*ninterface;
 double *robs = (double *)malloc(sizeof(double)*ntarget*3);//unique face 2 tetra array
 int *faceindex = (int *)malloc(sizeof(int)*ninterface);//unique face 2 tetra array

ct2=0;
 for (i=0;i<nf;i++) {
if (f2te[2*i+1]==-1){
	faceindex[ct2]=i;
   for (j=0;j<nquad;j++) {
       robs[  3*(j+nquad*ct2)]=qpt[3*j  ]*p[3*faces2[nfelemp*ff2[i]  ]  ]+
                               qpt[3*j+1]*p[3*faces2[nfelemp*ff2[i]+1]  ]+
			                         qpt[3*j+2]*p[3*faces2[nfelemp*ff2[i]+2]  ];
       robs[1+3*(j+nquad*ct2)]=qpt[3*j  ]*p[3*faces2[nfelemp*ff2[i]  ]+1]+
			                         qpt[3*j+1]*p[3*faces2[nfelemp*ff2[i]+1]+1]+
			                         qpt[3*j+2]*p[3*faces2[nfelemp*ff2[i]+2]+1];
       robs[2+3*(j+nquad*ct2)]=qpt[3*j  ]*p[3*faces2[nfelemp*ff2[i]  ]+2]+
			                         qpt[3*j+1]*p[3*faces2[nfelemp*ff2[i]+1]+2]+
			                         qpt[3*j+2]*p[3*faces2[nfelemp*ff2[i]+2]+2];
	 }
   ct2++;
}
else{
if(reg[f2te[2*i+1]]!=reg[f2te[2*i]]){
	faceindex[ct2]=i;
   for (j=0;j<nquad;j++) {
       robs[  3*(j+nquad*ct2)]=qpt[3*j  ]*p[3*faces2[nfelemp*ff2[i]  ]  ]+
                               qpt[3*j+1]*p[3*faces2[nfelemp*ff2[i]+1]  ]+
			                         qpt[3*j+2]*p[3*faces2[nfelemp*ff2[i]+2]  ];
       robs[1+3*(j+nquad*ct2)]=qpt[3*j  ]*p[3*faces2[nfelemp*ff2[i]  ]+1]+
			                         qpt[3*j+1]*p[3*faces2[nfelemp*ff2[i]+1]+1]+
			                         qpt[3*j+2]*p[3*faces2[nfelemp*ff2[i]+2]+1];
       robs[2+3*(j+nquad*ct2)]=qpt[3*j  ]*p[3*faces2[nfelemp*ff2[i]  ]+2]+
			                         qpt[3*j+1]*p[3*faces2[nfelemp*ff2[i]+1]+2]+
			                         qpt[3*j+2]*p[3*faces2[nfelemp*ff2[i]+2]+2];
	 }
   ct2++;
 }
}
}



 double *Eprimary = (double *)malloc(sizeof(double)*ntarget*3);//unique face 2 tetra array
computehprimary_(rs,js,robs,Eprimary,&ntarget,&nsource,&iprec);

double v1[3],v2[3],outn[3],nhat[3];
double dadtnhat,del;
double bak[nfelem*nquad];

for (j=0; j<nquad; j++) {

bak[nfelem*j+0]=qpt[3*j  ]*(3*qpt[3*j  ]-1)*(3*qpt[3*j  ]-2)*0.5;
bak[nfelem*j+1]=qpt[3*j+1]*(3*qpt[3*j+1]-1)*(3*qpt[3*j+1]-2)*0.5;
bak[nfelem*j+2]=qpt[3*j+2]*(3*qpt[3*j+2]-1)*(3*qpt[3*j+2]-2)*0.5;
bak[nfelem*j+3]=9*qpt[3*j  ]*qpt[3*j+1]*(3*qpt[3*j  ]-1)*0.5;
bak[nfelem*j+4]=9*qpt[3*j  ]*qpt[3*j+1]*(3*qpt[3*j+1]-1)*0.5;
bak[nfelem*j+5]=9*qpt[3*j  ]*qpt[3*j+2]*(3*qpt[3*j  ]-1)*0.5;
bak[nfelem*j+6]=9*qpt[3*j  ]*qpt[3*j+2]*(3*qpt[3*j+2]-1)*0.5;
bak[nfelem*j+7]=9*qpt[3*j+1]*qpt[3*j+2]*(3*qpt[3*j+1]-1)*0.5;
bak[nfelem*j+8]=9*qpt[3*j+1]*qpt[3*j+2]*(3*qpt[3*j+2]-1)*0.5;
bak[nfelem*j+9]=27*qpt[3*j  ]*qpt[3*j+1]*qpt[3*j+2];


}

for (j=0; j<ninterface; j++) {


    i=faceindex[j];//i now points to correct face
    v1[0]=p[3*faces2[nfelemp*ff2[i]  ]  ]-p[3*faces2[nfelemp*ff2[i]+2]  ];
    v1[1]=p[3*faces2[nfelemp*ff2[i]  ]+1]-p[3*faces2[nfelemp*ff2[i]+2]+1];
    v1[2]=p[3*faces2[nfelemp*ff2[i]  ]+2]-p[3*faces2[nfelemp*ff2[i]+2]+2];

    v2[0]=p[3*faces2[nfelemp*ff2[i]+1]  ]-p[3*faces2[nfelemp*ff2[i]+2]  ];
    v2[1]=p[3*faces2[nfelemp*ff2[i]+1]+1]-p[3*faces2[nfelemp*ff2[i]+2]+1];
    v2[2]=p[3*faces2[nfelemp*ff2[i]+1]+2]-p[3*faces2[nfelemp*ff2[i]+2]+2];

    outn[0]=p[3*faces2[nfelemp*ff2[i]+1]  ]-p[3*freenode[i]  ];
    outn[1]=p[3*faces2[nfelemp*ff2[i]+1]+1]-p[3*freenode[i]+1];
    outn[2]=p[3*faces2[nfelemp*ff2[i]+1]+2]-p[3*freenode[i]+2];
    nhat[0]=v1[1]*v2[2]-v1[2]*v2[1];
    nhat[1]=v1[2]*v2[0]-v1[0]*v2[2];
    nhat[2]=v1[0]*v2[1]-v1[1]*v2[0];
if (f2te[2*i+1]==-1){
del=reg[f2te[2*i]];
}
else {
del=(-reg[f2te[2*i+1]]+reg[f2te[2*i]]);
}

del=del*sgn(nhat[0]*outn[0]+nhat[1]*outn[1]+nhat[2]*outn[2]);
nhat[0]=nhat[0]*del;
nhat[1]=nhat[1]*del;
nhat[2]=nhat[2]*del;

for (k=0; k<nquad; k++) {
dadtnhat=(Eprimary[  3*(k+nquad*j)]*nhat[0]+
          Eprimary[1+3*(k+nquad*j)]*nhat[1]+
          Eprimary[2+3*(k+nquad*j)]*nhat[2]);

//rhs[3*j]=Eprimary[3*(k+nquad*j)];
//rhs[3*j+1]=Eprimary[3*(k+nquad*j)+1];
//rhs[3*j+2]=Eprimary[3*(k+nquad*j)+2];
    for (ii=0; ii<nfelem; ii++) rhs[faces2[nfelemp*ff2[i]+ii]]=rhs[faces2[nfelemp*ff2[i]+ii]]
        + bak[nfelem*k+ii]*dadtnhat*qwt[k]*0.5;

  }


  }

}