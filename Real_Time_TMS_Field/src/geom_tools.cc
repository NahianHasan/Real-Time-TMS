
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <math.h>
#include <complex.h>
#include <algorithm>
extern "C" void generatetetrapaths(const int nte,const  int* te2p,int* te2te);

extern "C" void generatemesh3rd(const int nte,const  int* te2p,const double* p,const int np,const  int ne,int* te2p2,double* p2);
extern "C" void generatemesh2nd(const int nte,const  int* te2p,const double* p,const int np,int* te2p2,double* p2);
extern "C" int generateedge_ct(const int nte,const  int* te2p);
extern "C" int generateface_ct(const int nte,const  int* te2p);


int compare_edges_geomtools(const void * a, const void * b)
{
  int i= *(int*)a - *(int*)b;
  if (i==0) i= *((int*)a+1) - *((int*)b+1);
  return (i);
}

int compare_faces_geomtools(const void * a, const void * b)
{
  int i= *(int*)a - *(int*)b;
  if (i==0) i= *((int*)a+1) - *((int*)b+1);
  if (i==0) i= *((int*)a+2) - *((int*)b+2);
  return (i);
}


void generatetetrapaths(const int nte,const  int* te2p,int* te2te)
{

  //all this just makes sorted triangles with lists of tetra pointers
  int faces[3];
int *faces2 = (int *)malloc(sizeof(int)*16*nte);
  const int sh[12]={1,2,3,2,3,0,3,0,1,0,1,2};
  int i,j,k,ii;

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
  qsort((void*)faces2, 4*nte, 4*sizeof(int), compare_faces_geomtools); //sorting triangles in ascending order

  for (i=0;i<4*nte-1;i++)  {
if ((faces2[4*i+0]==faces2[4*(i+1)+0]) && (faces2[4*i+1]==faces2[4*(i+1)+1]) && (faces2[4*i+2]==faces2[4*(i+1)+2])) {

te2te[faces2[4*i+3]]=floor(faces2[4*(i+1)+3]*0.25);
te2te[faces2[4*(i+1)+3]]=floor(faces2[4*i+3]*0.25);
}

}

  }


int generateface_ct(const int nte,const  int* te2p)
{
int nf;
  //all this just makes sorted triangles with lists of tetra pointers
  int faces[3];
int *faces2 = (int *)malloc(sizeof(int)*16*nte);
  const int sh[12]={1,2,3,2,3,0,3,0,1,0,1,2};
  int i,j,k,ii;

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
  qsort((void*)faces2, 4*nte, 4*sizeof(int), compare_faces_geomtools); //sorting triangles in ascending order
nf=1;
  for (i=1;i<4*nte;i++)  {
if ((faces2[4*i+0]==faces2[4*(i-1)+0]) && (faces2[4*i+1]==faces2[4*(i-1)+1]) && (faces2[4*i+2]==faces2[4*(i-1)+2])) {
}
else{
nf=nf+1;
}

}
return (nf);
}


int generateedge_ct(const int nte,const  int* te2p)
{
int ne;
  //all this just makes sorted triangles with lists of tetra pointers
  int faces[2];
int *faces2 = (int *)malloc(sizeof(int)*18*nte);

  const int sh[12]={0,1,0,2,0,3,1,2,1,3,2,3};
  int i,j,k,ii;

  //sort matrix entries
  for (k=0;k<nte;k++)  {
    for (j=0;j<6;j++)  {
      i=6*k+j;
      faces[0]=te2p[sh[2*j  ]+k*4];
      faces[1]=te2p[sh[2*j+1]+k*4];

      if (faces[0]>faces[1]){
      faces2[3*i+0]=faces[1]; //global id of triangle
      faces2[3*i+1]=faces[0]; //global id of triangle
	}
else{
      faces2[3*i+0]=faces[0]; //global id of triangle
      faces2[3*i+1]=faces[1]; //global id of triangle
}
      faces2[3*i+2]=i; //global id of triangle
    }
  }

qsort((void*)faces2, 6*nte, 3*sizeof(int), compare_edges_geomtools); //sorting triangles in ascending order
ne=1;
  for (i=1;i<6*nte;i++)  {
if ((faces2[3*i+0]==faces2[3*(i-1)+0]) && (faces2[3*i+1]==faces2[3*(i-1)+1]))  {
}
else {
ne=ne+1;
}

  }
return (ne);
}


void generatemesh2nd(const int nte,const  int* te2p,const double* p,const int np, int* te2p2,double* p2)
{

  //all this just makes sorted triangles with lists of tetra pointers
  int faces[2];
int *faces2 = (int *)malloc(sizeof(int)*18*nte);

  const int sh[12]={0,1,0,2,0,3,1,2,1,3,2,3};
  int i,j,k,ii;

  //sort matrix entries
  for (k=0;k<nte;k++)  {
    for (j=0;j<6;j++)  {
      i=6*k+j;
      faces[0]=te2p[sh[2*j  ]+k*4];
      faces[1]=te2p[sh[2*j+1]+k*4];

      if (faces[0]>faces[1]){
      faces2[3*i+0]=faces[1]; //global id of triangle
      faces2[3*i+1]=faces[0]; //global id of triangle
	}
else{
      faces2[3*i+0]=faces[0]; //global id of triangle
      faces2[3*i+1]=faces[1]; //global id of triangle
}
      faces2[3*i+2]=10*k+4+j; //global id of triangle
    }
  }

qsort((void*)faces2, 6*nte, 3*sizeof(int), compare_edges_geomtools); //sorting triangles in ascending order
//write points
  for (i=0;i<3*np;i++)  {
p2[i]=p[i];
}

int ct=np;
//write midpoints
i=0;
p2[3*ct+0]=(p[3*faces2[3*i+0]+0]+p[3*faces2[3*i+1]+0])*0.5;
p2[3*ct+1]=(p[3*faces2[3*i+0]+1]+p[3*faces2[3*i+1]+1])*0.5;
p2[3*ct+2]=(p[3*faces2[3*i+0]+2]+p[3*faces2[3*i+1]+2])*0.5;
ct=ct+1;

for (i=1;i<6*nte;i++)  {
if ((faces2[3*i+0]==faces2[3*(i-1)+0]) && (faces2[3*i+1]==faces2[3*(i-1)+1]))  {
}
else {
p2[3*ct+0]=(p[3*faces2[3*i+0]+0]+p[3*faces2[3*i+1]+0])*0.5;
p2[3*ct+1]=(p[3*faces2[3*i+0]+1]+p[3*faces2[3*i+1]+1])*0.5;
p2[3*ct+2]=(p[3*faces2[3*i+0]+2]+p[3*faces2[3*i+1]+2])*0.5;
ct=ct+1;
}
}

//write tetrahedrons
  for (i=0;i<nte;i++)  {
te2p2[10*i+0]=te2p[4*i+0];
te2p2[10*i+1]=te2p[4*i+1];
te2p2[10*i+2]=te2p[4*i+2];
te2p2[10*i+3]=te2p[4*i+3];
}


//write midtetrahedrons
ct=np;
i=0;
te2p2[faces2[3*i+2]]=ct;

for (i=1;i<6*nte;i++)  {
if ((faces2[3*i+0]==faces2[3*(i-1)+0]) && (faces2[3*i+1]==faces2[3*(i-1)+1]))  {
}
else {
ct=ct+1;
}

te2p2[faces2[3*i+2]]=ct;
}

}



void generatemesh3rd(const int nte,const  int* te2p,const double* p,const int np,const int ne,int* te2p2,double* p2)
{

  //all this just makes sorted triangles with lists of tetra pointers
  int faces[3];
int *faces2 = (int *)malloc(sizeof(int)*6*4*nte);

  const int sh[12]={0,1,0,2,0,3,1,2,1,3,2,3};
  int i,j,k,ii;

//write points
  for (i=0;i<3*np;i++)  {
p2[i]=p[i];
}

//write tetrahedrons
  for (i=0;i<nte;i++)  {
te2p2[20*i+0]=te2p[4*i+0];
te2p2[20*i+1]=te2p[4*i+1];
te2p2[20*i+2]=te2p[4*i+2];
te2p2[20*i+3]=te2p[4*i+3];
}


  //sort matrix entries
  for (k=0;k<nte;k++)  {
    for (j=0;j<6;j++)  {
      i=6*k+j;
      faces[0]=te2p[sh[2*j  ]+k*4];
      faces[1]=te2p[sh[2*j+1]+k*4];

      if (faces[0]>faces[1]){
      faces2[4*i+0]=faces[1]; //global id of triangle
      faces2[4*i+1]=faces[0]; //global id of triangle
      faces2[4*i+2]=20*k+4+2*j+1; //global id of triangle
      faces2[4*i+3]=20*k+4+2*j; //global id of triangle
	}
else{
      faces2[4*i+0]=faces[0]; //global id of triangle
      faces2[4*i+1]=faces[1]; //global id of triangle
      faces2[4*i+2]=20*k+4+2*j; //global id of triangle
      faces2[4*i+3]=20*k+4+2*j+1; //global id of triangle
}

    }
  }

qsort((void*)faces2, 6*nte, 4*sizeof(int), compare_edges_geomtools); //sorting triangles in ascending order

int ct=np;
//write midpoints
i=0;
p2[3*ct+0]=(2*p[3*faces2[4*i+0]+0]+p[3*faces2[4*i+1]+0])*0.33333333333333333;
p2[3*ct+1]=(2*p[3*faces2[4*i+0]+1]+p[3*faces2[4*i+1]+1])*0.33333333333333333;
p2[3*ct+2]=(2*p[3*faces2[4*i+0]+2]+p[3*faces2[4*i+1]+2])*0.33333333333333333;

p2[3*(ct+1)+0]=(p[3*faces2[4*i+0]+0]+2*p[3*faces2[4*i+1]+0])*0.33333333333333333;
p2[3*(ct+1)+1]=(p[3*faces2[4*i+0]+1]+2*p[3*faces2[4*i+1]+1])*0.33333333333333333;
p2[3*(ct+1)+2]=(p[3*faces2[4*i+0]+2]+2*p[3*faces2[4*i+1]+2])*0.33333333333333333;
ct=ct+2;

for (i=1;i<6*nte;i++)  {
if ((faces2[4*i+0]==faces2[4*(i-1)+0]) && (faces2[4*i+1]==faces2[4*(i-1)+1]))  {
}
else {
p2[3*ct+0]=(2*p[3*faces2[4*i+0]+0]+p[3*faces2[4*i+1]+0])*0.33333333333333333;
p2[3*ct+1]=(2*p[3*faces2[4*i+0]+1]+p[3*faces2[4*i+1]+1])*0.33333333333333333;
p2[3*ct+2]=(2*p[3*faces2[4*i+0]+2]+p[3*faces2[4*i+1]+2])*0.33333333333333333;

p2[3*(ct+1)+0]=(p[3*faces2[4*i+0]+0]+2*p[3*faces2[4*i+1]+0])*0.33333333333333333;
p2[3*(ct+1)+1]=(p[3*faces2[4*i+0]+1]+2*p[3*faces2[4*i+1]+1])*0.33333333333333333;
p2[3*(ct+1)+2]=(p[3*faces2[4*i+0]+2]+2*p[3*faces2[4*i+1]+2])*0.33333333333333333;
ct=ct+2;
}
}

//write midtetrahedrons
ct=np;
i=0;
te2p2[faces2[4*i+2]]=ct;
te2p2[faces2[4*i+3]]=ct+1;

for (i=1;i<6*nte;i++)  {
if ((faces2[4*i+0]==faces2[4*(i-1)+0]) && (faces2[4*i+1]==faces2[4*(i-1)+1]))  {
}
else {
ct=ct+2;
}

te2p2[faces2[4*i+2]]=ct;
te2p2[faces2[4*i+3]]=ct+1;
}



const int sh2[12]={0,1,2,0,1,3,0,2,3,1,2,3};


  //sort matrix entries
  for (k=0;k<nte;k++)  {
    for (j=0;j<4;j++)  {
      i=4*k+j;
      faces[0]=te2p[sh2[3*j  ]+k*4];
      faces[1]=te2p[sh2[3*j+1]+k*4];
      faces[2]=te2p[sh2[3*j+2]+k*4];

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
      faces2[4*i+3]=20*k+16+j; //global id of triangle
    }
  }
  qsort((void*)faces2, 4*nte, 4*sizeof(int), compare_faces_geomtools); //sorting triangles in ascending order

//define face points
ct=np+2*ne;
i=0;
p2[3*ct+0]=(p[3*faces2[4*i+0]+0]+p[3*faces2[4*i+1]+0]+p[3*faces2[4*i+2]+0])*0.33333333333333333;
p2[3*ct+1]=(p[3*faces2[4*i+0]+1]+p[3*faces2[4*i+1]+1]+p[3*faces2[4*i+2]+1])*0.33333333333333333;
p2[3*ct+2]=(p[3*faces2[4*i+0]+2]+p[3*faces2[4*i+1]+2]+p[3*faces2[4*i+2]+2])*0.33333333333333333;
ct=ct+1;

for (i=1;i<4*nte;i++)  {
if ((faces2[4*i+0]==faces2[4*(i-1)+0]) && (faces2[4*i+1]==faces2[4*(i-1)+1]) && (faces2[4*i+2]==faces2[4*(i-1)+2])) {
}
else{
p2[3*ct+0]=(p[3*faces2[4*i+0]+0]+p[3*faces2[4*i+1]+0]+p[3*faces2[4*i+2]+0])*0.33333333333333333;
p2[3*ct+1]=(p[3*faces2[4*i+0]+1]+p[3*faces2[4*i+1]+1]+p[3*faces2[4*i+2]+1])*0.33333333333333333;
p2[3*ct+2]=(p[3*faces2[4*i+0]+2]+p[3*faces2[4*i+1]+2]+p[3*faces2[4*i+2]+2])*0.33333333333333333;
ct=ct+1;
}

}
ct=np+2*ne;
i=0;
te2p2[faces2[4*i+3]]=ct;


for (i=1;i<4*nte;i++)  {
if ((faces2[4*i+0]==faces2[4*(i-1)+0]) && (faces2[4*i+1]==faces2[4*(i-1)+1]) && (faces2[4*i+2]==faces2[4*(i-1)+2])) {
}
else{
ct=ct+1;
}
te2p2[faces2[4*i+3]]=ct;

}



}
