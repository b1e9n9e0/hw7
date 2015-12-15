/*
author:
b1e9n9e0

Date:
14.12.15

Task: implement RK ode solver 4 5
*/

#include <iostream>
#include <cmath>

using namespace std; 

//Globals
static double mu = 0.012277471;

// f(R)
void f(double *f, const double*  const R, double r, double s, const double mu){

  r=sqrt((R[0]+mu)*(R[0]+mu)+R[1]*R[1]);
  s=sqrt((R[0]-1+mu)*(R[0]-1+mu)+R[1]*R[1]);

  f[0]=R[2];
  f[1]=R[3];
  f[2]=R[0]+2*R[3]-((1-mu)*(R[0]+mu)/(r*r*r))-mu*(R[0]-1+mu)/(s*s*s);
  f[3]=R[1]-2*R[2]-((1-mu)*(R[1])/(r*r*r))-mu*R[1]/(s*s*s);
}


int main(){

double dt=1e-3;
const double T= 19.5;
//const int N= T/dt;
double dtnew;
double t=0;
double tol = 1e-5;
double r;
double s;


double Rstep[4];
double Rmax;

double R[4];
R[0] = 0.994;
R[1] = 0;
R[2] = 0;
R[3] = -2.0158510637908;

double R4[4];
double R5[4];

double k1[4];
double k2[4];
double k3[4];
double k4[4];
double k5[4];
double k6[4];
double k7[4];
double rtemp[4];


int i=0;
while(T>t){
  i++;
  cout << t <<'\t'<< R[0]<< "\t" << R[1] << "\t" << dt << endl; 
  t+=dt;

  f(k1,R,r,s,mu);

  for(int j = 0; j < 4; j++){
  	rtemp[j] = R[j] + dt/5. * k1[j];
  }

  f(k2,rtemp,r,s,mu); 


  for(int j = 0; j < 4; j++){
  	rtemp[j] = R[j] + dt * (3./40.*k1[j] + 9./40.*k2[j]);
  }
    
  f(k3,rtemp,r,s,mu);

  for(int j = 0; j < 4; j++){
  	rtemp[j] = R[j] + dt * (44./45.*k1[j] - 56./15.*k2[j] + 32./9.*k3[j]);
  }   

  f(k4,rtemp,r,s,mu); 

  for(int j = 0; j < 4; j++){
  	rtemp[j] = R[j] + dt * (19372./6561.*k1[j] - 25360./2187.*k2[j] + 64448./6561.*k3[j] - 212./729.*k4[j]);
  }   

  f(k5,rtemp,r,s,mu);

  for(int j = 0; j < 4; j++){
  	rtemp[j] = R[j] + dt * (9017./3168.*k1[j] - 355./33.*k2[j] + 46732./5247.*k3[j] + 49./176.*k4[j] - 5103./18656.*k5[j]);
  }    
  
  f(k6,rtemp,r,s,mu);

  for(int j = 0; j < 4; j++){
  	rtemp[j] = R[j] + dt * (35./384.*k1[j] + 500./113.*k3[j] + 125./192.*k4[j] - 2187./6784.*k5[j] + 11./84.*k6[j]);
  }

  f(k7,rtemp,r,s,mu);
  

  //R in 5th order

  for(int j = 0; j < 4; j++){
  	R5[j] = R[j] + dt * (35./384.*k1[j] + 500./1113.*k3[j] + 125./192.*k4[j] - 2187./6784.*k5[j] + 11./48.*k6[j]);
  }

  //R in 4th order

  for(int j = 0; j < 4; j++){
  	R4[j] = R[j] + dt * (5179./57600.*k1[j] + 7571./16695.*k3[j] + 393./640.*k4[j] - 92097./339200.*k5[j] + 187./2100.*k6[j] - 1./40.*k7[j]);
  }
  
  //determine Stepwidth

  for(int j = 0; j < 4; j++){
  	Rstep[j]=abs(R4[j]-R5[j]);
  }
  
  Rmax=Rstep[0];

  for(int j=0; j<4; j++){
    if(Rstep[j]>Rmax)
      Rmax=Rstep[j];
      else continue;
    
  }

  dtnew=dt*pow((tol/Rmax),1./5.);
  dt=dtnew;
  
  //Starting with the 4th order RK results
  
  for(int j = 0; j < 4; j++){
    R[j]=R4[j];
  }

  
   
  
}


 return 0;
}
