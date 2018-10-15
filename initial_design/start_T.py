#/////calculate starting temperature.
#	double dbar=(n+1)*k*pow(3,(-1));
#	double c=0.5;
#	double p0=0.99;
#	double pf=0.01;
#	double t1=(pow(dbar,(-p))*n*(n-1))*pow(4*c*(p-1),(-1));
#	double t2=pow((1-c),(1-p));
#	double t3=pow((1+c),(1-p));
#	double phi0=pow((t1*(t2-t3)),(pow(p,(-1))));
#	double delt0=pow((pow(phi0,p)-pow((1-c)*dbar,(-p))+pow((1-c)*dbar-1,(-p))-pow((1+c)*dbar,(-p))+pow((1+c)*dbar+1,(-p))),pow(p,(-1)));
#	t0=-delt0*pow((log(p0)),(-1));
#	//ut<<"start"<<"\n";
#	printf("start");
#
#	for(int loop=0;loop<5;loop++)
#	{
#		nterp=nterp+100;
#		out<<"nterp="<<nterp<<"\n";
#		out<<"t0="<<t0<<"\n";


import os
import random
import numpy


def start_T(n, k, p):
    
    dbar = (n+1)*k*(3**(-1))
    c = 0.5
    p0 = 0.99
    pf = 0.01
    t1 = ((dbar**(-p))*n*(n-1))*((4*c*(p-1))**(-1))
    t2 = (1-c)**(1-p)
    t3 = (1+c)**(1-p)
    phi0 = (t1*(t2-t3))**(p**(-1))
    delt0 = ((phi0**p)-(((1-c)*dbar)**(-p))+(((1-c)*dbar-1)**(-p))-(((1+c)*dbar)**(-p))+(((1+c)*dbar+1)**(-p)))**(p**(-1))
    t0=-delt0*(numpy.log(p0)**(-1))

    return t0
