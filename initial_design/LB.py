#double LB(int n,int k,int p)
#{
#	double lb=0;
#	lb=3*pow((n+1)*k,(-1))*pow(n*(n-1)*pow(2,(-1)),(pow(p,(-1))));
#	return(lb);
#}


import os
import random

def LB(n, k, p):
    lb=0

    lb = 3*(((n+1)*k)**(-1))*(((n*(n-1))*0.5)**(p**(-1)))

    return lb
