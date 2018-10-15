#double UB(int n, int k,int p)
#{
#	double ub=0;
#	for(int i=1;i<(n-1);i++)
#	{
#		ub += (n-i)* pow((i*k),(-p));
#	}
#	ub = pow(ub,(pow(p,(-1))));
#	return(ub);
#}



import os
import random

def UB(n, k, p):
    ub=0

    for i in range(1, (n-1)):
        ub = ub + (n-i)*((i*k)**(-p))

    ub = ub**(p**(-1))

    return ub
