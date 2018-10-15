#double mcorr(int **A, int n,int k)
#{
#	double mcor=0;
#
# 	for(int j1=1;j1<k;j1++)
#	{
#		for(int j2=0;j2<j1;j2++)
#		{
#			mcor += pow(corrd(A,n,j1,j2),2);
#		}
#	}
#
#	mcor = 2 * mcor / (double)(k*(k-1));
#	mcor = pow(mcor,(0.5));
#	return(mcor);
#}


from initial_design.corrd import corrd
import os
import random

def mcorr(A, n, k):

    mcor=0

    for i in range(1, k):
        for j in range(0, i):
            mcor = mcor + (corrd(A, n, i, j)**2)


    mcor = 2 * mcor / (k*(k-1))
    mcor = mcor ** (0.5)

    return mcor
