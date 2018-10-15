#double corrd(int **A, int n,int a, int b)
#{
#	double cor=0;
#	double cov=0;
#	double var=0;
#	const double p2tm1 = pow(2,(-1));
#
#	for(int i1=0;i1<n;i1++)
#	{
#		cov += (A[a][i1]-(n+1)*p2tm1) * (A[b][i1]-(n+1)*p2tm1);
#		var += (A[0][i1]-(n+1)*p2tm1) * (A[0][i1]-(n+1)*p2tm1);
#	}
#	cor = cov * pow(var,(-1));
#	return(cor);
#}


import os
import random

def corrd(A, n, a, b):
    cor = 0
    cov = 0
    var = 0
    p2tm1 = 2**(-1)

    for i in range(0, n):
        cov = cov + ((A[i][a] - (n+1)*p2tm1) * (A[i][b] - (n+1)*p2tm1))
        var = var + ((A[i][0] - (n+1)*p2tm1) * (A[i][0] - (n+1)*p2tm1))
    cor = cov * (var**(-1))

    return cor
