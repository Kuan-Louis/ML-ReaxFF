#int element(int **A, int n, int k,int p)
#{
#   int ele=0;
#	double distan=0;
#	double maxscore=0;
#	double *dist;
#	dist=new double[n];
#	int pp1, pp2, pp3;
#	int i;
#
#	for(i=0; i<n; i++)	dist[i] = 0;
#
#	for(pp1=0;pp1<n-1;pp1++)
#	{
#		for(pp2=pp1+1;pp2<n;pp2++)
#		{
#			distan = 0;
#			for(pp3=0;pp3<k;pp3++)
#			{
#				distan += abs(A[pp3][pp1]-A[pp3][pp2]);
#			}
#			dist[pp1] += pow(distan,(-p));
#			dist[pp2] += pow(distan,(-p));
#		}
#	}
#
#	for(i=0; i<n; i++)
#	{
#		if(dist[i]>maxscore)
#        {
#			maxscore = dist[i];
#			ele=i;
#		}
#	}
#
#	delete [] dist;
#	return ele;
#}


import os
import random
import numpy

def element(A, n, k, p):
    ele = 0
    distan=0
    maxscore=0
    dist = numpy.zeros(shape=(n))

    for pp1 in range(0, n):
        for pp2 in range(pp1+1, n):
            distan = 0
            for pp3 in range(0, k):
                distan = distan + abs(A[pp1][pp3] - A[pp2][pp3])
            dist[pp1] = dist[pp1] + distan**(-p)
            dist[pp2] = dist[pp2] + distan**(-p)
    for i in range(0, n):
        if (dist[i] > maxscore):
            maxscore = dist[i]
            ele = i
    del dist

    return ele
