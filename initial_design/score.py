#double score(int **A, int n, int k,int p) // try to input a LHD metrix A[][]
#{
#//	int A[3][6];
#
#	//int n=6;  // number of points
#	//int k=3;  // number of dimention
#	int i;
#    int dim=n*(n-1)/2;
#	double score=0;
#	int count=0;
#	int c=0;
#
#	int *d;
#	int *j;
#	d=new int[dim];
#	j=new int[dim];
#
#	for(i=0;i<dim;i++)
#	{
#		d[i] = 0;
#		j[i] = 1;
#	}
#	count= c = 0;
# 	for(int k1=0;k1<(n-1);k1++)
#	{
#		for(int k2=(k1+1);k2<n;k2++)
#		{
#			d[count++] = 0;
#		//	cout<< "count="<<count<<"\n ";
#			for(int k3=0;k3<k;k3++)
#			{
#				d[count-1] += abs(A[k3][k1]-A[k3][k2]);
#			}
#			for(int k4=0;k4<(count-1);k4++)
#			{
#				if(d[count-1]== d[k4])
#				{
#					j[k4]++;
#                    count--;
#				}
#			}
#		}
#
#	}
#	for(i=0; i<count; i++)
#	{
#		score += (j[i])* pow(d[i],(-p));
#	}
#
#	score=pow(score,(pow(p,(-1))));
#
#	delete [] d;
#	delete [] j;
#
#	return(score);
#}

from initial_design.corrd import corrd
import os
import random
import numpy


def score(A, n, k, p):

    dim=int(n*(n-1)/2)
    score=0
    count=0
    c=0

    di = numpy.zeros(shape=(dim))
    jey = numpy.ones(shape=(dim))


    for k1 in range(0, (n-1)):
        for k2 in range((k1+1), n):
            di[count]=0
            count = count + 1

            for k3 in range(0, k):
                di[count-1] = di[count-1] + abs(A[k1][k3]-A[k2][k3])
            for k4 in range(0, (count-1)):
                if di[count-1] == di[k4]:
                    jey[k4] = jey[k4] + 1
                    count = count - 1

    for i in range(0, count):
        score = score + (jey[i])* (di[i]**(-p))

    score=(score**((p**(-1))))


    del di
    del jey

    return score
