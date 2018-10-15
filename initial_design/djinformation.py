#void djinformation(int **A, int n, int k, int* Mminf) // try to input a LHD metrix A[][]
#{
#    const int dim= (int)(n*(n-1)*0.5);
#	double score=0;
#
#	int *d;
#	int *j;
#	d=new int[dim];
#	j=new int[dim];
#
#	for(int i=0;i<dim;i++)
#	{
#		d[i] = 0;
#		j[i] = 1;
#	}
#
#	int count=0;
#	int c=0;
#
#  	for(int k1=0;k1<(n-1);k1++)
#	{
#		for(int k2=(k1+1);k2<n;k2++)
#		{
#			d[count]=0;
#			count++;
#
#			for(int k3=0;k3<k;k3++)
#			{
#				d[count-1] += abs(*(*(A+k3)+k1)-*(*(A+k3)+k2));
#			}
#			for(int k4=0;k4<(count-1);k4++)
#			{
#				if(d[count-1]==d[k4])
#				{
#					j[k4]++;
#                    count--;
#				}
#			}
#		}
#	}
#	int Mmd=100000;
#	int Mmj=0;
#
#	for(int k5=0;k5<count;k5++)
#	{
#		if(Mmd>d[k5])
#		{
#			Mmd= d[k5];
#			Mmj= j[k5];
#		}
#	}
#
#	Mminf[0] = Mmd;
#	Mminf[1] = Mmj;
#
#	delete [] d;
#	delete [] j;
#
#}

from initial_design.corrd import corrd
import os
import random
import numpy


def djinformation(A, n, k, Mminf):

    dim= int((n*(n-1)*0.5))
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
    Mmd=100000
    Mmj=0

    for k5 in range(0, count):
        if Mmd > di[k5]:
            Mmd= di[k5]
            Mmj= jey[k5]

    Mminf[0] = Mmd
    Mminf[1] = Mmj

    return Mminf

    del di
    del jey
