#int **LHD(int n, int k, double seed)
#{
#	int te;
#	int **LHD;
#	int *r;
#	r=new int [n];
#    LHD=new int*[k];
#	for(int iii=0;iii<k;iii++)
#	{
#		LHD[iii]=new int[n];
#	}
#     for(int j2=0;j2<n;j2++)  // first dimention is 1,2,3,.....
#	{
#		*(*(LHD+0)+j2)=(j2+1);
#//        cout<<*(*(LHD+0)+j2)<<"  ";
#	}
#//	cout<<"\n";
#	for(int j1=0;j1<(k-1);j1++)
#	{
#	   for(int cc=0;cc<n;cc++)
#	   {
#		*(r+cc)=cc+1;
#	   }
#	   for(int c=0;c<n;c++)
#	   {
#			seed=seed+j1*c;
#			seed=seed+10;
#			te=rc(n-c,seed);
#			*(*(LHD+j1+1)+c)=*(r+te);
#//		cout<<*(*(LHD+j1+1)+c)<<"  "; // print sample
#		    for(int c1=0;c1<(n-c);c1++)
#			{
#			    if(c1>=te)
#				{
#				*(r+c1)=*(r+c1+1);
#				}
#			}
#	   }//cout<<"\n";
#	}
#  return(LHD);
#   for(int iiii=0;iiii<k;iiii++)
#  {
#	delete [] LHD[iiii];
#   }
#   delete [] LHD;
#   delete [] r;
#}

from initial_design.corrd import corrd
from initial_design.rc import rc
import os
import random
import numpy


def LHD(n, k, seed, MAX):

    r = numpy.zeros(shape=(n))
    LHD = numpy.zeros(shape=(n, k))

    for j2 in range(0, n):
        LHD[j2][0] = j2+1
    for j1 in range(0, (k-1)):
        for cc in range(0, n):
            r[cc] = cc + 1
        for c in range(0, n):
            seed = seed + j1*c
            seed = seed + 10
            te = rc((n-c), seed, MAX)
            LHD[c][j1+1] = r[te]
            for c1 in range(0, (n-c)-1):
                if c1 >= te:
                    r[c1]=r[c1+1]
    return LHD

    del LHD
    del r
