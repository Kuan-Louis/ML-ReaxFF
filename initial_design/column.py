#int column(int **A, int n,int k)
#{
#	int col=0;
#	int i;
#	double maxcorr=0;
#
#
#	double *correlation;
#	correlation=new double[k];
#
 #	for(i=0;i<k;i++)
#	correlation[i] = 0;
#
#	for(int p1=0;p1<k;p1++)
#	{
#		for(int p2=0;p2<k;p2++)
#		{
#			if(p2!=p1)
#			{
#				correlation[p1] += corrd(A,n,p1,p2);
#			}
#		}
#
#		if (correlation[p1] < 0 ) correlation[p1] *= -1;
 #
#		if(correlation[p1] >= maxcorr)
#		{
#			maxcorr = correlation[p1];
#            col=p1;
#		}
#	}
#	delete [] correlation;
#
#	return(col);
#
#}

from initial_design.corrd import corrd
import os
import numpy

def column(A, n, k):

    col=0
    maxcorr=0
    correlation = numpy.zeros(shape=(k))

    for p1 in range(0, k):
        for p2 in range(0, k):
            if (p2!=p1):
                correlation[p1] = correlation[p1] + corrd(A, n, p1, p2)
        if correlation[p1] < 0:
            correlation[p1] = -1*correlation[p1]

        if correlation[p1] >= maxcorr:
            maxcorr = correlation[p1]
            col = p1


    del correlation

    return col
