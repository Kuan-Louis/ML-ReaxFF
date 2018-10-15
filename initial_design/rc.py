#int rc(int n,double seed ) // choose randomly from 0 to (n-1)
#{
#   int r;
#   double u; // assume I can find number from U(0,1)
#   autoseed(x) = rand()*x+x
#   u = (double) (rand()*autoseed(seed)%MAX)/MAX; //rand()%100 means in the range of 0 to 99
#   r=(int)(n*u);
#   return(r);
#}
import os
import random
import numpy

def rc(n, seed, MAX):
    u = (numpy.random.rand(1)[0])
    rc = int(n*u)

    return rc
