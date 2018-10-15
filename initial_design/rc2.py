#int rc2(int n,int del,double seed)
#{
#   int rctwo;
#
#   rctwo= rc( n-1, seed);
#   if (rctwo >= del)  rctwo++;
#
#   return(rctwo);
#}

from initial_design.rc import rc
import os
import random
import numpy

def rc2(n, seed, MAX, dell):

    rctwo = rc(n-1, seed, MAX)
    if (rctwo >= dell):
        rctwo = rctwo + 1

    return rctwo
