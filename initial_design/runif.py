#double runif(double seed)
#{
#   double runif;
#   runif= (double) (rand()*autoseed(seed)%MAX)/MAX;
#   return(runif);
#}


import os
import random

def runif(seed, MAX):

    runif = random.randrange(0, MAX, 1) / MAX

    return runif
