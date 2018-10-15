import sys
sys.path.append("/Users/mert/python_packages")

import initial_design
from initial_design.column import column
from initial_design.corrd import corrd
from initial_design.djinformation import djinformation
from initial_design.element import element
from initial_design.LB import LB
from initial_design.LHD import LHD
from initial_design.mcorr import mcorr
from initial_design.rc import rc
from initial_design.rc2 import rc2
from initial_design.runif import runif
from initial_design.score import score
from initial_design.start_T import start_T
from initial_design.UB import UB

from pyDOE import *
from collections import defaultdict
import os
import random
import time
import numpy



#Parallel part
from joblib import Parallel, delayed
import multiprocessing




#if __name__ == "__main__":
def main(loop):

    MAX = (2**31)-1
    seed = numpy.random.rand(1)[0]
    k = 20
    n = 20
    nsearch = 2
    p = 100
    tfac = 0.95
    eps = 0.0000001
    nterp = 20 #parameter for simulated annealing
    bestyet=100000
    crit=0
    cor=0
    newcrit = 0
    newcritbest= 0
    critbest= 0
    corbest= 0
    nbestyet = 0
    itotal = 0
    temp = 0

    Mminf = numpy.zeros(shape=(2))
    xtry = numpy.zeros(shape=(n, k))
    xbest = numpy.zeros(shape=(n, k))
    x = numpy.zeros(shape=(n, k))

    t0 = start_T(n, k, p)
    print("start")




    nterp = nterp + 100
    print("nterp="+str(nterp))
    print("t0="+str(t0))

    x = LHD(n, k, seed, MAX)

    #repeated search loop
    for isearch in range(1, nsearch):
        seed=numpy.random.rand(1)[0]
        #initialize best design

        starta = time.time()
        #xbest[:][:] = x[:][:]
        #xtry[:][:] = x[:][:]
        xbest = x.copy()
        xtry = x.copy()
        #for n2 in range(0, k):
        #    for n1 in range(0, n):
        #        xbest[n1][n2] = x[n1][n2]
        #    xtry[n1][n2] = x[n1][n2]
        crit = score(x, n, k, p)
        cor = mcorr(x, n, k)
        critbest = crit
        corbest = cor
        newcritbest = (corbest**2)+(critbest-LB(n, k, p))*((UB(n, k, p)-LB(n, k, p))**(-1))

        # initialize temperatures and counts
        temp = t0
        itotal = 0
        ichange = 1

        #variable temperature loop
        while ichange == 1:
            ichange = 0
            startb = time.time()
            #constant temperature loop
            ipert = 1
            while ipert < nterp:
                itotal = itotal + 1
                #change two component in a column


                seed=numpy.random.rand(1)[0]
                ind = column(x, n, k)
                ind = ind - 1
                seed=numpy.random.rand(1)[0]
                tran1 = element(x, n, k, p)
                seed=numpy.random.rand(1)[0]
                tran2 = rc2(n, seed, MAX, tran1)

                #perturb x to xtry
                xtry[tran2][ind+1] = x[tran1][ind+1]
                xtry[tran1][ind+1] = x[tran2][ind+1]

                crittry = score(xtry, n, k, p)
                cortry = mcorr(xtry, n, k)
                newcrittry = (cortry**2)+(crittry-LB(n, k, p))*((UB(n, k, p)-LB(n, k, p))**(-1))

                #is try better than best?

                if (newcrittry < newcritbest):
                    ichange = 1



                    #for ii in range(0, k):
                    #    for jj in range(0, n):
                    #        xbest[jj][ii] = xtry[jj][ii]

                    xbest[:][:] = xtry[:][:]

                    x[tran1][ind+1] = xtry[tran1][ind+1]
                    x[tran2][ind+1] = xtry[tran2][ind+1]
                    critbest = crittry
                    newcritbest = newcrittry
                    ipert = 1
                    crit = crittry
                    cor = cortry
                    newcrit = newcrittry



                else:
                    ipert = ipert + 1

                    if  newcrittry < newcrit:
                        x[tran1][ind+1] = xtry[tran1][ind+1]
                        x[tran2][ind+1] = xtry[tran2][ind+1]
                        ichange = 1
                        crit = crittry
                        cor = cortry
                        newcrit = newcrittry

                    else:
                        delta1 = crittry - crit
                        delta2 = cortry - cor
                        rang = UB(n, k, p) - LB(n, k, p)
                        prob = numpy.exp(-delta1*(temp**(-1))-delta2*rang*(temp**(-1)))
                        seed = seed + isearch + ipert
                        q=runif(seed, MAX)
                        if prob >= q:
                            x[tran1][ind+1] = xtry[tran1][ind+1]
                            x[tran2][ind+1] = xtry[tran2][ind+1]
                            ichange = 1
                            crit = crittry
                            cor = cortry
                            newcrit = newcrittry
                        else:
                            xtry[tran1][ind+1] = x[tran1][ind+1]
                            xtry[tran2][ind+1] = x[tran2][ind+1]

            temp = temp*tfac
            endb = time.time()
            print("While loop time: "+str((endb-startb)))

        if newcritbest < (bestyet + eps):
            nbestyet = nbestyet + 1
        if newcritbest < (bestyet - eps):
            bestyet = newcritbest
            nbestyet = 1

        print(xbest)
        #for i in numpy.nditer(xbest):
        #    print(i)
        enda = time.time()
        print("isearch loop time: "+str((enda-starta)))

        combine = (mcorr(xbest, n, k)**2)+(score(xbest, n, k, p)-LB(n, k, p))*((UB(n, k, p)-LB(n, k, p))**(-1))
        Mminf = djinformation(xbest, n, k, Mminf)
        print("d="+str(Mminf[0])+" ,j="+str(Mminf[1]))
        print("score="+str(newcritbest))
        print("isearch="+str(isearch))
        print("Number of design examed="+str(itotal))
        print("phi_p="+str(score(xbest, n, k, p)))
        print("correlation="+str(mcorr(xbest, n, k))+"\n")

    del xbest
    del xtry
    del x
    del Mminf

if __name__ == "__main__":
    cycle = range(5)
    start = time.time()
    num_cores = 4
    Parallel(n_jobs=num_cores)(delayed(main)(loop) for loop in cycle)
    #main()
    end = time.time()
    print("Time: "+str((end-start)))
