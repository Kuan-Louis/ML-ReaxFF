import sys
sys.path.append("/storage/home/mys12/python_packages")

import ffield_generator
from ffield_generator.parameter import assign_parameter
from ffield_generator.params_input_train import params_input_train
from ffield_generator.create_folder import create_folder
from ffield_generator.params_output import params_output
from ffield_generator.error_output import error_output
from ffield_generator.parameter_generator import parameter_generator
from ffield_generator.generate_training_folder import generate_training_folder
from ffield_generator.remove_files import remove_files
from ffield_generator.geometry_check import geometry_check
from ffield_generator.check_iappen import check_iappen


from pyDOE import *
from collections import defaultdict
import os
import random

#Parallel part
from joblib import Parallel, delayed
import multiprocessing


cycle = range(500)
parameters = params_input_train()
values = parameter_generator(parameters, len(cycle))
length = len(parameters)


def main(j):
    
    #taking the information about which parameters will be altered

    try:
           
        folder_name = 'ffield-'+ str(j)
        generate_training_folder(folder_name)
        os.chdir(folder_name)
       # check_iappen()
        
        for i in range(0,length):

            value = str("{:.4f}".format(values[j][i]))            
            assign_parameter(parameters[i][0], parameters[i][1], parameters[i][2], value)
            params_output(parameters[i][0],parameters[i][1],parameters[i][2], value)
            
        os.system('./exe')
        

        error_output(length)

        
        os.remove('ffieldss')
        os.remove('control')
        os.remove('dipole.out')
        os.remove('exe')
        os.remove('output.*')
        os.remove('geo')
        os.remove('iopt')
        os.remove('mol*')
        os.remove('summary.txt')
        os.remove('run.log')
        os.remove('fort.13')
        os.remove('fort.20')
        os.remove('fort.35')
        os.remove('fort.71')
	    os.remove('fort.9')
	    os.remove('fort.91')
	    os.remove('fort.98')


        os.chdir('..')
    
    except:
        pass

#num_cores = multiprocessing.cpu_count()
num_cores = 20
Parallel(n_jobs=num_cores)(delayed(main)(j) for j in cycle)
