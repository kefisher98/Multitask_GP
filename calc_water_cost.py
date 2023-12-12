import os
import sys
import random
import pandas
import numpy  as np
import q_chem

path = os.path.dirname(os.path.abspath(__file__))

# tasks =====================================================

if len(sys.argv)>2:
    my_task_id = int(sys.argv[1])
    num_tasks  = int(sys.argv[2])
else:
    my_task_id = 0
    num_tasks  = 1

# water =====================================================

# path to extended xyz file
xyz_file = path + "/data/3b_data.xyz"

# path to index file
inds_file = path + "/data/water_cost_inds.csv"

# path to output file
output_file = path + "data/water_cost.out"

# split indices among processes
inds       = pandas.read_csv(inds_file)
inds       = inds.iloc[:,0].to_numpy()
my_inds    = inds[my_task_id:len(inds):num_tasks]




for num in my_inds:   
    trimer, dimers = common.get_trimer_dimer(file,num)
    common.run_dft_3b(num, 'SCAN-D3ZERO2B', trimer, dimers, output_file, be_quiet=true)
    common.run_ccsd_3b(num, trimer, dimers, output_file, be_quiet=true)

