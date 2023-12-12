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

# organic ====================================================

# path to extended xyz file
xyz_file = path + "/data/organic_molecules.xyz"

# path to index file
inds_file = path + "/data/organic_cost_inds.csv"

# path to output file
output_file = path + "data/organic_cost.out"

# split indices among processes
inds       = pandas.read_csv(inds_file)
inds       = inds.iloc[:,0].to_numpy()
my_inds    = inds[my_task_id:len(inds):num_tasks]


for num in my_inds:
    molecule = common.get_molecule(file,num)
    common.run_dft(num, 'PBE0', molecule, output_file, be_quiet=true)
    common.run_dft(num, 'PBE', molecule, output_file, be_quiet=true)
    common.run_ccsd(num, molecule, output_file, be_quiet=true)

