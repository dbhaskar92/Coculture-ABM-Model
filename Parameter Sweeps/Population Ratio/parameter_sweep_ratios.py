import numpy as np
import scipy.io as sio

# Generate mat files to perform parameter sweep, varying population ratio
# Proliferation is turned off

## Helper functions

def idx2param(i, j, k):

    a_vals = [0.00, 0.01, 0.05, 0.09, 0.13, 0.17, 0.21, 0.25]
    
    return (a_vals[i], a_vals[j], a_vals[k])
    
    
def idx2sim(i, j, k):
    
    # i: Adhesion B-B
    # j: Adhesion B-O
    # k: Adhesion O-O

    return 1 + i + 8*j + 64*k
    
def sim2idx(n):
    
    k = int((n-1)/64)
    j = int((n - 1 - 64*k)/8)
    i = int(n - 1 - 64*k - 8*j)
    
    assert(n == idx2sim(i, j, k))
    
    return (i, j, k)

# Not included: 73, 90, 193 (orange clusters)
sim_subset = [1, 5, 36, 171, 197, 374, 500, 512]  

pop_ratios = [0.1, 0.25, 0.5, 0.75, 0.9]
num_particles = 200
pol_val = 0.005
    
task_id = 0

for sim_id in sim_subset:
    
    (i, j, k) = sim2idx(sim_id)
    (BB_adh, BO_adh, OO_adh) = idx2param(i, j, k)
    
    print(repr(sim_id) + " : (" + repr(BB_adh) + "," + repr(BO_adh) + "," + repr(OO_adh) + ")")
    
    for p_ratio in pop_ratios:
    
        task_id += 1
        
        cell_population_prop = [p_ratio, 1.0 - p_ratio]
        
        save_data = dict()
        save_data['n'] = num_particles
        
        # Color scheme
        # Red == Orange
        # Green == Blue
        
        save_data['RR'] = OO_adh
        save_data['RG'] = BO_adh
        save_data['GG'] = BB_adh
        
        save_data['pol_val'] = pol_val
        save_data['cell_pop_prop'] = cell_population_prop
        
        sio.savemat('ParamSweep_' + repr(task_id) + '.mat', save_data) 
