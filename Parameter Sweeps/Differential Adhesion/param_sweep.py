import numpy as np
import scipy.io as sio

num_particles = 200
cell_population_prop = [0.4, 0.6]
pol_val = 0.005
adh_vals = np.linspace(0.01, 0.25, 4)

task_id = 0
for RR_adh in adh_vals:
    for RG_adh in adh_vals:
        for GG_adh in adh_vals:
            task_id += 1
            save_data = dict()
            save_data['n'] = num_particles
            save_data['RR'] = RR_adh
            save_data['RG'] = RG_adh
            save_data['GG'] = GG_adh
            save_data['pol_val'] = pol_val
            save_data['cell_pop_prop'] = cell_population_prop
            sio.savemat('ParamSweep_' + repr(task_id) + '.mat', save_data)
