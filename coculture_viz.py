import os
import math
import numpy as np
import pandas as pd
import multiprocessing

import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from optparse import OptionParser
from joblib import Parallel, delayed

parser = OptionParser()
parser.add_option("-j", "--job", dest="jid", action = "store", type = "string", help = "run job", metavar = "JID")

(options, args) = parser.parse_args()
if options.jid:
  job_id = options.jid
else:
  print("Error: Job ID not specified")
  exit()

def get_frame_data(sim_dir, itr):
	c_iter = str(itr*100).zfill(7)
	posMat = pd.read_csv(os.sep.join([sim_dir, 'Pos_'+c_iter+'.dat']), header=None)
	position = [complex(pos.replace('i','j')) for pos in posMat.values[0]]
	velMat = pd.read_csv(os.sep.join([sim_dir, 'Velocity_'+c_iter+'.dat']), header=None)
	velocity = list()
	for vel in velMat.values[0]:
		if repr(vel)=='0':
			velocity.append(complex('0+0j'))
		else:
			velocity.append(complex(vel.replace('i','j')))
	nbdMat = pd.read_csv(os.sep.join([sim_dir, 'Neighbors_'+c_iter+'.dat']), header=None)
	num_neighbors = [int(nbd) for nbd in nbdMat.values]
	edges = None
	if os.path.getsize(os.sep.join([sim_dir, 'EdgeX_'+c_iter+'.dat'])) > 0:
		edgexMat = pd.read_csv(os.sep.join([sim_dir, 'EdgeX_'+c_iter+'.dat']), header=None)
		edgeyMat = pd.read_csv(os.sep.join([sim_dir, 'EdgeY_'+c_iter+'.dat']), header=None)
		edges = zip(list(edgexMat[0]),list(edgexMat[1]),list(edgeyMat[0]),list(edgeyMat[1]))
	typeMat = pd.read_csv(os.sep.join([sim_dir, 'Types_'+c_iter+'.dat']), header=None)
	cell_types = [int(ctype) for ctype in typeMat.values]
	return (position, velocity, num_neighbors, edges, cell_types)

def visualize_frame(sim_dir, itr, disp):
	itr_s = str(itr).zfill(5)
	data = get_frame_data(sim_dir, itr)
	position = data[0]
	velocity = data[1]
	num_neighbors = data[2]
	edges = data[3]
	cell_types = data[4]
	cdata_neighbors = list()
	for nbd in num_neighbors:
		if nbd == 0:
			cdata_neighbors.append('indianred')
		elif nbd < 4:
			cdata_neighbors.append('royalblue')
		else:
			cdata_neighbors.append('seagreen')
	cdata_types = list()
	for ctype in cell_types:
		if ctype == 1:
			cdata_types.append('red')
		elif ctype == 2:
			cdata_types.append('green')
		else:
			cdata_types.append('blue')
	plt.figure(figsize=(3,3), dpi=300)
	plt.scatter(np.real(position), np.imag(position), color=cdata_types, s=8)
	if edges is not None:
		for e in edges:
			plt.plot([e[0], e[1]], [e[2], e[3]], marker=None, linewidth=0.75, color='gray', alpha=0.15)
	plt.xticks([])
	plt.yticks([])
	plt.xlim([-10.1, 10.1])
	plt.ylim([-10.1, 10.1])
	if disp:
		plt.show()
	else:
		plt.savefig(sim_dir+'_plots'+os.sep+itr_s+'.png', bbox_inches='tight')
	plt.close()

sim_dir = 'ParamSweep_' + job_id + '_Output'
max_iter = 50000
os.mkdir(sim_dir + '_plots')

num_cores = multiprocessing.cpu_count()
cores_used = min(5, num_cores)
print("Number of cores: " + repr(num_cores))

Parallel(n_jobs=cores_used)(delayed(visualize_frame)(sim_dir, frame, False) for frame in range(0, max_iter+1, 500));
