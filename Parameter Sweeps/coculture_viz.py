import os
import math
import numpy as np
import pandas as pd
import multiprocessing
import glob 
import re
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from joblib import Parallel, delayed

def extract_number(f):
    s = f[-11:-4]
    return (int(s),f)

def get_frame_data(sim_dir, itr):
	c_iter = str(itr*50000).zfill(7)
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
	itr_s = str(itr).zfill(4)
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
			cdata_types.append('darkorange')
		elif ctype == 2:
			cdata_types.append('blue')
		else:
			cdata_types.append('black')
	plt.figure(figsize=(3,3), dpi=300)
	if edges is not None:
		for e in edges:
			plt.plot([e[0], e[1]], [e[2], e[3]], marker=None, linewidth=0.75, color='gray', alpha=0.15, zorder = 1)
	plt.scatter(np.real(position), np.imag(position), color=cdata_types, s=8, zorder = 2, edgecolors = 'k',linewidths = 0.2)
	plt.xticks([])
	plt.yticks([])
	plt.xlim([-10.05, 10.05])
	plt.ylim([-10.05, 10.05])
	#plt.margins(0,0)
	#plt.gca().xaxis.set_major_locator(plt.NullLocator())
	#plt.gca().yaxis.set_major_locator(plt.NullLocator())
	if disp:
		plt.show()
	else:
		#plt.tight_layout()
		plt.savefig(sim_dir+'_plots'+os.sep+itr_s+'.png', bbox_inches='tight')#,pad_inches = 0
	plt.close()


simlist = glob.glob("ParamSweep_*_Output")
for name in simlist:
	print(name)
	sim_dir = name
	max_iter = 100
	if(not os.path.isdir(name + '_plots')):# and os.path.isfile(name + '/Pos_5000000.dat') and os.path.isfile(name + '/Types_5000000.dat')):
		os.mkdir(sim_dir+'_plots')
	list_of_files = glob.glob(name + '/Pos_*.dat')
	file = max(list_of_files, key=extract_number)
	iternum1,_ = extract_number(file)
	iternum1=int(iternum1/50000)
	list_of_files = glob.glob(name + '/Neighbors_*.dat')
	file = max(list_of_files, key=extract_number)
	iternum2, _ = extract_number(file)
	iternum2 = int(iternum2 / 50000)
	list_of_files = glob.glob(name + '/Types_*.dat')
	file = max(list_of_files, key=extract_number)
	iternum3, _ = extract_number(file)
	iternum3= int(iternum3 / 50000)
	visualize_frame(sim_dir, min([iternum1,iternum2,iternum3]), False)
        #if(os.path.isfile(name + '/Pos_5000000.dat') and os.path.isfile(name + '/Types_5000000.dat')):
        #        if(not os.path.isfile(name + '_plots'+ os.path.sep + '0100.png')):
        #                        visualize_frame(sim_dir,max_iter,False)
                #if(not os.path.isfile(name + '_plots'+ os.path.sep + '0000.png')):
                #                        visualize_frame(sim_dir,0,False)
                #if(not os.path.isfile(name + '_plots'+ os.path.sep + '0001.png')):
                #                        visualize_frame(sim_dir,1,False)
                #if(not os.path.isfile(name + '_plots'+ os.path.sep + '0005.png')):
                #                        visualize_frame(sim_dir,5,False)
                #if(not os.path.isfile(name + '_plots'+ os.path.sep + '0010.png')):
                #                        visualize_frame(sim_dir,10,False)
                #if(not os.path.isfile(name + '_plots'+ os.path.sep + '0050.png')):
                #                        visualize_frame(sim_dir,50,False)
                #if(not os.path.isfile(name + '_plots'+ os.path.sep + '0020.png')):
                #                        visualize_frame(sim_dir,20,False)
                #if(not os.path.isfile(name + '_plots'+ os.path.sep + '0040.png')):
                #                        visualize_frame(sim_dir,40,False)
                #if(not os.path.isfile(name + '_plots'+ os.path.sep + '0060.png')):
                #                       visualize_frame(sim_dir,60,False)
                #if(not os.path.isfile(name + '_plots'+ os.path.sep + '0080.png')):
                #                        visualize_frame(sim_dir,80,False)

                        
