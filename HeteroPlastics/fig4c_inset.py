######################################################
# 第9个分支及近端刺激树突神经元示意图
######################################################

from __future__ import division

from brian2 import *
import numpy as np
import matplotlib.pyplot as plt
import json
import copy as cp

import os, sys
mod_path = os.path.abspath(os.path.join('..','Model'))
sys.path.append(mod_path)

from oo_Parameters import *
from oo_equations_AMPAplast import *
from oo_initScripts import set_init_nrn, set_init_syn
from MakeNeuron_AMPAplast import *
from MorphologyData import *


start_scope()

######################################################
## Load Morpho
######################################################
#morph = '../Model/Branco2010_Morpho.swc'
#morph_data = BrancoData
#
morph = '../Model/Acker2008.swc'
morph_data = AckerData
test_model = BRIANModel(morph)
neuron = test_model.makeNeuron_Ca(morph_data)
neuron.run_regularly('Mgblock = 1./(1.+ exp(-0.062*vu2)/3.57)',dt=defaultclock.dt)

print('Neurons created...')
loc1 = 'basal'   #'tuft','apical','basal'
homoloc = 'distal'   #'proximal','distal'
abranch = 9

if loc1 == 'tuft':
    distComps = distal_Acker_tuft
    proxComps = proximal_Acker_tuft
elif loc1 == 'apical':
    distComps = distal_Acker_apical
    proxComps = proximal_Acker_apical
elif loc1 == 'basal':
    distComps = distal_Acker_basal
    proxComps = proximal_Acker_basal
else:
    print('Error!')
    sys.exit(1)
branchNr = len(proxComps)
print('branchNr: ',branchNr)

compset = list(range(proxComps[abranch],distComps[abranch]+1))
print('分支包含的分室数：',len(compset))
if homoloc == 'proximal':
    homocomps = [compset[0]]
    heterocomps = compset[1:]
else:
    homocomps = [compset[-1]]
    heterocomps = compset[:-1]

#####################################################
# Plot
#####################################################
fig, ax = plt.subplots(figsize=(1., 1.))
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
if loc1 == 'tuft':
    ax.set_xlim(-50, 310)
    ax.set_ylim(300, 500)
else: 	
    ax.set_xlim(-1, 40)
    ax.set_ylim(-6, -1)
ax.set_aspect("equal")
test_model.show_segm_3ROI(fig=fig,ax=ax,linewdbasic=0.5,segm_range1=homocomps,colorv1='m',widthv1=2,   \
    segm_range2=heterocomps,colorv2='c',widthv2=1,segm_range3=[],colorv3='b',widthv3=0)
#test_model.show_segm_byWidth(fig=fig,ax=ax,linewdbasic=0.5, idle_weight_th=idle_weight, segm_range1=compFWR,width_to_show1=wFWR,colorv1='r',   \
#    segm_range2=compBWR,width_to_show2=wBWR,colorv2='b', segm_range3=[],width_to_show3=[],colorv3='k', widthv_=0)
plt.savefig('IMG/'+'fig4c_inset'+'.eps', format='eps', dpi=1000)   #tuft,dist,nooverlap
show()
test_model.lines={}



sys.exit(1)


