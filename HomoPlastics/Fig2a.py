from __future__ import division

######################################################
# Morphology data setup in Acker&Antic2008 and Branco&Hausser 2011 ***
######################################################


from brian2 import *
import numpy as np
import matplotlib.pylab as plt
import json
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
morph = '../Model/Acker2008.swc'
morph_data = AckerData
test_model = BRIANModel(morph)

loc1 = 'basal'   #'tuft','apical','basal'
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
d_compartm = proxComps+distComps
nrIn = len(d_compartm)
#nrIn = len(proxComps)
#d_compartm = []
#for ii in range(nrIn):
#    d_compartm.append(int(0.5*(distComps[ii]+proxComps[ii])))
#####################################################
# Plot
#####################################################
fig, ax = plt.subplots(figsize=(2., 2.))
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
if loc1 == 'tuft':
    ax.set_xlim(-50, 310)
    ax.set_ylim(300, 500)
else: 	
    ax.set_xlim(-240, 110)
    ax.set_ylim(-200, 350)
ax.set_aspect("equal")
test_model.show_segm_3ROI(fig=fig,ax=ax,linewdbasic=0.5,segm_range1=d_compartm,colorv1='m',widthv1=1,segm_range2=[],colorv2='b',widthv2=2)
savefig('./IMG/'+'Fig2a.eps', format='eps', dpi=1000)
show()
test_model.lines={}

sys.exit(1)







#sn,pc,dc = test_model.calc_LayerDegree(segm_range=[655,671,707,711],layer=-1,degree=0,dist_limit=[],diam_limit=[])
#sn,pc,dc = test_model.calc_LayerDegree(segm_range=morph_data['apical'],layer=-1,degree=0,dist_limit=[300,1e4],diam_limit=[0,1])
#sn,pc,dc = test_model.calc_LayerDegree(segm_range=morph_data['apical'],layer=-1,degree=0,dist_limit=[160,1e4],diam_limit=[0,1])
#sn,pc,dc = test_model.calc_LayerDegree(segm_range=morph_data['basal'],layer=-1,degree=0,dist_limit=[],diam_limit=[])
#sn,pc,dc = test_model.calc_LayerDegree(segm_range=morph_data['apical'],layer=-1,degree=0,dist_limit=[0,150],diam_limit=[])  #dist_limit=[0,300] for AckerData
sn,pc,dc = test_model.calc_LayerDegree(segm_range=morph_data['apical'],layer=-1,degree=0,dist_limit=[],diam_limit=[])
print('Sections name = ',sn)
print('Distal compartments = ',dc)
print('Proximal compartments = ',pc)

datafile = open('../Model/'+'AckerData.txt','w')
#datafile = open('BrancoData.txt','w')
json.dump(sn,datafile)
json.dump(pc,datafile)
json.dump(dc,datafile)
datafile.close()