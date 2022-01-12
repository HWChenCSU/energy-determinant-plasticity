######################################################
# 随机生成两个任务的突触连接分室，共生成20个随机连接
######################################################

from __future__ import division

from brian2 import *
import numpy as np
import matplotlib.pyplot as plt
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
#
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

comps = []
for ii in range(branchNr):
    comps = comps+list(range(proxComps[ii],distComps[ii]+1))
compsArr = np.array(comps)
compsLen = len(compsArr)
halfLen = int((compsLen+1)/2.0)
print('compsNr: ',compsLen)	

for ii in range(10):
    randComps = np.random.permutation(compsArr)
    compsF = randComps[:halfLen]
    compsB = randComps[halfLen:]
    print(compsF)
    print(compsB)

    titlestr = 'DataFBComps/'+loc1+'_'+str(ii)
    data1 = open(titlestr+'_compsF'+'.txt','w')
    data2 = open(titlestr+'_compsB'+'.txt','w')
    json.dump(compsF.tolist(),data1)
    json.dump(compsB.tolist(),data2)
    data1.close()
    data2.close()

#    fig, ax = plt.subplots(figsize=(2., 2.))
#    #ax.set_xlim(-240, 110)
#    #ax.set_ylim(-200, 350)
#    ax.set_xlim(-50, 310)
#    ax.set_ylim(270, 500)
#    ax.set_aspect("equal")
#    test_model.show_segm_3ROI(fig=fig,ax=ax,linewdbasic=0.5,segm_range1=[],colorv1='y',widthv1=1,segm_range2=compsF,colorv2='r',widthv2=1,segm_range3=compsB,colorv3='b',widthv3=1)
#    #surfix = './IMG/'+'Fig9ct'+'('+loc1+'_'+syn_arrange1+'_'+overlap_type+')'+'.eps'
#    #plt.savefig(surfix, format='eps', dpi=1000)   #tuft,dist,nooverlap
#    show()
#    test_model.lines={}

