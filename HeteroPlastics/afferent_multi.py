######################################################
# 多个分支、远端近端、多个频率下某个分支从近端到远端的突触权值分布数据
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

synmodel = 'Chen'   # synmodel = 'Chen' , synmodel = 'Clopath', synmodel = 'nonPlast'
print('突触可塑性模型：',synmodel)
expNr = 0   #实验编号，指两个任务对应的初始随机连接
loc1 = 'basal'   #'tuft','apical','basal'
print('树突神经元的区域：',loc1)

titlestr = 'DataFBComps/'+loc1+'_'+str(expNr)+'_'
data1 = open(titlestr+'compsF'+'.txt','r')
data2 = open(titlestr+'compsB'+'.txt','r')
compFWR = json.load(data1)
compBWR = json.load(data2)
data1.close()
data2.close()
	
nrFWR = len(compFWR)  # nr of compartments for forward running
nrBWR = len(compBWR)  # nr of compartments for backward running
#print('nrFWR=',nrFWR,', ','nrBWR=',nrBWR)
#print('compFWR=',compFWR)
#print('compBWR=',compBWR)
allcomps = compFWR + compBWR
allcompsArr = np.array(allcomps)
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
print('区域包含的分支数：',branchNr)

#homoloc = 'proximal'   #'proximal','distal'
homolocs = ['proximal', 'distal']  #'proximal','distal'
#abranch = 0
abranchs = range(branchNr)
nr_clst = 2   # nr of synapses per compartment
print('每个分室的突触数：',nr_clst)
#signal_rate = 10*Hz # activation rate of the pools
signal_rates = np.array([100])*Hz   #signal_rates = np.array([0.5,1,5,10,20,30,40,50,100])*Hz
t_stim = 100*ms # length of activation   50*ms
buffertime = 300*ms # rest time between two activations   150*ms
init_weight = 0.5  # initial weight

Theta_low = morph_data['thetalow']*mV

#####################################################
# Input Neuron
#####################################################
V_rest = 0.*mV
V_thresh = 0.5*mV

# Equations input neuron
eqs_in = ''' 
dv/dt = (V_rest-v)/ms: volt
v2 = rand()<(1.0*rate_v*dt) :1 (constant over dt)
rate_v :Hz
ds_trace/dt = -s_trace/taux :1
''' 

#####################################################
# Create neuron
#####################################################
N_input = NeuronGroup(nr_clst*(nrFWR+nrBWR), eqs_in, threshold='v+v2*2*V_thresh>V_thresh', 
                      reset='v=V_rest;s_trace+=x_reset*(taux/ms)',method='linear')#

test_model = BRIANModel(morph)
neuron = test_model.makeNeuron_Ca(morph_data)
neuron.run_regularly('Mgblock = 1./(1.+ exp(-0.062*vu2)/3.57)',dt=defaultclock.dt)

print('Neurons created...')

#####################################################
# create Synapses
#####################################################
if synmodel == 'Clopath':	
    Syn_1 = Synapses(N_input,neuron,
                    model= eq_1_plastAMPA,    
                    on_pre = eq_2_plastAMPA,    
                    method='heun'
                    )
elif synmodel == 'Chen':
    Syn_1 = Synapses(N_input,neuron,
                    model= chen_1_plastAMPA,    
                    on_pre = chen_2_plastAMPA,    
                    method='heun'
                    )
else:
    Syn_1 = Synapses(N_input,neuron,
                    model= eq_1_nonPlast,    
                    on_pre = eq_2_nonPlast,    
                    method='heun'
                    )
				
for rr in range(nrFWR):
    Syn_1.connect(i=range(rr*nr_clst,(rr+1)*nr_clst),j=neuron[compFWR[rr]:compFWR[rr]+1])
for rr in range(nrBWR):
    Syn_1.connect(i=range((rr+nrFWR)*nr_clst,(rr+nrFWR+1)*nr_clst),j=neuron[compBWR[rr]:compBWR[rr]+1])
print('Synapses created...')

print('------------------')
print('Simulating...')
for abranch in abranchs:
    print('第几个分支：',abranch)
    print('分支包含的分室数：',distComps[abranch]-proxComps[abranch]+1)
    for homoloc in homolocs:
        print('刺激位置：',homoloc)
        for signal_rate in signal_rates:
            print('刺激频率：',signal_rate)
            #####################################################
            # Initial Values
            #####################################################
            set_init_syn(Syn_1,init_weight)
            nr_syn = len(Syn_1.wampa[:])
            set_init_nrn(neuron,Theta_low)
            N_input.v = V_rest
            N_input.rate_v = 0*Hz

            compset = list(range(proxComps[abranch],distComps[abranch]+1))            
            compsind = []
            for acomp in compset:
                inputind = allcomps.index(acomp)
                compsind.append(inputind)
            if homoloc == 'proximal':
                homocomps = [compset[0]]
                heterocomps = compset[1:]
            else:
                homocomps = [compset[-1]]
                heterocomps = compset[:-1]
            homosyns = []
            for acomp in homocomps:
                inputind = allcomps.index(acomp)
                homosyns = homosyns + list(range(inputind*nr_clst,(inputind+1)*nr_clst))

            #####################################################
            # Run
            #####################################################
            #run(MEt0)
            for iii in range(40):
                N_input.rate_v[homosyns] = signal_rate
                run(t_stim)
                N_input.rate_v[homosyns] = 0*Hz
                run(buffertime)

            #####################################################
            # Weight distribution in a branch
            #####################################################
            if synmodel == 'Chen':
                wbranch = np.zeros(len(compset))
                Erbranch = np.zeros(len(compset))
                Efbranch = np.zeros(len(compset))
                PEbranch = np.zeros(len(compset))
                MEmaxbranch = np.zeros(len(compset))
                MEdampbranch = np.zeros(len(compset))
                wbranch_Free = np.zeros(len(compset))
                Erbranch_Free = np.zeros(len(compset))
                Efbranch_Free = np.zeros(len(compset))
                PEbranch_Free = np.zeros(len(compset))
                for jj in range(len(compset)):
                    compind = compsind[jj]
                    wbranch[jj] = np.mean(Syn_1.wampa[compind*nr_clst:(compind+1)*nr_clst])
                    Erbranch[jj] = np.mean(Syn_1.Erest[compind*nr_clst:(compind+1)*nr_clst])
                    Efbranch[jj] = np.mean(Syn_1.Efire[compind*nr_clst:(compind+1)*nr_clst])
                    PEbranch[jj] = np.mean(Syn_1.PE[compind*nr_clst:(compind+1)*nr_clst])
                    MEmaxbranch[jj] = np.mean(Syn_1.MEmax[compind*nr_clst:(compind+1)*nr_clst])
                    MEdampbranch[jj] = np.mean(Syn_1.MEdamp[compind*nr_clst:(compind+1)*nr_clst])
                    wbranch_Free[jj] = np.mean(Syn_1.wampa_Free[compind*nr_clst:(compind+1)*nr_clst])
                    Erbranch_Free[jj] = np.mean(Syn_1.Erest_Free[compind*nr_clst:(compind+1)*nr_clst])
                    Efbranch_Free[jj] = np.mean(Syn_1.Efire_Free[compind*nr_clst:(compind+1)*nr_clst])
                    PEbranch_Free[jj] = np.mean(Syn_1.PE_Free[compind*nr_clst:(compind+1)*nr_clst])
                if homoloc == 'distal':
                    wbranch = wbranch[::-1]
                    Erbranch = Erbranch[::-1]
                    Efbranch = Efbranch[::-1]
                    PEbranch = PEbranch[::-1]
                    MEmaxbranch = MEmaxbranch[::-1]
                    MEdampbranch = MEdampbranch[::-1]
                    wbranch_Free = wbranch_Free[::-1]
                    Erbranch_Free = Erbranch_Free[::-1]
                    Efbranch_Free = Efbranch_Free[::-1]
                    PEbranch_Free = PEbranch_Free[::-1]
            elif synmodel == 'Clopath':
                wbranch = np.zeros(len(compset))
                for jj in range(len(compset)):
                    compind = compsind[jj]
#                    wmean = np.mean(Syn_1.wampa[compind*nr_clst:(compind+1)*nr_clst])
#                    wbranch[jj] = wmean + 15.*(wmean-init_weight)
                    wbranch[jj] = np.mean(Syn_1.wampa[compind*nr_clst:(compind+1)*nr_clst])
                if homoloc == 'distal':
                    wbranch = wbranch[::-1]
            else:
                wbranch = np.zeros(len(compset))
                for jj in range(len(compset)):
                    compind = compsind[jj]
                    wbranch[jj] = np.mean(Syn_1.wampa[compind*nr_clst:(compind+1)*nr_clst])
                if homoloc == 'distal':
                    wbranch = wbranch[::-1]	
				
            #------------------------------	
            titlestr = 'Data/'+synmodel+'_'+loc1+'_'+str(nr_clst)+'_'+str(init_weight)+'_'+str(abranch)+'_'+homoloc+'_'+str(signal_rate/Hz)
            data0 = open(titlestr+'_wbranch.txt','w')		
            json.dump(wbranch.tolist(),data0)
            data0.close()
            if synmodel == 'Chen':
                data0 = open(titlestr+'_wbranch_Free.txt','w')		
                json.dump(wbranch_Free.tolist(),data0)
                data0.close()
                data0 = open(titlestr+'_Erbranch.txt','w')		
                json.dump(Erbranch.tolist(),data0)
                data0.close()
                data0 = open(titlestr+'_Efbranch.txt','w')		
                json.dump(Efbranch.tolist(),data0)
                data0.close()
                data0 = open(titlestr+'_PEbranch.txt','w')		
                json.dump(PEbranch.tolist(),data0)
                data0.close()
                data0 = open(titlestr+'_MEmaxbranch.txt','w')		
                json.dump(MEmaxbranch.tolist(),data0)
                data0.close()
                data0 = open(titlestr+'_MEdampbranch.txt','w')		
                json.dump(MEdampbranch.tolist(),data0)
                data0.close()				
                data0 = open(titlestr+'_Erbranch_Free.txt','w')		
                json.dump(Erbranch_Free.tolist(),data0)
                data0.close()
                data0 = open(titlestr+'_Efbranch_Free.txt','w')		
                json.dump(Efbranch_Free.tolist(),data0)
                data0.close()
                data0 = open(titlestr+'_PEbranch_Free.txt','w')		
                json.dump(PEbranch_Free.tolist(),data0)
                data0.close()				

sys.exit(1)


