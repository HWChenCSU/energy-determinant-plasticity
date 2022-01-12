######################################################
# Reproduce the rate dependence of plasticity as observed in Sjostrom et al. 2001
######################################################

from __future__ import division

from brian2 import *
#prefs.codegen.target = "numpy"
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
loc1 = 'basal'   #'tuft','apical','basal'
print('loc1: ',loc1)
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
# Create Neurons and Synapses
#####################################################    
synmodel = 'Chen'   # synmodel = 'Chen' , synmodel = 'Clopath', synmodel = 'nonPlast'
print('synmodel: ',synmodel)
Theta_low = morph_data['thetalow']*mV

V_rest = -70.*mV
tau_in = 8.*ms
V_thresh = -45.*mV
C = 200.*pF # membrane capacitance   
eqs_in = ''' 
    dv/dt = (V_rest-v)/tau_in + Idrive/C: volt
    Idrive : amp
    ds_trace/dt = -s_trace/taux :1
    '''     
N_input = NeuronGroup(2*nrIn, eqs_in, threshold='v>V_thresh', 
                      reset='v=V_rest;s_trace+=x_reset*(taux/ms)', method='linear')#    
test_model = BRIANModel(morph)
neuron = test_model.makeNeuron_Ca(morph_data)
neuron.run_regularly('Mgblock = 1./(1.+ exp(-0.062*vu2)/3.57)',dt=defaultclock.dt)
neuron2 = test_model.makeNeuron_Ca(morph_data) 
neuron2.run_regularly('Mgblock = 1./(1.+ exp(-0.062*vu2)/3.57)',dt=defaultclock.dt)
print('Neurons created...')    
     
if synmodel == 'Clopath':	
    Syn_1 = Synapses(N_input,neuron,
                    model= eq_1_plastAMPA,    
                    on_pre = eq_2_plastAMPA,    
                    method='heun'
                    )
    Syn_2 = Synapses(N_input,neuron2,
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
    Syn_2 = Synapses(N_input,neuron2,
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
    Syn_2 = Synapses(N_input,neuron2,
                    model= eq_1_nonPlast,    
                    on_pre = eq_2_nonPlast,    
                    method='heun'
                    )
for jj in range(nrIn):
    Syn_1.connect(i=jj,j=neuron[d_compartm[jj]:d_compartm[jj]+1])
    Syn_2.connect(i=nrIn+jj,j=neuron2[d_compartm[jj]:d_compartm[jj]+1])
print('Synapses created...')    

#####################################################
# Sim parameters
#####################################################   
# rates for the protocol as in sjostrom2001
hz_array =  np.array([0.1,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.])     #np.array([0.1,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.])
# number of pairings
reps = 5   #5
ME_Ascale = 60.0/5.0   #相当于 reps = ME_Ascale*reps = 60

init_weight = 0.5
ME_A = 0.02
ME_Vrhigh = -60*mV
ME_Ar = 0.2
MEmaxRatio = 175.0
MEtau = 2.0*second

MESmax0 = 25.0
#####################################################
# Simulation
#####################################################   
for zzz in range(nrIn):	#range(nrIn)
    print('Start compartment: '+synmodel+'_'+loc1+'_'+str(init_weight)+'_'+str(ME_A)+'_'+str(ME_Vrhigh/mV)+'_'+str(ME_Ar)+'_'+str(MEmaxRatio)+'_'+str(MEtau/second)+'_'+str(d_compartm[zzz]))
    print('> '+str(np.round(100*(zzz+1)/nrIn,2))+'%')
    weight_change1 = np.zeros(shape(hz_array))
    weight_change2 = np.zeros(shape(hz_array))
    Erest1 = np.zeros(shape(hz_array))
    Efire1 = np.zeros(shape(hz_array))
    PE1 = np.zeros(shape(hz_array))
    MEmax1 = np.zeros(shape(hz_array))
    MEdamp1 = np.zeros(shape(hz_array))
    Erest2 = np.zeros(shape(hz_array))
    Efire2 = np.zeros(shape(hz_array))
    PE2 = np.zeros(shape(hz_array))
    MEmax2 = np.zeros(shape(hz_array))
    MEdamp2 = np.zeros(shape(hz_array))	
    print('Start running ...')          
    for jj in range(size(hz_array)):
        pair_interval = 1000./hz_array[jj]-13.
        print('-> '+str(hz_array[jj])+'Hz')
        set_init_syn(Syn_1,init_weight)
        set_init_syn(Syn_2,init_weight)               
        set_init_nrn(neuron,Theta_low)
        set_init_nrn(neuron2,Theta_low)
        N_input.v = V_rest
        N_input.s_trace = 0.
		
        #run(100*ms)      
        # Pairings
        for ii in range(reps):
            neuron.I = 0.*pA
            neuron2.I = 0.*pA
            N_input.Idrive = 0.*mA
            ###### 1st SPIKE
            neuron2.main.I = 1000.*pA
            N_input.Idrive[zzz] = 2000.*pA
            run(3*ms) 
            neuron2.I = 0.*pA
            N_input.Idrive = 0.*mA
            run(7*ms) 
            ###### 2nd SPIKE
            neuron.main.I = 1000.*pA    
            N_input.Idrive[nrIn+zzz] = 2000.*pA
            run(3*ms) 
            neuron.I = 0.*pA
            N_input.Idrive = 0.*mA
            ######
            run(pair_interval*ms) 

        #store weight changes
        if synmodel == 'Chen':
            weight_change1[jj] = Syn_1.wampa[zzz]
            weight_change2[jj] = Syn_2.wampa[zzz]			
            Erest1[jj] = Syn_1.Erest[zzz]
            Efire1[jj] = Syn_1.Efire[zzz]
            PE1[jj] = Syn_1.PE[zzz]
            MEmax1[jj] = Syn_1.MEmax[zzz]
            MEdamp1[jj] = Syn_1.MEdamp[zzz] 
            Erest2[jj] = Syn_2.Erest[zzz]
            Efire2[jj] = Syn_2.Efire[zzz]
            PE2[jj] = Syn_2.PE[zzz]
            MEmax2[jj] = Syn_2.MEmax[zzz]
            MEdamp2[jj] = Syn_2.MEdamp[zzz]
        elif synmodel == 'Clopath':
            weight_change1[jj] = Syn_1.wampa[zzz] + 15.*(Syn_1.wampa[zzz]-init_weight)
            weight_change2[jj] = Syn_2.wampa[zzz] + 15.*(Syn_2.wampa[zzz]-init_weight)
        else:
            weight_change1[jj] = Syn_1.wampa[zzz]
            weight_change2[jj] = Syn_2.wampa[zzz]

    print('Finished running!')    

    titlestr = 'DataRateDependence/'+synmodel+'_'+loc1+'_'+str(init_weight)+'_'+str(ME_A)+'_'+str(ME_Vrhigh/mV)+'_'+str(ME_Ar)+'_'+str(MEmaxRatio)+'_'+str(MEtau/second)+'_'+str(d_compartm[zzz]) 
    data1 = open(titlestr+'_w1.txt','w')
    data2 = open(titlestr+'_w2.txt','w')
    json.dump(weight_change1.tolist(),data1)
    json.dump(weight_change2.tolist(),data2)
    data1.close()
    data2.close()
    if synmodel == 'Chen':	
        data1 = open(titlestr+'_Er1.txt','w')
        data2 = open(titlestr+'_Ef1.txt','w')
        data3 = open(titlestr+'_PE1.txt','w')
        data4 = open(titlestr+'_MEmax1.txt','w')
        data5 = open(titlestr+'_MEdamp1.txt','w')
        json.dump(Erest1.tolist(),data1)
        json.dump(Efire1.tolist(),data2)
        json.dump(PE1.tolist(),data3)
        json.dump(MEmax1.tolist(),data4)
        json.dump(MEdamp1.tolist(),data5)
        data1.close()
        data2.close()
        data3.close()
        data4.close()
        data5.close()			
        data1 = open(titlestr+'_Er2.txt','w')
        data2 = open(titlestr+'_Ef2.txt','w')
        data3 = open(titlestr+'_PE2.txt','w')
        data4 = open(titlestr+'_MEmax2.txt','w')
        data5 = open(titlestr+'_MEdamp2.txt','w')
        json.dump(Erest2.tolist(),data1)
        json.dump(Efire2.tolist(),data2)
        json.dump(PE2.tolist(),data3)
        json.dump(MEmax2.tolist(),data4)
        json.dump(MEdamp2.tolist(),data5)
        data1.close()
        data2.close()
        data3.close()
        data4.close()
        data5.close()					
	
