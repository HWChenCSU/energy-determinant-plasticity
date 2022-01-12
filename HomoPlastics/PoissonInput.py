######################################################
# plot weight changes after presynaptic stimulation of synapses using Poisson process
######################################################

from __future__ import division

from brian2 import *
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
nr_clst = 1   # nr of synapses per compartment   nr_clst = 10

#####################################################
# Create Neurons and Synapses
#####################################################    
synmodel = 'Chen'   # synmodel = 'Chen' , synmodel = 'Clopath', synmodel = 'nonPlast'
print('synmodel: ',synmodel)
Theta_low = morph_data['thetalow']*mV

V_rest = 0.*mV
V_thresh = 0.5*mV
eqs_in = ''' 
    dv/dt = (V_rest-v)/ms: volt
    v2 = rand()<rate_v*dt :1  (constant over dt)
    rate_v :Hz
    ds_trace/dt = -s_trace/taux :1
    '''  
N_input = NeuronGroup(nrIn*nr_clst, eqs_in, threshold='v+v2*2*V_thresh>V_thresh', 
                  reset='v=V_rest;s_trace+=x_reset*(taux/ms)',method='linear')
test_model = BRIANModel(morph)
neuron = test_model.makeNeuron_Ca(morph_data)
neuron.run_regularly('Mgblock = 1./(1.+ exp(-0.062*vu2)/3.57)',dt=defaultclock.dt)
print('Neurons created...')    
     
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
for jj in range(nrIn):
    Syn_1.connect(i=range(jj*nr_clst,(jj+1)*nr_clst),j=neuron[d_compartm[jj]:d_compartm[jj]+1])
print('Synapses created...')    

#####################################################
# Sim parameters
#####################################################   
hz_array = np.array([1.,3.,5.,10.,20.,30.,40.,50.])
#stim_time = 200*ms # stimulation time
reps = 5
ME_Ascale = 4.0   # 相当于 reps = ME_Ascale*reps = 900次脉冲

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
    print('Start compartment: '+synmodel+'_'+loc1+'_'+str(ME_Ascale)+'_'+str(nr_clst)+'_'+str(init_weight)+'_'+str(ME_A)+'_'+str(ME_Vrhigh/mV)+'_'+str(ME_Ar)+'_'+str(MEmaxRatio)+'_'+str(MEtau/second)+'_'+str(d_compartm[zzz]))
    print('> '+str(np.round(100*(zzz+1)/nrIn,2))+'%')
    weight_change1 = np.zeros(shape(hz_array))
    Erest1 = np.zeros(shape(hz_array))
    Efire1 = np.zeros(shape(hz_array))
    PE1 = np.zeros(shape(hz_array))
    MEmax1 = np.zeros(shape(hz_array))
    MEdamp1 = np.zeros(shape(hz_array))
    weight_change1_Free = np.zeros(shape(hz_array))
    Erest1_Free = np.zeros(shape(hz_array))
    Efire1_Free = np.zeros(shape(hz_array))
    PE1_Free = np.zeros(shape(hz_array))
    print('Start running ...')          
    for jj in range(size(hz_array)):
        print('-> '+str(hz_array[jj])+'Hz')
        set_init_syn(Syn_1,init_weight)
        set_init_nrn(neuron,Theta_low)
        N_input.v = V_rest
        N_input.s_trace = 0.
        N_input.rate_v = 0*Hz		
        #run(100*ms)
        ###### activate inputs
        N_input.rate_v[range(zzz*nr_clst,(zzz+1)*nr_clst)] = hz_array[jj]*Hz
        stim_time = (reps/hz_array[jj])*second
        run(stim_time) 
        ###### deactivate inputs
        N_input.rate_v = 0*Hz
        #store weight changes
        if synmodel == 'Chen':
            weight_change1[jj] = np.mean(Syn_1.wampa[zzz*nr_clst:(zzz+1)*nr_clst])
            Erest1[jj] = np.mean(Syn_1.Erest[zzz*nr_clst:(zzz+1)*nr_clst])
            Efire1[jj] = np.mean(Syn_1.Efire[zzz*nr_clst:(zzz+1)*nr_clst])	
            PE1[jj] = np.mean(Syn_1.PE[zzz*nr_clst:(zzz+1)*nr_clst])
            MEmax1[jj] = np.mean(Syn_1.MEmax[zzz*nr_clst:(zzz+1)*nr_clst])
            MEdamp1[jj] = np.mean(Syn_1.MEdamp[zzz*nr_clst:(zzz+1)*nr_clst])
            weight_change1_Free[jj] = np.mean(Syn_1.wampa_Free[zzz*nr_clst:(zzz+1)*nr_clst])
            Erest1_Free[jj] = np.mean(Syn_1.Erest_Free[zzz*nr_clst:(zzz+1)*nr_clst])
            Efire1_Free[jj] = np.mean(Syn_1.Efire_Free[zzz*nr_clst:(zzz+1)*nr_clst])	
            PE1_Free[jj] = np.mean(Syn_1.PE_Free[zzz*nr_clst:(zzz+1)*nr_clst])
        elif synmodel == 'Clopath':
#            weight_change1[jj] = np.mean(Syn_1.wampa[zzz*nr_clst:(zzz+1)*nr_clst]) + 15.*(np.mean(Syn_1.wampa[zzz*nr_clst:(zzz+1)*nr_clst])-init_weight)
            weight_change1[jj] = np.mean(Syn_1.wampa[zzz*nr_clst:(zzz+1)*nr_clst])
        else:
            weight_change1[jj] = np.mean(Syn_1.wampa[zzz*nr_clst:(zzz+1)*nr_clst])

    print('Finished running!')    
        
    titlestr = 'DataPoissonInput/'+synmodel+'_'+loc1+'_'+str(ME_Ascale)+'_'+str(nr_clst)+'_'+str(init_weight)+'_'+str(ME_A)+'_'+str(ME_Vrhigh/mV)+'_'+str(ME_Ar)+'_'+str(MEmaxRatio)+'_'+str(MEtau/second)+'_'+str(d_compartm[zzz])
    data1 = open(titlestr+'_w1.txt','w')
    json.dump(weight_change1.tolist(),data1)
    data1.close()
    data1 = open(titlestr+'_w1_Free.txt','w')
    json.dump(weight_change1_Free.tolist(),data1)
    data1.close()
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
        data1 = open(titlestr+'_Er1_Free.txt','w')
        data2 = open(titlestr+'_Ef1_Free.txt','w')
        data3 = open(titlestr+'_PE1_Free.txt','w')
        json.dump(Erest1_Free.tolist(),data1)
        json.dump(Efire1_Free.tolist(),data2)
        json.dump(PE1_Free.tolist(),data3)
        data1.close()
        data2.close()
        data3.close()