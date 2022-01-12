######################################################
# basal区域所有分支、远近端、不同频率刺激下，Chen可塑性模型运算的权值和势能分布图形
######################################################

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import json

import os, sys
mod_path = os.path.abspath(os.path.join('..','Model'))
sys.path.append(mod_path)
from oo_Parameters import *
from MorphologyData import *

loc1 = 'basal'   #'tuft','apical','basal'
print('树突神经元的区域：',loc1)
signal_rate = 100.0 # activation rate of the pools
print('刺激频率：',signal_rate)
init_weight = 0.5  # initial weight
nr_clst = 2   # nr of synapses per compartment

morph = '../Model/Acker2008.swc'
morph_data = AckerData
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

homolocs = ['proximal', 'distal']  #'proximal','distal'
abranchs = range(branchNr)

#------------------------------	
#w_exp_mean = np.array([0.18,-0.1,-0.08,-0.06,-0.04,-0.06,-0.08])   #实验数据来自[301_36_5]Fig3Ac，其中0,1,2取平均值对应横坐标0
w_exp_mean = np.array([0.32,-0.1,-0.08,-0.06,-0.04,-0.06,-0.08])   #实验数据来自[301_36_5]Fig3Ac，其中0,1,2取最大值对应横坐标0
w_exp_std = np.array([0.08,0.08,0.01,0.04,0.02,0.04,0.05])
PE_Chen_mean = np.zeros(len(w_exp_mean))
PE_Chen_std = np.zeros(len(w_exp_mean))
MEmax_Chen_mean = np.zeros(len(w_exp_mean))
MEmax_Chen_std = np.zeros(len(w_exp_mean))
MEdamp_Chen_mean = np.zeros(len(w_exp_mean))
MEdamp_Chen_std = np.zeros(len(w_exp_mean))
Er_Chen_mean = np.zeros(len(w_exp_mean))
Er_Chen_std = np.zeros(len(w_exp_mean))
Ef_Chen_mean = np.zeros(len(w_exp_mean))
Ef_Chen_std = np.zeros(len(w_exp_mean))
w_Chen_mean = np.zeros(len(w_exp_mean))
w_Chen_std = np.zeros(len(w_exp_mean))
 
groupNr = len(w_exp_mean)-1
listChenPE = [[] for _ in range(len(w_exp_mean))]
listChenMEmax = [[] for _ in range(len(w_exp_mean))]
listChenMEdamp = [[] for _ in range(len(w_exp_mean))]
listChenEr = [[] for _ in range(len(w_exp_mean))]
listChenEf = [[] for _ in range(len(w_exp_mean))]
listChenW = [[] for _ in range(len(w_exp_mean))]
synmodel = 'Chen'   # synmodel = 'Chen' , synmodel = 'Clopath', synmodel = 'nonPlast'
print('突触可塑性模型：',synmodel)
for abranch in abranchs:
    for homoloc in homolocs:
        titlestr = 'Data/'+synmodel+'_'+loc1+'_'+str(nr_clst)+'_'+str(init_weight)+'_'+str(abranch)+'_'+homoloc+'_'+str(signal_rate)
        data0 = open(titlestr+'_PEbranch.txt','r')
        PEbranch = json.load(data0)
        data0.close()
        data0 = open(titlestr+'_MEmaxbranch.txt','r')
        MEmaxbranch = json.load(data0)
        data0.close()
        data0 = open(titlestr+'_MEdampbranch.txt','r')
        MEdampbranch = json.load(data0)
        data0.close()
        data0 = open(titlestr+'_Erbranch.txt','r')
        Erbranch = json.load(data0)
        data0.close()
        data0 = open(titlestr+'_Efbranch.txt','r')
        Efbranch = json.load(data0)
        data0.close()
        data0 = open(titlestr+'_wbranch.txt','r')
        wbranch = json.load(data0)
        data0.close()
        listChenPE[0] = listChenPE[0]+[PEbranch[0]]
        listChenMEmax[0] = listChenMEmax[0]+[MEmaxbranch[0]]
        listChenMEdamp[0] = listChenMEdamp[0]+[MEdampbranch[0]]
        listChenEr[0] = listChenEr[0]+[Erbranch[0]]
        listChenEf[0] = listChenEf[0]+[Efbranch[0]]
        listChenW[0] = listChenW[0]+[wbranch[0]]
        groupInt = int((len(wbranch)-1)/groupNr)
#        groupInt = int(round((len(wbranch)-1)/groupNr))
        for ii in range(groupNr):
            beginInd = ii*groupInt+1
            endInd = beginInd+groupInt
            if ii==groupNr-1:
                listChenPE[ii+1] = listChenPE[ii+1]+PEbranch[beginInd:]
                listChenMEmax[ii+1] = listChenMEmax[ii+1]+MEmaxbranch[beginInd:]
                listChenMEdamp[ii+1] = listChenMEdamp[ii+1]+MEdampbranch[beginInd:]
                listChenEr[ii+1] = listChenEr[ii+1]+Erbranch[beginInd:]
                listChenEf[ii+1] = listChenEf[ii+1]+Efbranch[beginInd:]
                listChenW[ii+1] = listChenW[ii+1]+wbranch[beginInd:]
            else:
                listChenPE[ii+1] = listChenPE[ii+1]+PEbranch[beginInd:endInd]
                listChenMEmax[ii+1] = listChenMEmax[ii+1]+MEmaxbranch[beginInd:endInd]
                listChenMEdamp[ii+1] = listChenMEdamp[ii+1]+MEdampbranch[beginInd:endInd]
                listChenEr[ii+1] = listChenEr[ii+1]+Erbranch[beginInd:endInd]
                listChenEf[ii+1] = listChenEf[ii+1]+Efbranch[beginInd:endInd]
                listChenW[ii+1] = listChenW[ii+1]+wbranch[beginInd:endInd]

for jj in range(len(w_exp_mean)):
    PE_Chen_mean[jj] = np.mean(listChenPE[jj])
    PE_Chen_std[jj] = np.std(listChenPE[jj])    #/np.sqrt(len(listChenPE[jj]))
    MEmax_Chen_mean[jj] = np.mean(listChenMEmax[jj])
    MEmax_Chen_std[jj] = np.std(listChenMEmax[jj])    #/np.sqrt(len(listChenMEmax[jj]))
    MEdamp_Chen_mean[jj] = np.mean(listChenMEdamp[jj])
    MEdamp_Chen_std[jj] = np.std(listChenMEdamp[jj])    #/np.sqrt(len(listChenMEdamp[jj]))
    Er_Chen_mean[jj] = np.mean(listChenEr[jj])
    Er_Chen_std[jj] = np.std(listChenEr[jj])    #/np.sqrt(len(listChenEr[jj]))
    Ef_Chen_mean[jj] = np.mean(listChenEf[jj])
    Ef_Chen_std[jj] = np.std(listChenEf[jj])    #/np.sqrt(len(listChenEf[jj]))		
    w_Chen_mean[jj] = np.mean(listChenW[jj])-init_weight
    w_Chen_std[jj] = np.std(listChenW[jj])    #/np.sqrt(len(listChenW[jj]))

#------------------------------	
PE_Chen_mean_Free = np.zeros(len(w_exp_mean))
PE_Chen_std_Free = np.zeros(len(w_exp_mean))
Er_Chen_mean_Free = np.zeros(len(w_exp_mean))
Er_Chen_std_Free = np.zeros(len(w_exp_mean))
Ef_Chen_mean_Free = np.zeros(len(w_exp_mean))
Ef_Chen_std_Free = np.zeros(len(w_exp_mean))
w_Chen_mean_Free = np.zeros(len(w_exp_mean))
w_Chen_std_Free = np.zeros(len(w_exp_mean))
 
groupNr = len(w_exp_mean)-1
listChenPE_Free = [[] for _ in range(len(w_exp_mean))]
listChenEr_Free = [[] for _ in range(len(w_exp_mean))]
listChenEf_Free = [[] for _ in range(len(w_exp_mean))]
listChenW_Free = [[] for _ in range(len(w_exp_mean))]
synmodel = 'Chen'   # synmodel = 'Chen' , synmodel = 'Clopath', synmodel = 'nonPlast'
print('突触可塑性模型：',synmodel)
for abranch in abranchs:
    for homoloc in homolocs:
        titlestr = 'Data/'+synmodel+'_'+loc1+'_'+str(nr_clst)+'_'+str(init_weight)+'_'+str(abranch)+'_'+homoloc+'_'+str(signal_rate)
        data0 = open(titlestr+'_PEbranch_Free.txt','r')
        PEbranch_Free = json.load(data0)
        data0.close()
        data0 = open(titlestr+'_Erbranch_Free.txt','r')
        Erbranch_Free = json.load(data0)
        data0.close()
        data0 = open(titlestr+'_Efbranch_Free.txt','r')
        Efbranch_Free = json.load(data0)
        data0.close()
        data0 = open(titlestr+'_wbranch_Free.txt','r')
        wbranch_Free = json.load(data0)
        data0.close()
        listChenPE_Free[0] = listChenPE_Free[0]+[PEbranch_Free[0]]
        listChenEr_Free[0] = listChenEr_Free[0]+[Erbranch_Free[0]]
        listChenEf_Free[0] = listChenEf_Free[0]+[Efbranch_Free[0]]
        listChenW_Free[0] = listChenW_Free[0]+[wbranch_Free[0]]
        groupInt = int((len(wbranch_Free)-1)/groupNr)
#        groupInt = int(round((len(wbranch)-1)/groupNr))
        for ii in range(groupNr):
            beginInd = ii*groupInt+1
            endInd = beginInd+groupInt
            if ii==groupNr-1:
                listChenPE_Free[ii+1] = listChenPE_Free[ii+1]+PEbranch_Free[beginInd:]
                listChenEr_Free[ii+1] = listChenEr_Free[ii+1]+Erbranch_Free[beginInd:]
                listChenEf_Free[ii+1] = listChenEf_Free[ii+1]+Efbranch_Free[beginInd:]
                listChenW_Free[ii+1] = listChenW_Free[ii+1]+wbranch_Free[beginInd:]
            else:
                listChenPE_Free[ii+1] = listChenPE_Free[ii+1]+PEbranch_Free[beginInd:endInd]
                listChenEr_Free[ii+1] = listChenEr_Free[ii+1]+Erbranch_Free[beginInd:endInd]
                listChenEf_Free[ii+1] = listChenEf_Free[ii+1]+Efbranch_Free[beginInd:endInd]
                listChenW_Free[ii+1] = listChenW_Free[ii+1]+wbranch_Free[beginInd:endInd]

for jj in range(len(w_exp_mean)):
    PE_Chen_mean_Free[jj] = np.mean(listChenPE_Free[jj])
    PE_Chen_std_Free[jj] = np.std(listChenPE_Free[jj])    #/np.sqrt(len(listChenPE_Free[jj]))
    Er_Chen_mean_Free[jj] = np.mean(listChenEr_Free[jj])
    Er_Chen_std_Free[jj] = np.std(listChenEr_Free[jj])    #/np.sqrt(len(listChenEr_Free[jj]))
    Ef_Chen_mean_Free[jj] = np.mean(listChenEf_Free[jj])
    Ef_Chen_std_Free[jj] = np.std(listChenEf_Free[jj])    #/np.sqrt(len(listChenEf_Free[jj]))		
    w_Chen_mean_Free[jj] = np.mean(listChenW_Free[jj])-init_weight
    w_Chen_std_Free[jj] = np.std(listChenW_Free[jj])    #/np.sqrt(len(listChenW_Free[jj]))


	
#------------------------------					
#fig, ax = plt.subplots(1, 1, figsize=(1.4, 1.0))
#fig.subplots_adjust(left=0.2, right=0.9, bottom=0.15, top=0.95)
fig, ax = plt.subplots(1, 1, figsize=(2.0, 1.4))
fig.subplots_adjust(left=0.3, right=0.9, bottom=0.3, top=0.9)
plt.sca(ax)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both',direction='out',length=2,labelcolor=(0.5,0.5,0.5))
#ax.set_prop_cycle(color=['red', 'olive', 'blue', 'green', 'black', 'brown', 'orange', 'skyblue', 'springgreen', 'gray', 'darkblue', 'magenta'])
#ax.set_xlim(-0.05, 1.05)
ax.set_ylim(-30, 30)
plt.xticks(np.arange(len(w_exp_mean)),fontsize=6,fontname='Times New Roman')   #[0,1,2,3,4,5,6],
plt.yticks(fontsize=6,fontname='Times New Roman')   #[-100,-50,0,50,100],
PEmax = sign(PE_Chen_mean)*(MEdamp_Chen_mean*MEmax_Chen_mean+MESmax0)
ax.plot(np.arange(len(w_exp_mean)),PEmax,'--r',linewidth=1.,label='$\itP_{max}$')
ax.plot(0,PEmax[0],'mo',ms=2)
ax.plot(range(1,len(w_exp_mean)),PEmax[1:],'co',ms=2)
ax.plot(np.arange(len(w_exp_mean)),PE_Chen_mean,'-k',lw=1.,label='$\itP$')
ax.plot(0,PE_Chen_mean[0],'mo',ms=2)   #,mfc='w'
ax.plot(range(1,len(w_exp_mean)),PE_Chen_mean[1:],'co',ms=2)   #,mfc='w'
#ax.fill_between(np.arange(len(w_exp_mean)),PE_Chen_mean+PE_Chen_std,PE_Chen_mean-PE_Chen_std,color=(0.95,0.95,0.95),alpha=0.1)   #color=(1.0,0.9,1.0)
ax.plot(np.arange(len(w_exp_mean)),PE_Chen_mean_Free,'-b',lw=1.)    #,label='$\itP$'
ax.plot(0,PE_Chen_mean_Free[0],'mo',ms=2)   #,mfc='w'
ax.plot(range(1,len(w_exp_mean)),PE_Chen_mean_Free[1:],'co',ms=2)   #,mfc='w'
ax.fill_between(np.arange(len(w_exp_mean)),PE_Chen_mean_Free+PE_Chen_std_Free,PE_Chen_mean_Free-PE_Chen_std_Free,color=(0.95,0.95,1.0),alpha=0.1)   #color=(1.0,0.9,1.0)
#ax.legend(loc='lower right',fontsize=4,frameon=False)
#ax.legend(loc='lower left',prop={'family':'Times New Roman','size':4},frameon=True,framealpha=0.1)
#plt.xlabel('Distance from stimulation site',fontdict={'family': 'Times New Roman','size':8})
#plt.ylabel('Potential energy\n(fJ/$\mu m^2$)',fontdict={'family': 'Times New Roman','size':8})
plt.savefig('IMG/'+'fig4fa.eps', format='eps', dpi=1000)
plt.show()


#fig, ax = plt.subplots(1, 1, figsize=(1.4, 1.0))
#fig.subplots_adjust(left=0.2, right=0.9, bottom=0.15, top=0.95)
fig, ax = plt.subplots(1, 1, figsize=(2.0, 1.4))
fig.subplots_adjust(left=0.3, right=0.9, bottom=0.3, top=0.9)
plt.sca(ax)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both',direction='out',length=2,labelcolor=(0.5,0.5,0.5))
#ax.set_prop_cycle(color=['red', 'olive', 'blue', 'green', 'black', 'brown', 'orange', 'skyblue', 'springgreen', 'gray', 'darkblue', 'magenta'])
#ax.set_xlim(-0.05, 1.05)
ax.set_ylim(-50,10)
plt.xticks(np.arange(len(w_exp_mean)),fontsize=6,fontname='Times New Roman')   #[0,1,2,3,4,5,6],
plt.yticks(fontsize=6,fontname='Times New Roman')   #[0,50,100,150,200],
#ax.errorbar(np.arange(len(w_exp_mean)),w_Chen_mean,yerr=w_Chen_std,fmt='--ok',lw=1.,ms=3,mfc='w',     \
#            ecolor='gray',elinewidth=0.3,capsize=0.5,capthick=0.3,label='Chen')
ax.plot(np.arange(len(w_exp_mean)),Er_Chen_mean,'--k',lw=1.,label='$\itP_{bas}$')
ax.plot(0,Er_Chen_mean[0],'mo',ms=2)   #,mfc='w'
ax.plot(range(1,len(w_exp_mean)),Er_Chen_mean[1:],'co',ms=2)   #,mfc='w'
#ax.fill_between(np.arange(len(w_exp_mean)),Er_Chen_mean+Er_Chen_std,Er_Chen_mean-Er_Chen_std,color=(0.95,0.95,0.95),alpha=0.1)   #color=(1.0,0.9,1.0)
ax.plot(np.arange(len(w_exp_mean)),Ef_Chen_mean,'-k',lw=1.,label='$\itP_{sup}$')
ax.plot(0,Ef_Chen_mean[0],'mo',ms=2)   #,mfc='w'
ax.plot(range(1,len(w_exp_mean)),Ef_Chen_mean[1:],'co',ms=2)   #,mfc='w'
#ax.fill_between(np.arange(len(w_exp_mean)),Ef_Chen_mean+Ef_Chen_std,Ef_Chen_mean-Ef_Chen_std,color=(0.95,0.95,0.95),alpha=0.1)   #color=(1.0,0.9,1.0)
ax.plot(np.arange(len(w_exp_mean)),Er_Chen_mean_Free,'--b',lw=1.)    #,label='$\itP_{bas}$'
ax.plot(0,Er_Chen_mean_Free[0],'mo',ms=2)   #,mfc='w'
ax.plot(range(1,len(w_exp_mean)),Er_Chen_mean_Free[1:],'co',ms=2)   #,mfc='w'
#ax.fill_between(np.arange(len(w_exp_mean)),Er_Chen_mean_Free+Er_Chen_std_Free,Er_Chen_mean_Free-Er_Chen_std_Free,color=(0.95,0.95,1.0),alpha=0.1)   #color=(1.0,0.9,1.0)
ax.plot(np.arange(len(w_exp_mean)),Ef_Chen_mean_Free,'-b',lw=1.)    #,label='$\itP_{sup}$'
ax.plot(0,Ef_Chen_mean_Free[0],'mo',ms=2)   #,mfc='w'
ax.plot(range(1,len(w_exp_mean)),Ef_Chen_mean_Free[1:],'co',ms=2)   #,mfc='w'
#ax.fill_between(np.arange(len(w_exp_mean)),Ef_Chen_mean_Free+Ef_Chen_std_Free,Ef_Chen_mean_Free-Ef_Chen_std_Free,color=(0.95,0.95,1.0),alpha=0.1)   #color=(1.0,0.9,1.0)
#ax.legend(bbox_to_anchor=(1.0, 1.05),loc='upper right',fontsize=6,borderaxespad=0.)
#ax.legend(loc='lower right',fontsize=4,frameon=False)
#ax.legend(loc='lower right',prop={'family':'Times New Roman','size':4},frameon=True,framealpha=0.1)
#plt.xlabel('Distance from stimulation site',fontdict={'family': 'Times New Roman','size':8})
#plt.ylabel('Potential energy\n(fJ/$\mu m^2$)',fontdict={'family': 'Times New Roman','size':8})
plt.savefig('IMG/'+'fig4fb.eps', format='eps', dpi=1000)
plt.show()

#fig, ax = plt.subplots(1, 1, figsize=(1.4, 1.0))
#fig.subplots_adjust(left=0.2, right=0.9, bottom=0.15, top=0.95)
fig, ax = plt.subplots(1, 1, figsize=(2.0, 1.4))
fig.subplots_adjust(left=0.3, right=0.9, bottom=0.3, top=0.9)
plt.sca(ax)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both',direction='out',length=2,labelcolor=(0.5,0.5,0.5))
#ax.set_prop_cycle(color=['red', 'olive', 'blue', 'green', 'black', 'brown', 'orange', 'skyblue', 'springgreen', 'gray', 'darkblue', 'magenta'])
#ax.set_xlim(-0.05, 1.05)
ax.set_ylim(-0.25, 1.2)
plt.xticks(np.arange(len(w_exp_mean)),fontsize=6,fontname='Times New Roman')   #[0,1,2,3,4,5,6],
plt.yticks(fontsize=6,fontname='Times New Roman')   #[-0.25,0,0.25,0.5],
ax.errorbar(np.arange(len(w_exp_mean)),w_exp_mean,yerr=w_exp_std,fmt=':ok',lw=1.,ms=2,   \
            ecolor='gray',elinewidth=0.3,capsize=0.5,capthick=0.3,label='Exp')
#ax.errorbar(np.arange(len(w_exp_mean)),w_Chen_mean,yerr=w_Chen_std,fmt='--ok',lw=1.,ms=3,mfc='w',     \
#            ecolor='gray',elinewidth=0.3,capsize=0.5,capthick=0.3,label='Chen')
ax.plot(np.arange(len(w_exp_mean)),w_Chen_mean,'-k',lw=1.,label='Our model')
ax.plot(0,w_Chen_mean[0],'mo',ms=2)   #,mfc='w'
ax.plot(range(1,len(w_exp_mean)),w_Chen_mean[1:],'co',ms=2)   #,mfc='w'
ax.fill_between(np.arange(len(w_exp_mean)),w_Chen_mean+w_Chen_std,w_Chen_mean-w_Chen_std,color=(0.95,0.95,0.95),alpha=0.1)   #color=(1.0,0.9,1.0)
ax.plot(np.arange(len(w_exp_mean)),w_Chen_mean_Free,'-b',lw=1.)    #,label='Chen'
ax.plot(0,w_Chen_mean_Free[0],'mo',ms=2)   #,mfc='w'
ax.plot(range(1,len(w_exp_mean)),w_Chen_mean_Free[1:],'co',ms=2)   #,mfc='w'
#ax.fill_between(np.arange(len(w_exp_mean)),w_Chen_mean_Free+w_Chen_std_Free,w_Chen_mean_Free-w_Chen_std_Free,color=(0.95,0.95,1.0),alpha=0.1)   #color=(1.0,0.9,1.0)
#ax.legend(bbox_to_anchor=(1.0, 1.05),loc='upper right',fontsize=6,borderaxespad=0.)
ax.legend(loc='upper right',prop={'family':'Times New Roman','size':5},frameon=True,framealpha=0.1)
plt.axhline(0,linestyle='--',color='k',linewidth=0.2)
plt.xlabel('Distance from stimulation site',fontsize=8,fontname='Times New Roman')
#plt.ylabel('Weight change',fontsize=8,fontname='Times New Roman')
plt.savefig('IMG/'+'fig4fc.eps', format='eps', dpi=1000)
plt.show()


sys.exit(1)

