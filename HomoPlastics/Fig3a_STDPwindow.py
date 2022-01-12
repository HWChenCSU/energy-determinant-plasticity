
from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import json
import os, sys
mod_path = os.path.abspath(os.path.join('..','Model'))
sys.path.append(mod_path)
from oo_Parameters import *
from MorphologyData import *

#start_scope()

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
dt_array =  np.array([1.,2.5,5.,7.5,10.,12.5,15.,17.5,20.])  # interspike intervals
x_values = np.array(list(reversed(-1.*dt_array)))
x_values = np.append(x_values,dt_array)
nrDt = dt_array.size

synmodel = 'Chen'   # synmodel = 'Chen' , synmodel = 'Clopath', synmodel = 'nonPlast'
print('synmodel: ',synmodel)

init_weight = 0.5
ME_A = 0.02
ME_Vrhigh = -60*mV
ME_Ar = 0.2
MEmaxRatio = 175.0
MEtau = 2.0*second

ChenW = np.zeros((nrIn,2*nrDt))
ChenEr = np.zeros((nrIn,2*nrDt))
ChenEf = np.zeros((nrIn,2*nrDt))
ChenMEdamp = np.zeros((nrIn,2*nrDt))
ChenMEmax = np.zeros((nrIn,2*nrDt))
ChenPE = np.zeros((nrIn,2*nrDt))
for zzz in range(nrIn):
    titlestr = 'DataSTDPwindow/'+synmodel+'_'+loc1+'_'+str(init_weight)+'_'+str(ME_A)+'_'+str(ME_Vrhigh/mV)+'_'+str(ME_Ar)+'_'+str(MEmaxRatio)+'_'+str(MEtau/second)+'_'+str(d_compartm[zzz])
    data1 = open(titlestr+'_w1.txt','r')
    data2 = open(titlestr+'_w2.txt','r')
    pw1 = json.load(data1)
    pw2 = json.load(data2)
    pw = np.append(np.array(list(reversed(pw2))),pw1)
    ChenW[zzz,:] = pw
    data1.close()
    data2.close()
    data1 = open(titlestr+'_Er1.txt','r')
    data2 = open(titlestr+'_Er2.txt','r')
    pw1 = json.load(data1)
    pw2 = json.load(data2)
    pw = np.append(np.array(list(reversed(pw2))),pw1)
    ChenEr[zzz,:] = pw
    data1.close()
    data2.close()
    data1 = open(titlestr+'_Ef1.txt','r')
    data2 = open(titlestr+'_Ef2.txt','r')
    pw1 = json.load(data1)
    pw2 = json.load(data2)
    pw = np.append(np.array(list(reversed(pw2))),pw1)
    ChenEf[zzz,:] = pw
    data1.close()
    data2.close()
    data1 = open(titlestr+'_MEdamp1.txt','r')
    data2 = open(titlestr+'_MEdamp2.txt','r')
    pw1 = json.load(data1)
    pw2 = json.load(data2)
    pw = np.append(np.array(list(reversed(pw2))),pw1)
    ChenMEdamp[zzz,:] = pw
    data1.close()
    data2.close()
    data1 = open(titlestr+'_MEmax1.txt','r')
    data2 = open(titlestr+'_MEmax2.txt','r')
    pw1 = json.load(data1)
    pw2 = json.load(data2)
    pw = np.append(np.array(list(reversed(pw2))),pw1)
    ChenMEmax[zzz,:] = pw
    data1.close()
    data2.close()
    data1 = open(titlestr+'_PE1.txt','r')
    data2 = open(titlestr+'_PE2.txt','r')
    pw1 = json.load(data1)
    pw2 = json.load(data2)
    pw = np.append(np.array(list(reversed(pw2))),pw1)
    ChenPE[zzz,:] = pw
    data1.close()
    data2.close()
	
ChenWmean = 100.0*np.mean(ChenW,axis=0)/init_weight
ChenWstd = 100.0*np.std(ChenW,axis=0)    #/np.sqrt(ChenW.shape[0])
ChenErmean = np.mean(ChenEr,axis=0)
ChenErstd = np.std(ChenEr,axis=0)    #/np.sqrt(ChenEr.shape[0])
ChenEfmean = np.mean(ChenEf,axis=0)
ChenEfstd = np.std(ChenEf,axis=0)    #/np.sqrt(ChenEf.shape[0])
ChenMEdampmean = np.mean(ChenMEdamp,axis=0)
ChenMEdampstd = np.std(ChenMEdamp,axis=0)    #/np.sqrt(ChenMEdamp.shape[0])
ChenMEmaxmean = np.mean(ChenMEmax,axis=0)
ChenMEmaxstd = np.std(ChenMEmax,axis=0)    #/np.sqrt(ChenMEmax.shape[0])
ChenPEmean = np.mean(ChenPE,axis=0)
ChenPEstd = np.std(ChenPE,axis=0)/np.sqrt(ChenPE.shape[0])
ChenSmax = sign(ChenPE)*(ChenMEdamp*ChenMEmax+MESmax0)
ChenSmaxmean = np.mean(ChenSmax,axis=0)

ChenW_Free = np.zeros((nrIn,2*nrDt))
ChenEr_Free = np.zeros((nrIn,2*nrDt))
ChenEf_Free = np.zeros((nrIn,2*nrDt))
ChenPE_Free = np.zeros((nrIn,2*nrDt))
for zzz in range(nrIn):
    titlestr = 'DataSTDPwindow_Free/'+synmodel+'_'+loc1+'_'+str(init_weight)+'_'+str(ME_A)+'_'+str(ME_Vrhigh/mV)+'_'+str(ME_Ar)+'_'+str(d_compartm[zzz])
    data1 = open(titlestr+'_w1.txt','r')
    data2 = open(titlestr+'_w2.txt','r')
    pw1 = json.load(data1)
    pw2 = json.load(data2)
    pw = np.append(np.array(list(reversed(pw2))),pw1)
    ChenW_Free[zzz,:] = pw
    data1.close()
    data2.close()
    data1 = open(titlestr+'_Er1.txt','r')
    data2 = open(titlestr+'_Er2.txt','r')
    pw1 = json.load(data1)
    pw2 = json.load(data2)
    pw = np.append(np.array(list(reversed(pw2))),pw1)
    ChenEr_Free[zzz,:] = pw
    data1.close()
    data2.close()
    data1 = open(titlestr+'_Ef1.txt','r')
    data2 = open(titlestr+'_Ef2.txt','r')
    pw1 = json.load(data1)
    pw2 = json.load(data2)
    pw = np.append(np.array(list(reversed(pw2))),pw1)
    ChenEf_Free[zzz,:] = pw
    data1.close()
    data2.close()
    data1 = open(titlestr+'_PE1.txt','r')
    data2 = open(titlestr+'_PE2.txt','r')
    pw1 = json.load(data1)
    pw2 = json.load(data2)
    pw = np.append(np.array(list(reversed(pw2))),pw1)
    ChenPE_Free[zzz,:] = pw
    data1.close()
    data2.close()
	
ChenWmean_Free = 100.0*np.mean(ChenW_Free,axis=0)/init_weight
ChenWstd_Free = 100.0*np.std(ChenW_Free,axis=0)    #/np.sqrt(ChenW_Free.shape[0])
ChenErmean_Free = np.mean(ChenEr_Free,axis=0)
ChenErstd_Free = np.std(ChenEr_Free,axis=0)    #/np.sqrt(ChenEr_Free.shape[0])
ChenEfmean_Free = np.mean(ChenEf_Free,axis=0)
ChenEfstd_Free = np.std(ChenEf_Free,axis=0)    #/np.sqrt(ChenEf_Free.shape[0])
ChenPEmean_Free = np.mean(ChenPE_Free,axis=0)
ChenPEstd_Free = np.std(ChenPE_Free,axis=0)    #/np.sqrt(ChenPE_Free.shape[0])

expdt = np.array([-10,10])
expWmean = 100.0*(np.array([-0.34,0.29])+1.0)
expWstd = 100.0*np.array([0.1,0.14])
##------- STDP: [257_7]Tab.1 from [257_8]Fig.1D and Fig.7B -------
##firing rate(Hz): np.array([0.1,10,20,40,50])
## dt(ms): 10
## wmean: np.array([-0.04,0.14,0.29,0.53,0.56])
## wSEM: np.array([0.05,0.1,0.14,0.11,0.26])
## dt(ms): -10
## wmean: np.array([-0.29,-0.41,-0.34,0.56,0.75])
## wSEM: np.array([0.08,0.11,0.1,0.32,0.19])
##-------
##firing rate(Hz): 20
## dt(ms): np.array([-10,10])
## wmean: np.array([-0.34,0.29])
## wSEM: np.array([0.1,0.14])

#fig, ax = plt.subplots(1, 1, figsize=(1.4, 1.0))
fig, ax = plt.subplots(1, 1, figsize=(2.0, 1.4))
fig.subplots_adjust(left=0.3, right=0.9, bottom=0.3, top=0.9)
#fig, ax = plt.subplots(1, 1, figsize=(2, 1.4))
#fig.subplots_adjust(left=.15, right=.95, bottom=.15, top=.95)
plt.sca(ax)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(direction='out',length=2,labelcolor=(0.5,0.5,0.5))
ax.set_xlim(-22, 22)
#ax.set_ylim(-0.52, 1.02)
plt.xticks([-20,-10,0,10,20],fontsize=6,fontname='Times New Roman')
plt.yticks(fontsize=6,fontname='Times New Roman')
ax.plot(x_values,ChenSmaxmean,'--r',linewidth=1.,label='$\itP_{max}$')
ax.plot(x_values,ChenPEmean,'-k',linewidth=1.,label='$\itP$')
#ax.plot(x_values,ChenPEmean+ChenPEstd,'-k',linewidth=0.1)
#ax.plot(x_values,ChenPEmean-ChenPEstd,'-k',linewidth=0.1)
#ax.fill_between(x_values,ChenPEmean+ChenPEstd,ChenPEmean-ChenPEstd,color=(0.95,0.95,0.95),alpha=0.1)
ax.plot(x_values,ChenPEmean_Free,'-b',linewidth=1.)    #
ax.plot(x_values,ChenPEmean_Free+ChenPEstd_Free,'-b',linewidth=0.1)
ax.plot(x_values,ChenPEmean_Free-ChenPEstd_Free,'-b',linewidth=0.1)
ax.fill_between(x_values,ChenPEmean_Free+ChenPEstd_Free,ChenPEmean_Free-ChenPEstd_Free,color=(0.95,0.95,1.0),alpha=0.1)
#ax.legend(loc='lower left',prop={'family':'Times New Roman','size':4},frameon=True,framealpha=0.1)
#legend1 = ax.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left',prop={'family':'Times New Roman','size':4},borderaxespad=0.,title='$\itMER_{max}$ [fJ/($\mu m^2$*s)]')
#legend1.get_title().set_fontsize(fontsize = 6)
#legend1.get_title().set_fontname(fontname='Times New Roman')
#plt.axhline(0,linestyle='--',color='k',linewidth=0.2)
#plt.xlabel('Time intervals $\Delta t$\n(ms)',fontdict={'family': 'Times New Roman','size':8})
plt.ylabel('Potential energy\n(fJ/$\mu m^2$)',fontdict={'family': 'Times New Roman','size':8})
plt.savefig('./IMG/'+'Fig3aa.eps', format='eps', dpi=1000)
plt.show()	

fig, ax = plt.subplots(1, 1, figsize=(2.0, 1.4))
fig.subplots_adjust(left=0.3, right=0.9, bottom=0.3, top=0.9)
#fig, ax = plt.subplots(1, 1, figsize=(2, 1.4))
#fig.subplots_adjust(left=.15, right=.95, bottom=.15, top=.95)
plt.sca(ax)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(direction='out',length=2,labelcolor=(0.5,0.5,0.5))
ax.set_xlim(-22, 22)
#ax.set_ylim(-0.52, 0.52)
plt.xticks([-20,-10,0,10,20],fontsize=6,fontname='Times New Roman')#
plt.yticks(fontsize=6,fontname='Times New Roman')
ax.plot(x_values,ChenErmean,'--k',linewidth=1.,label='$\itP_{bas}$')
#ax.plot(x_values,ChenErmean+ChenErstd,'--k',linewidth=0.1)
#ax.plot(x_values,ChenErmean-ChenErstd,'--k',linewidth=0.1)
#ax.fill_between(x_values,ChenErmean+ChenErstd,ChenErmean-ChenErstd,color=(0.9,0.9,0.9),alpha=0.1)
ax.plot(x_values,ChenEfmean,'-k',linewidth=1.,label='$\itP_{sup}$')
#ax.plot(x_values,ChenEfmean+ChenEfstd,'-k',linewidth=0.1)
#ax.plot(x_values,ChenEfmean-ChenEfstd,'-k',linewidth=0.1)
#ax.fill_between(x_values,ChenEfmean+ChenEfstd,ChenEfmean-ChenEfstd,color=(0.9,0.9,0.9),alpha=0.1)
ax.plot(x_values,ChenErmean_Free,'--b',linewidth=1.)    #,label='$\itP_{bas}$'
#ax.plot(x_values,ChenErmean+ChenErstd,'--k',linewidth=0.1)
#ax.plot(x_values,ChenErmean-ChenErstd,'--k',linewidth=0.1)
#ax.fill_between(x_values,ChenErmean+ChenErstd,ChenErmean-ChenErstd,color=(0.9,0.9,0.9),alpha=0.1)
ax.plot(x_values,ChenEfmean_Free,'-b',linewidth=1.)    #,label='$\itP_{sup}$'
#ax.plot(x_values,ChenEfmean+ChenEfstd,'-k',linewidth=0.1)
#ax.plot(x_values,ChenEfmean-ChenEfstd,'-k',linewidth=0.1)
#ax.fill_between(x_values,ChenEfmean+ChenEfstd,ChenEfmean-ChenEfstd,color=(0.9,0.9,0.9),alpha=0.1)
ax.legend(loc='lower left',prop={'family':'Times New Roman','size':5},frameon=True,framealpha=0.1)
#plt.axhline(0,linestyle='--',color='k',linewidth=0.2)
#plt.xlabel('Time intervals $\Delta t$\n(ms)',fontsize=8,fontname='Times New Roman')
plt.ylabel('Potential energy\n(fJ/$\mu m^2$)',fontsize=8,fontname='Times New Roman')
plt.savefig('./IMG/'+'Fig3ab.eps', format='eps', dpi=1000)
plt.show()	

fig, ax = plt.subplots(1, 1, figsize=(2.0, 1.4))
fig.subplots_adjust(left=0.3, right=0.9, bottom=0.3, top=0.9)
#fig, ax = plt.subplots(1, 1, figsize=(2, 1.4))
#fig.subplots_adjust(left=.15, right=.95, bottom=.15, top=.95)
plt.sca(ax)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(direction='out',length=2,labelcolor=(0.5,0.5,0.5))
ax.set_xlim(-22, 22)
#ax.set_ylim(-2, 202)
plt.xticks([-20,-10,0,10,20],fontsize=6,fontname='Times New Roman')#
plt.yticks(fontsize=6,fontname='Times New Roman')    #[0,25,50,75,100,125,150,175,200],
ax.plot(x_values,ChenWmean,'-k',linewidth=1.,label='Our model')
ax.plot(x_values,ChenWmean+ChenWstd,'-k',linewidth=0.1)
ax.plot(x_values,ChenWmean-ChenWstd,'-k',linewidth=0.1)
ax.fill_between(x_values,ChenWmean+ChenWstd,ChenWmean-ChenWstd,color=(0.95,0.95,0.95),alpha=0.1)
ax.plot(x_values,ChenWmean_Free,'-b',linewidth=1.)    #,label='Chen'
#ax.plot(x_values,ChenWmean+ChenWstd,'-k',linewidth=0.1)
#ax.plot(x_values,ChenWmean-ChenWstd,'-k',linewidth=0.1)
#ax.fill_between(x_values,ChenWmean+ChenWstd,ChenWmean-ChenWstd,color=(0.9,0.9,0.9),alpha=0.1)
ax.errorbar(expdt,expWmean,yerr=expWstd,color='black',ecolor='gray',elinewidth=0.3,capsize=0.5,capthick=0.3,fmt='o',ms=0.1,label='Exp')
ax.legend(loc='upper left',prop={'family':'Times New Roman','size':5},frameon=True,framealpha=0.1)
plt.axhline(100,linestyle='--',color='k',linewidth=0.2)
plt.xlabel('Time intervals $\Delta t$ (ms)',fontsize=8,fontname='Times New Roman')
plt.ylabel('Normalised weight\n(%)',fontsize=8,fontname='Times New Roman')
plt.savefig('./IMG/'+'Fig3ac.eps', format='eps', dpi=1000)
plt.show()	


