
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
hz_array = np.array([1.,3.,5.,10.,20.,30.,40.,50.]) 
nrHz = hz_array.size

synmodel = 'Chen'   # synmodel = 'Chen' , synmodel = 'Clopath', synmodel = 'nonPlast'
print('synmodel: ',synmodel)
ME_Ascale = 4.0
nr_clst = 1

init_weight = 0.5
ME_A = 0.02
ME_Vrhigh = -60*mV
ME_Ar = 0.2
MEmaxRatio = 175.0
MEtau = 2.0*second

ChenW = np.zeros((nrIn,nrHz))
ChenEr = np.zeros((nrIn,nrHz))
ChenEf = np.zeros((nrIn,nrHz))
ChenMEdamp = np.zeros((nrIn,nrHz))
ChenMEmax = np.zeros((nrIn,nrHz))
ChenPE = np.zeros((nrIn,nrHz))
for zzz in range(nrIn):
    titlestr = 'DataPoissonInput/'+synmodel+'_'+loc1+'_'+str(ME_Ascale)+'_'+str(nr_clst)+'_'+str(init_weight)+'_'+str(ME_A)+'_'+str(ME_Vrhigh/mV)+'_'+str(ME_Ar)+'_'+str(MEmaxRatio)+'_'+str(MEtau/second)+'_'+str(d_compartm[zzz])
    data1 = open(titlestr+'_w1.txt','r')
    ChenW[zzz,:] = json.load(data1)
    data1.close()
    data1 = open(titlestr+'_Er1.txt','r')
    ChenEr[zzz,:] = json.load(data1)
    data1.close()
    data1 = open(titlestr+'_Ef1.txt','r')
    ChenEf[zzz,:] = json.load(data1)
    data1.close()
    data1 = open(titlestr+'_MEdamp1.txt','r')
    ChenMEdamp[zzz,:] = json.load(data1)
    data1.close()
    data1 = open(titlestr+'_MEmax1.txt','r')
    ChenMEmax[zzz,:] = json.load(data1)
    data1.close()
    data1 = open(titlestr+'_PE1.txt','r')
    ChenPE[zzz,:] = json.load(data1)
    data1.close()
	
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
ChenPEstd = np.std(ChenPE,axis=0)    #/np.sqrt(ChenPE.shape[0])

ChenW_Free = np.zeros((nrIn,nrHz))
ChenEr_Free = np.zeros((nrIn,nrHz))
ChenEf_Free = np.zeros((nrIn,nrHz))
ChenPE_Free = np.zeros((nrIn,nrHz))
for zzz in range(nrIn):
    titlestr = 'DataPoissonInput/'+synmodel+'_'+loc1+'_'+str(ME_Ascale)+'_'+str(nr_clst)+'_'+str(init_weight)+'_'+str(ME_A)+'_'+str(ME_Vrhigh/mV)+'_'+str(ME_Ar)+'_'+str(MEmaxRatio)+'_'+str(MEtau/second)+'_'+str(d_compartm[zzz])
    data1 = open(titlestr+'_w1_Free.txt','r')
    ChenW_Free[zzz,:] = json.load(data1)
    data1.close()
    data1 = open(titlestr+'_Er1_Free.txt','r')
    ChenEr_Free[zzz,:] = json.load(data1)
    data1.close()
    data1 = open(titlestr+'_Ef1_Free.txt','r')
    ChenEf_Free[zzz,:] = json.load(data1)
    data1.close()
    data1 = open(titlestr+'_PE1_Free.txt','r')
    ChenPE_Free[zzz,:] = json.load(data1)
    data1.close()
	
ChenWmean_Free = 100.0*np.mean(ChenW_Free,axis=0)/init_weight
ChenWstd_Free = 100.0*np.std(ChenW_Free,axis=0)    #/np.sqrt(ChenW_Free.shape[0])
ChenErmean_Free = np.mean(ChenEr_Free,axis=0)
ChenErstd_Free = np.std(ChenEr_Free,axis=0)    #/np.sqrt(ChenEr_Free.shape[0])
ChenEfmean_Free = np.mean(ChenEf_Free,axis=0)
ChenEfstd_Free = np.std(ChenEf_Free,axis=0)    #/np.sqrt(ChenEf_Free.shape[0])
ChenPEmean_Free = np.mean(ChenPE_Free,axis=0)
ChenPEstd_Free = np.std(ChenPE_Free,axis=0)    #/np.sqrt(ChenPE_Free.shape[0])


# [301_50]p.2,Fig.2, and Fig.4A
expfrq = np.array([1,3,5,10,50])
expW1mean = np.array([88.,81.,84.,101.,121.])
expW1std = np.array([1.,2.,6.,4.,7.])

#fig, ax = plt.subplots(1, 1, figsize=(1.4, 1.0))
fig, ax = plt.subplots(1, 1, figsize=(2.0, 1.4))
fig.subplots_adjust(left=0.3, right=0.9, bottom=0.3, top=0.9)
#fig, ax = plt.subplots(1, 1, figsize=(2, 1.4))
#fig.subplots_adjust(left=.15, right=.95, bottom=.15, top=.95)
plt.sca(ax)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(direction='out',length=2,labelcolor=(0.5,0.5,0.5))
ax.set_xlim(-2, 52)
#ax.set_ylim(-0.52, 1.02)
plt.xticks([0,10,20,30,40,50],fontsize=6,fontname='Times New Roman')
plt.yticks(fontsize=6,fontname='Times New Roman')
ax.plot(hz_array,sign(ChenPEmean)*(ChenMEdampmean*ChenMEmaxmean+MESmax0),'--r',linewidth=1.,label='$\itP_{max}$')
ax.plot(hz_array,ChenPEmean,'-k',linewidth=1.,label='$\itP$')
#ax.plot(hz_array,ChenPEmean+ChenPEstd,'-k',linewidth=0.1)
#ax.plot(hz_array,ChenPEmean-ChenPEstd,'-k',linewidth=0.1)
#ax.fill_between(hz_array,ChenPEmean+ChenPEstd,ChenPEmean-ChenPEstd,color=(0.9,0.9,0.9),alpha=0.1)
ax.plot(hz_array,ChenPEmean_Free,'-b',linewidth=1.)    #,label='$\itP$'
ax.plot(hz_array,ChenPEmean_Free+ChenPEstd_Free,'-b',linewidth=0.1)
ax.plot(hz_array,ChenPEmean_Free-ChenPEstd_Free,'-b',linewidth=0.1)
ax.fill_between(hz_array,ChenPEmean_Free+ChenPEstd_Free,ChenPEmean_Free-ChenPEstd_Free,color=(0.95,0.95,1.0),alpha=0.1)
ax.legend(loc='lower right',prop={'family':'Times New Roman','size':5},frameon=True,framealpha=0.1)
#legend1 = ax.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left',prop={'family':'Times New Roman','size':4},borderaxespad=0.,title='$\itMER_{max}$ [fJ/($\mu m^2$*s)]')
#legend1.get_title().set_fontsize(fontsize = 6)
#legend1.get_title().set_fontname(fontname='Times New Roman')
#plt.axhline(0,linestyle='--',color='k',linewidth=0.2)
#plt.xlabel('Input frequency (Hz)',fontdict={'family': 'Times New Roman','size':8})
#plt.ylabel('Potential energy\n(fJ/$\mu m^2$)',fontdict={'family': 'Times New Roman','size':8})
plt.savefig('./IMG/'+'Fig3ba.eps', format='eps', dpi=1000)
plt.show()	

fig, ax = plt.subplots(1, 1, figsize=(2.0, 1.4))
fig.subplots_adjust(left=0.3, right=0.9, bottom=0.3, top=0.9)
#fig, ax = plt.subplots(1, 1, figsize=(2, 1.4))
#fig.subplots_adjust(left=.15, right=.95, bottom=.15, top=.95)
plt.sca(ax)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(direction='out',length=2,labelcolor=(0.5,0.5,0.5))
ax.set_xlim(-2, 52)
#ax.set_ylim(-0.52, 1.02)
plt.xticks([0,10,20,30,40,50],fontsize=6,fontname='Times New Roman')
plt.yticks(fontsize=6,fontname='Times New Roman')
ax.plot(hz_array,ChenErmean,'--k',linewidth=1.,label='$\itP_{bas}$')
#ax.plot(hz_array,ChenErmean+ChenErstd,'--k',linewidth=0.1)
#ax.plot(hz_array,ChenErmean-ChenErstd,'--k',linewidth=0.1)
#ax.fill_between(hz_array,ChenErmean+ChenErstd,ChenErmean-ChenErstd,color=(0.9,0.9,0.9),alpha=0.1)
ax.plot(hz_array,ChenEfmean,'-k',linewidth=1.,label='$\itP_{sup}$')
#ax.plot(hz_array,ChenEfmean+ChenEfstd,'-k',linewidth=0.1)
#ax.plot(hz_array,ChenEfmean-ChenEfstd,'-k',linewidth=0.1)
#ax.fill_between(hz_array,ChenEfmean+ChenEfstd,ChenEfmean-ChenEfstd,color=(0.9,0.9,0.9),alpha=0.1)
ax.plot(hz_array,ChenErmean_Free,'--b',linewidth=1.)    #,label='$\itP_{bas}$'
#ax.plot(hz_array,ChenErmean+ChenErstd,'--k',linewidth=0.1)
#ax.plot(hz_array,ChenErmean-ChenErstd,'--k',linewidth=0.1)
#ax.fill_between(hz_array,ChenErmean+ChenErstd,ChenErmean-ChenErstd,color=(0.9,0.9,0.9),alpha=0.1)
ax.plot(hz_array,ChenEfmean_Free,'-b',linewidth=1.)    #,label='$\itP_{sup}$'
#ax.plot(hz_array,ChenEfmean+ChenEfstd,'-k',linewidth=0.1)
#ax.plot(hz_array,ChenEfmean-ChenEfstd,'-k',linewidth=0.1)
#ax.fill_between(hz_array,ChenEfmean+ChenEfstd,ChenEfmean-ChenEfstd,color=(0.9,0.9,0.9),alpha=0.1)
#ax.legend(loc='lower left',prop={'family':'Times New Roman','size':4},frameon=True,framealpha=0.1)
#plt.axhline(0,linestyle='--',color='k',linewidth=0.2)
#plt.xlabel('Input frequency (Hz)',fontdict={'family': 'Times New Roman','size':8})
#plt.ylabel('Potential energy\n(fJ/$\mu m^2$)',fontdict={'family': 'Times New Roman','size':8})
plt.savefig('./IMG/'+'Fig3bb.eps', format='eps', dpi=1000)
plt.show()	

fig, ax = plt.subplots(1, 1, figsize=(2.0, 1.4))
fig.subplots_adjust(left=0.3, right=0.9, bottom=0.3, top=0.9)
#fig, ax = plt.subplots(1, 1, figsize=(2, 1.4))
#fig.subplots_adjust(left=.15, right=.95, bottom=.15, top=.95)
plt.sca(ax)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(direction='out',length=2,labelcolor=(0.5,0.5,0.5))
ax.set_xlim(-2, 52)
#ax.set_ylim(-2, 202)
plt.xticks([0,10,20,30,40,50],fontsize=6,fontname='Times New Roman')
plt.yticks(fontsize=6,fontname='Times New Roman')    #[0,25,50,75,100,125,150,175,200],
ax.plot(hz_array,ChenWmean,'-k',linewidth=1.,label='Chen')
ax.plot(hz_array,ChenWmean+ChenWstd,'-k',linewidth=0.1)
ax.plot(hz_array,ChenWmean-ChenWstd,'-k',linewidth=0.1)
ax.fill_between(hz_array,ChenWmean+ChenWstd,ChenWmean-ChenWstd,color=(0.95,0.95,0.95),alpha=0.1)
ax.plot(hz_array,ChenWmean_Free,'-b',linewidth=1.)    #,label='Chen'
#ax.plot(hz_array,ChenWmean+ChenWstd,'-k',linewidth=0.1)
#ax.plot(hz_array,ChenWmean-ChenWstd,'-k',linewidth=0.1)
#ax.fill_between(hz_array,ChenWmean+ChenWstd,ChenWmean-ChenWstd,color=(0.9,0.9,0.9),alpha=0.1)
ax.errorbar(expfrq,expW1mean,yerr=expW1std,color='black',ecolor='gray',elinewidth=0.3,capsize=0.5,capthick=0.3,fmt='o',ms=0.1,label='Exp')
#ax.legend(loc='upper left',prop={'family':'Times New Roman','size':4},frameon=False,framealpha=0.1)
plt.axhline(100,linestyle='--',color='k',linewidth=0.2)
plt.xlabel('Input frequency (Hz)',fontsize=8,fontname='Times New Roman')
#plt.ylabel('Normalised weight\n(%)',fontsize=8,fontname='Times New Roman')
plt.savefig('./IMG/'+'Fig3bc.eps', format='eps', dpi=1000)
plt.show()	
