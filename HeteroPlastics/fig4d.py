######################################################
# 第9分支、远端、不同频率刺激下，采用Chen可塑性模型，从远端到近端的势能和权值分布的图形显示
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
abranch = 9
print('分支编号：',abranch)
signal_rate = 100.0 # activation rate of the pools
print('刺激频率：',signal_rate)
init_weight = 0.5  # initial weight
synmodel = 'Chen'   # synmodel = 'Chen' , synmodel = 'Clopath', synmodel = 'nonPlast'
print('突触可塑性模型：',synmodel)
homoloc = 'distal'   #'proximal','distal'
print('分支位置：',homoloc)
nr_clst = 2   # nr of synapses per compartment

#------------------------------	
titlestr = 'Data/'+synmodel+'_'+loc1+'_'+str(nr_clst)+'_'+str(init_weight)+'_'+str(abranch)+'_'+homoloc+'_'+str(signal_rate)
data0 = open(titlestr+'_PEbranch.txt','r')
PEbranch = np.array(json.load(data0))	
data0.close()
data0 = open(titlestr+'_MEmaxbranch.txt','r')
MEmaxbranch = np.array(json.load(data0))	
data0.close()
data0 = open(titlestr+'_MEdampbranch.txt','r')
MEdampbranch = np.array(json.load(data0))	
data0.close()
data0 = open(titlestr+'_Erbranch.txt','r')
Erbranch = np.array(json.load(data0))	
data0.close()
data0 = open(titlestr+'_Efbranch.txt','r')
Efbranch = np.array(json.load(data0))	
data0.close()
data0 = open(titlestr+'_wbranch.txt','r')
wbranch = np.array(json.load(data0))
data0.close()

titlestr = 'Data/'+synmodel+'_'+loc1+'_'+str(nr_clst)+'_'+str(init_weight)+'_'+str(abranch)+'_'+homoloc+'_'+str(signal_rate)
data0 = open(titlestr+'_PEbranch_Free.txt','r')
PEbranch_Free = np.array(json.load(data0))	
data0.close()
data0 = open(titlestr+'_Erbranch_Free.txt','r')
Erbranch_Free = np.array(json.load(data0))	
data0.close()
data0 = open(titlestr+'_Efbranch_Free.txt','r')
Efbranch_Free = np.array(json.load(data0))	
data0.close()
data0 = open(titlestr+'_wbranch_Free.txt','r')
wbranch_Free = np.array(json.load(data0))
data0.close()
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
plt.xticks([0,1,2,3,4,5,6],fontsize=6,fontname='Times New Roman')   #[0,1,2,3,4,5,6],
plt.yticks(fontsize=6,fontname='Times New Roman')   #[-100,-50,0,50,100],
PEmax = sign(PEbranch)*(MEdampbranch*MEmaxbranch+MESmax0)
ax.plot(range(len(PEmax)),PEmax,'--r',linewidth=1.,label='$\itP_{max}$')
ax.plot(0,PEmax[0],'mo',ms=2)
ax.plot(range(1,len(PEmax)),PEmax[1:],'co',ms=2)
ax.plot(range(len(PEbranch)),PEbranch,'-k',linewidth=1.,label='$\itP$')
ax.plot(0,PEbranch[0],'mo',ms=2)
ax.plot(range(1,len(PEbranch)),PEbranch[1:],'co',ms=2)
ax.plot(range(len(PEbranch_Free)),PEbranch_Free,'-b',linewidth=1.)    #,label='$\itP$'
ax.plot(0,PEbranch_Free[0],'mo',ms=2)
ax.plot(range(1,len(PEbranch_Free)),PEbranch_Free[1:],'co',ms=2)
#ax.legend(loc='lower right',fontsize=4,frameon=False)
ax.legend(loc='lower right',prop={'family':'Times New Roman','size':5},frameon=True,framealpha=0.1)
#plt.xlabel('Distance from stimulation site',fontdict={'family': 'Times New Roman','size':8})
#plt.ylabel('Potential energy\n(fJ/$\mu m^2$)',fontdict={'family': 'Times New Roman','size':8})
plt.savefig('IMG/'+'fig4da.eps', format='eps', dpi=1000)
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
ax.set_ylim(-55, 10)
plt.xticks([0,1,2,3,4,5,6],fontsize=6,fontname='Times New Roman')   #[0,1,2,3,4,5,6],
plt.yticks(fontsize=6,fontname='Times New Roman')   #[-100,-50,0,50,100],
ax.plot(range(len(Erbranch)),Erbranch,'--k',linewidth=1.,label='$\itP_{bas}$')
ax.plot(0,Erbranch[0],'mo',ms=2)
ax.plot(range(1,len(Erbranch)),Erbranch[1:],'co',ms=2)
ax.plot(range(len(Efbranch)),Efbranch,'-k',linewidth=1.,label='$\itP_{sup}$')
ax.plot(0,Efbranch[0],'mo',ms=2)
ax.plot(range(1,len(Efbranch)),Efbranch[1:],'co',ms=2)
ax.plot(range(len(Erbranch_Free)),Erbranch_Free,'--b',linewidth=1.)    #,label='$\itP_{bas}$'
ax.plot(0,Erbranch_Free[0],'mo',ms=2)
ax.plot(range(1,len(Erbranch_Free)),Erbranch_Free[1:],'co',ms=2)
ax.plot(range(len(Efbranch_Free)),Efbranch_Free,'-b',linewidth=1.)    #,label='$\itP_{sup}$'
ax.plot(0,Efbranch_Free[0],'mo',ms=2)
ax.plot(range(1,len(Efbranch_Free)),Efbranch_Free[1:],'co',ms=2)
#ax.legend(loc='lower right',fontsize=4,frameon=False)
ax.legend(loc='lower right',prop={'family':'Times New Roman','size':5},frameon=True,framealpha=0.1)
#plt.xlabel('Distance from stimulation site',fontdict={'family': 'Times New Roman','size':8})
#plt.ylabel('Potential energy\n(fJ/$\mu m^2$)',fontdict={'family': 'Times New Roman','size':8})
plt.savefig('IMG/'+'fig4db.eps', format='eps', dpi=1000)
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
ax.set_ylim(-0.2, 1.2)
plt.xticks([0,1,2,3,4,5,6],fontsize=6,fontname='Times New Roman')   #[0,1,2,3,4,5,6],
plt.yticks(fontsize=6,fontname='Times New Roman')   #[-0.25,0,0.25,0.5],
wbranch = wbranch-init_weight
ax.plot(range(len(wbranch)),wbranch,'-k',linewidth=1.,label='Chen')
ax.plot(0,wbranch[0],'mo',ms=2)
ax.plot(range(1,len(wbranch)),wbranch[1:],'co',ms=2)
wbranch_Free = wbranch_Free-init_weight
ax.plot(range(len(wbranch_Free)),wbranch_Free,'-b',linewidth=1.)    #,label='Chen'
ax.plot(0,wbranch_Free[0],'mo',ms=2)
ax.plot(range(1,len(wbranch_Free)),wbranch_Free[1:],'co',ms=2)
plt.axhline(0,linestyle='--',color='k',linewidth=0.2)
plt.xlabel('Distance from stimulation site',fontdict={'family': 'Times New Roman','size':8})
#plt.ylabel('Weight change',fontsize=8,fontname='Times New Roman')
plt.savefig('IMG/'+'fig4dc.eps', format='eps', dpi=1000)
plt.show()



sys.exit(1)


