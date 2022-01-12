# This script contains the Brianmodel class
# Calling makeneuron_ca() on a Brianmodel object will create a biophysical neuron
# Multiple other functions allow for plotting, animating, ...



from __future__ import division

#folder with parameters, equations and morphology

import os, sys
mod_path = os.path.abspath(os.path.join('..','Model'))
sys.path.append(mod_path)
from copy import deepcopy 
import itertools as itools  

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg'
import matplotlib.colors as colorz
import matplotlib.cm as clrm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import brian2 as br2
from brian2 import uF, cm, um, ohm, ms, siemens, mV, nA, us,psiemens

# This is the 3D plotting toolkit
from mpl_toolkits.mplot3d import Axes3D

#import parameters and equations for neuron
from oo_Parameters import *
from oo_equations_AMPAplast import *
from MorphologyData import *
from Visualisation_functions import *
from oo_initScripts import set_init_nrn

br2.start_scope()
br2.defaultclock.dt = defaultclock.dt



class BRIANModel(object):
    """
    Neuron object in brian2
    """
    def __init__(self, swc_model):
        """
        Parameters
        ----------
        swc_model: a char
            path of the file containing the neuron model in .swc format
        """

        # Brian morphology
        self.morpho = br2.Morphology.from_file(swc_model)
        morpho = self.morpho
        # Store compartment numbers       
        self.segment,self.segment_swc,self.compStart,self.compEnd = get_swc(swc_model) 
        # Initialise an dictionary for distances to the soma per compartment
        self.distances = {}
        # Initialise an dictionary for lines to plot the neuron
        self.lines = {}


        # Add the first section as soma
        self.sections = {morpho.type: [self.morpho[0], 0, 0]}
        
        # Set a name and distances for the soma
        self.sections['soma'][0].name = 'soma'
        self.sections['soma'][0].f_x = self.morpho[0].x/meter
        self.sections['soma'][0].f_y = self.morpho[0].y/meter
        self.sections['soma'][0].f_z = self.morpho[0].z/meter
        self.sections['soma'][0].dist = 0
        self.distances['soma'] = [0.]
        
        # Initialize the dendrites numerotation
        dend_b = 0

        # Register soma's children in a sections dictionary
        for sec in morpho.children:

            # Create an attribut "name" for all children of the soma
            if str(sec.type) == "dend":
                sec.name = sec.type[:4]+"_"+str(dend_b)
                dend_b += 1
            else:
                sec.name = sec.type

            # store final coordinates of the parent (=soma) segment
            sec.f_x = self.morpho[0].x[0]/meter
            sec.f_y = self.morpho[0].y[0]/meter
            sec.f_z = self.morpho[0].z[0]/meter
            sec.dist = self.distances['soma'][0]
            
            # add distances to the parent
            self.distances = calc_dist(self.distances, sec)

            # get the coordinates for all compartments in this section
            xn = sec.x/meter
            yn = sec.y/meter
            zn = sec.z/meter
            # get first coordinates (and make integer)
            a=(int(round(xn[0]*1e9)),int(round(yn[0]*1e9)),int(round(zn[0]*1e9))) 
            # id for the section (they correspond to lnum in .swc)
            line_num = self.segment[a]
            # add id and section to the 'sections' dictionary
            self.sections[sec.name] = [sec,line_num,line_num]            
                       
            
        # Initialize the level value
        level = [sec for sec in morpho.children]
        while level != []:
            for i, sec in enumerate(level):
                for j, child in enumerate(sec.children):
                    
                    # Create an attribut "name" for all children of sec
                    name = sec.name + str(j)
                    child.name = name
                    # Store parent final coordinates
                    child.f_x = sec.x[-1]/meter
                    child.f_y = sec.y[-1]/meter
                    child.f_z = sec.z[-1]/meter
                    # Store distances to the soma
                    child.dist = self.distances[sec.name][-1]
                    self.distances = calc_dist(self.distances, child)
                    # Get the coordinates for all compartments in this section
                    xn = child.x/meter
                    yn = child.y/meter
                    zn = child.z/meter
                    # get first coordinates (and make integer)
                    a=(int(round(xn[0]*1e9)),int(round(yn[0]*1e9)),int(round(zn[0]*1e9)))
                    
                    # id for the section (corresponds to lnum in .swc)
                    line_num = self.segment[a]
                    # add id and section to the 'sections' dictionary
                    self.sections[name] = [child, line_num,line_num]            
            level = [sec.children for sec in level]
            # Flatten the list at this level
            level = [sublist for sl in level for sublist in sl]
            
    ################################################################################
    # THE FUNCTION BELOW CAN BE CALLED TO CREATE A BIOPHYSICAL NEURON
    ################################################################################

    def makeNeuron_Ca(self,morphodata):
        """return spatial neuron"""
        # Set Biophysics
        neuron = self.biophysics(morphodata)
        return neuron
        
    def biophysics(self,morpho_data):
        """Inserting biophysics"""
                    
        neuron = br2.SpatialNeuron(morphology=self.morpho, model=eqs, \
            Cm=Capacit, Ri=R_axial, threshold  = "v/mV>0", refractory = "v/mV > -10",
            threshold_location = 0, reset = 's_trace += x_reset*(taux/ms)',method='heun') #
		
        # define the different parts of the neuron
        N_soma = neuron[morpho_data['soma'][0]:morpho_data['soma'][-1]+1]
        N_axon = neuron[morpho_data['axon'][0]:morpho_data['axon'][-1]+1]
        N_basal = neuron[morpho_data['basal'][0]:morpho_data['basal'][-1]+1]
        N_apical = neuron[morpho_data['apical'][0]:morpho_data['apical'][-1]+1]
        Theta_low = morpho_data['thetalow']*mV
        
        # insert leak conductance
        neuron.gLeak = g_leak
               
        # noise
        neuron.noise_sigma = 0*pA # initial value membrane voltage
        neuron.noise_avg = 0*pA # initial value membrane voltage
        N_soma.noise_sigma = noise_std # initial value membrane voltage
        N_soma.noise_avg = noise_mean # initial value membrane voltage    
        
        ####################
        # ACTIVE CHANNELS
        ####################
        
        # Na channels soma, axon, apical dendrites
        N_soma.gNav = somaNa
        N_axon.gNav = axonNa
        N_apical.gNav = apicalNa
        neuron.thi1 = thi1_all
        N_axon.thi1 = thi1_axn
        neuron.thi2 = thi2_all
        N_axon.thi2 = thi2_axn
        
        #Kv channels
        N_soma.gKv = somagKv
        N_basal.gKv = dendgKv
        N_apical.gKv = dendgKv
        N_axon.gKv = axongKv
        
        #Ca channels sina
        N_soma.gCav = ratio_ca*somaCa
        N_soma.gIt = (1-ratio_ca)*somaCa
        
        #Ka channels soma
        N_soma.gKa_prox = somaKap
        
        #Ka channels dendrites, Na channels basal dendrites, Ca channels dendrites, axon initial segment
        for sec in self.sections:
            secNr = self.sections[sec][2]
            seclen = len(self.sections[sec][0].x)
            
            #BASAL
            if secNr in morpho_data['basal']:      
                # decreasing Na channels
                gNa_diff = 0.5*np.array(self.distances[sec][:])*psiemens/um**2
                neuron[secNr:secNr+seclen].gNav = np.multiply(basalNa - gNa_diff,basalNa - gNa_diff>0 ) 
                
                # increasing Ka channels
                gKa_diff = 0.7*np.array(self.distances[sec][:])*psiemens/um**2
                ratio_A = np.multiply(1. - (1./300.)*np.array(self.distances[sec][:]),1. - (1./300.)*np.array(self.distances[sec][:])>0)
                neuron[secNr:secNr+seclen].gKa_prox = ratio_A*np.multiply(basalKa + gKa_diff,basalKa + gKa_diff>0 )
                neuron[secNr:secNr+seclen].gKa_dist = (1.-ratio_A)*np.multiply(basalKa + gKa_diff,basalKa + gKa_diff>0 )
                
                # Ca channels
                neuron[secNr:secNr+seclen].gCav = dendCa*ratio_ca*(np.array(self.distances[sec][:])>30) + somaCa*ratio_ca*(np.array(self.distances[sec][:])<=30)
                neuron[secNr:secNr+seclen].gIt = dendCa*(1.-ratio_ca)*(np.array(self.distances[sec][:])>30) + somaCa*(1.-ratio_ca)*(np.array(self.distances[sec][:])<=30) 
                
                
                #spines
                addSpines = np.array(self.distances[sec][:]) > spinedist
                noSpines = np.array(self.distances[sec][:]) <= spinedist                
                neuron[secNr:secNr+seclen].gLeak = noSpines*g_leak + addSpines*g_leak_dend
                neuron[secNr:secNr+seclen].Cm = noSpines*Capacit + addSpines*Capacit_dend
                
            #APICAL
            if secNr in morpho_data['apical']:       
                #ratio of Ka channels
                ratio_A = np.multiply(1. - (1./300.)*np.array(self.distances[sec][:]),1. - (1./300.)*np.array(self.distances[sec][:])>0)
                neuron[secNr:secNr+seclen].gKa_prox = ratio_A*apicalKa
                neuron[secNr:secNr+seclen].gKa_dist = (1.-ratio_A)*apicalKa
                
                # Ca channels
                neuron[secNr:secNr+seclen].gCav = dendCa*ratio_ca*(np.array(self.distances[sec][:])>30) + somaCa*ratio_ca*(np.array(self.distances[sec][:])<=30)
                neuron[secNr:secNr+seclen].gIt = dendCa*(1.-ratio_ca)*(np.array(self.distances[sec][:])>30) + somaCa*(1.-ratio_ca)*(np.array(self.distances[sec][:])<=30) 
                
                #spines
                addSpines = np.array(self.distances[sec][:]) > spinedist
                noSpines = np.array(self.distances[sec][:]) <= spinedist                
                neuron[secNr:secNr+seclen].gLeak = noSpines*g_leak + addSpines*g_leak_dend
                neuron[secNr:secNr+seclen].Cm = noSpines*Capacit + addSpines*Capacit_dend
                
            #AXON
            if secNr in morpho_data['axon']:                       
                #KL current
                addKL = np.array(self.distances[sec][:]) > 35           
                neuron[secNr:secNr+seclen].gKL = addKL*axongL
                neuron[1:6].gKv = [40.,100.,500.,500.,500.]*psiemens/um**2 
                neuron[1:6].gKL =  [20.,35.,125.,250.,0]*psiemens/um**2 
                neuron[1:6].gNav = 4*np.array([8000.,7000.,5000.,5000.,5000.])*psiemens/um**2
#                neuron[1:6].gNav = [8000.,7000.,5000.,5000.,5000.]*psiemens/um**2
                neuron[1:3].gCav = somaCa*ratio_ca
                neuron[1:3].gIt = somaCa*(1.-ratio_ca)
                                
        
        # SET INITIAL VALUES
        set_init_nrn(neuron,Theta_low)
                
        
        return neuron
   
    ################################################################################
    # The functions below are mainly for visualization of the neuron morphology
    ################################################################################
        
    def print_dist(self, sec_number):
        ''' print the distance and diameter of a section to the soma'''
        for sec in self.sections:
                for ii in range(len(self.sections[sec][0].x)):
                    if (self.sections[sec][2]+ii == sec_number):
#                        print( 'Section '+ str(sec_number)+ ', part of: '+str(sec))
#                        print( 'Distance to soma: '+ str(self.distances[sec][ii]))
#                        print( 'Diameter: '+ str(self.sections[sec][0].diameter[ii]*1.e6))
                        sectiondistance = self.distances[sec][ii]
                        sectiondiameter = self.sections[sec][0].diameter[ii]*1.e6
        return [sectiondistance,sectiondiameter]
        
    def save_dist_vs_nr(self,maxNr):
        ''' save the distance and diameter in function of section nr'''
        dist_nr = np.zeros(maxNr)
        diam_nr = np.zeros(maxNr)
        print('saving distances')
        for sec in self.sections:
                for ii in range(len(self.sections[sec][0].x)):
                    dist_nr[self.sections[sec][2]+ii] = self.distances[sec][ii]
                    diam_nr[self.sections[sec][2]+ii] = self.sections[sec][0].diameter[ii]*1.e6
        return dist_nr,diam_nr
    
    def calc_distCompartments(self,sec_range,distances):
        ''' calculate the terminal ends of dendrites '''
        term_vec = [0]
        dist_vec = distances[sec_range]
        for jj in range(len(dist_vec)-1):
            if np.abs(dist_vec[jj+1]-dist_vec[jj])>20:
                term_vec = np.append(term_vec,np.array([sec_range[0]+jj]),axis=0)
        return term_vec,dist_vec

    def calc_morphdata(self,segm_range=[0]):
        ''' calculate proximal and distal compartments of thin terminal dendrites(sections), by huanwen chen '''
        dist_comps = []
        prox_comps = []
        for sec in self.sections:
            if str(sec) != 'soma' and str(sec) != 'axon':
                if bool(self.sections[sec][0].children):
                    continue
                if self.sections[sec][2] in segm_range:
                    prox_comps += [self.sections[sec][2]]
                    dist_comps += [self.sections[sec][2]+len(self.sections[sec][0].x)-1]
        return prox_comps, dist_comps  

    def calc_LayerDegree(self,segm_range=[0],layer=-1,degree=-1,dist_limit=[],diam_limit=[]):
        ''' calculate section name, proximal and distal compartments of dendrites(sections) in a layer and degree,
            layer=-1 means all layers(1-n), degree=-1 means all degrees(0-m). 
            dist_limit is limit of distance of a section, diam_limit is limit of diameter of a section. by huanwen chen 
        '''
        secs_name = []
        prox_comps = []
        dist_comps = []
        for sec in self.sections:
            sn = []
            pcNr = []
            dcNr = []
            if (str(sec) != 'soma') and (str(sec) != 'axon') and (self.sections[sec][2] in segm_range):
                if (layer == -1) or (len(sec) - len('soma') == layer):
                    layerID = True
                else:				
                    layerID = False
                if (degree == -1) or (len(self.sections[sec][0].children) == degree):
                    degreeID = True
                else:				
                    degreeID = False				
                if (dist_limit == []) or (dist_limit[0] <= self.distances[sec][0] <= dist_limit[1]):				
                    distID = True
                else:				
                    distID = False				
                if (diam_limit == []) or (diam_limit[0] <= self.sections[sec][0].diameter[0]/um <= diam_limit[1]):				
                    diamID = True				
                else:
                    diamID = False

                if layerID and degreeID and distID and diamID:
                    sn = [sec]
                    pcNr = [self.sections[sec][2]]
                    dcNr = [self.sections[sec][2]+len(self.sections[sec][0].x)-1]				

            secs_name += sn
            prox_comps += pcNr
            dist_comps += dcNr
        return secs_name,prox_comps,dist_comps
														
    def show_shape3d(self, fov=400, fig=None, ax=None):
        """Show a 3D plot of the morphology"""
        if fig is None:
            fig = plt.figure()
            ax = fig.add_subplot(111,projection='3d')
            # Set fov
            ax.set_xlim(-fov, fov)
            ax.set_ylim(-fov, fov)
            ax.set_zlim(-fov, fov)
            ax.set_aspect("equal")

#        data = self.sections["soma"][1]
        for sec in self.sections.values():
            self.lines = add_line3d(ax, self.lines, sec[0])
            
#        cmap = plt.get_cmap('Spectral')
#        cmap = plt.get_cmap('summer')
##        print np.max(data)
#        cNorm = colorz.Normalize(vmin=0, vmax=1)
#        scalarMap = clrm.ScalarMappable(norm=cNorm, cmap=cmap)
#        for sec in self.sections:
#            for ii in range(len(self.lines[sec])):
##                print ii
#                vsec = self.sections[sec][1][ii]
#                cval = scalarMap.to_rgba(vsec)
#                self.lines[sec][ii].set_color(cval)

            
    def show_segm(self, fov=500, fig=None, ax=None, segm_range = 0,colorv = 'r', segm_range2 = [0],colorv2 = 'k'):
        """Show a 2D plot of the morphology, highlight sections in range 'segm_range' """
        if fig is None:
            fig, ax = plt.subplots()
            # Set fov
            ax.set_xlim(-240, 110)
            ax.set_ylim(-200, 350)
            ax.set_aspect("equal")
        
        #add lines
        for sec in self.sections.values():
            self.lines = add_line(ax, self.lines, sec[0])
        
        # change colors
        for sec in self.sections:
            for ii in range(len(self.sections[sec][0].x)):
                if (self.sections[sec][2]+ii in segm_range):
                    cval = colorv
                    self.lines[sec][ii].set_linewidth(2)   #3
                elif (self.sections[sec][2]+ii in segm_range2):
                    cval = colorv2
                    self.lines[sec][ii].set_linewidth(2)   #3
                else:
                    cval = 'black'
                self.lines[sec][ii].set_color(cval)
                
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis='both',which='both',bottom='off',left='off',right='off',  top='off', labelbottom='off', labelleft='off')

#        savefig('./'+'Neuron'+'.eps', format='eps', dpi=1000)
#        fig.suptitle(str(int(distMin))+'um to '+str(int(distMax))+' um')

    def show_segm_3ROI(self, fov=500, fig=None, ax=None, linewdbasic=1.5,segm_range1 = [],colorv1 = 'k',widthv1 = 0, segm_range2 = [],colorv2 = 'k',widthv2 = 0, segm_range3 = [],colorv3 = 'k',widthv3 = 0,widthv_=0):
        """Show a 2D plot of the morphology, highlight sections in 3 ROIs. by huanwen chen """
        if fig is None:
            fig, ax = plt.subplots()
            # Set fov
            ax.set_xlim(-240, 110)
            ax.set_ylim(-200, 350)
            ax.set_aspect("equal")
        
        #add lines
        for sec in self.sections.values():
            self.lines = add_line(ax, self.lines, sec[0],linewdbasic=linewdbasic)
        
        # change colors
        for sec in self.sections:
            for ii in range(len(self.sections[sec][0].x)):
                if (self.sections[sec][2]+ii in segm_range1):
                    cval = colorv1
                    wval = widthv1
                elif (self.sections[sec][2]+ii in segm_range2):
                    cval = colorv2
                    wval = widthv2
                elif (self.sections[sec][2]+ii in segm_range3):
                    cval = colorv3
                    wval = widthv3
                else:
                    cval = 'black'
                    wval = widthv_
                self.lines[sec][ii].set_color(cval)
                if (wval > 0):
                    self.lines[sec][ii].set_linewidth(wval)
                
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis='both',which='both',bottom='off',left='off',right='off',  top='off', labelbottom='off', labelleft='off')

#        savefig('./'+'Neuron'+'.eps', format='eps', dpi=1000)
#        fig.suptitle(str(int(distMin))+'um to '+str(int(distMax))+' um')

    def show_segm_byWidth(self, fov=500, fig=None, ax=None, linewdbasic=0.5, idle_weight_th=0.5, segm_range1=[],width_to_show1=[],colorv1='k', segm_range2=[],width_to_show2=[],colorv2='k', segm_range3=[],width_to_show3=[],colorv3='k', widthv_=0):
        """Show a 2D plot of the morphology, highlight sections in 3 ROIs. by huanwen chen """
        if fig is None:
            fig, ax = plt.subplots()
            # Set fov
            ax.set_xlim(-240, 110)
            ax.set_ylim(-200, 350)
            ax.set_aspect("equal")

        segm_range1 = list(segm_range1)
        width_to_show1 = list(width_to_show1)
        segm_range2 = list(segm_range2)
        width_to_show2 = list(width_to_show2)		
        segm_range3 = list(segm_range3)
        width_to_show3 = list(width_to_show3)
		
        #add lines
        for sec in self.sections.values():
            self.lines = add_line(ax, self.lines, sec[0],linewdbasic=linewdbasic)
        
        # change colors
        for sec in self.sections:
            for ii in range(len(self.sections[sec][0].x)):
                currSegm = self.sections[sec][2]+ii
                if (currSegm in segm_range1):
                    segmindex = segm_range1.index(currSegm)
                    wval = width_to_show1[segmindex]
#                    if wval>linewdbasic:
#                        cval = colorv1
#                    else:
#                        cval = 'black'
#                        wval = linewdbasic
                    cval = colorv1
                    if wval < idle_weight_th:
                        wval = linewdbasic
                        cval = 'black'
                    else:
                        wval = 1.0*wval
                elif (currSegm in segm_range2):
                    segmindex = segm_range2.index(currSegm)
                    wval = width_to_show2[segmindex]
#                    if wval>linewdbasic:
#                        cval = colorv2
#                    else:
#                        cval = 'black'
#                        wval = linewdbasic
                    cval = colorv2
                    if wval < idle_weight_th:
                        wval = linewdbasic
                        cval = 'black'
                    else:
                        wval = 1.0*wval
                elif (currSegm in segm_range3):
                    segmindex = segm_range3.index(currSegm)
                    wval = width_to_show3[segmindex]
#                    if wval>linewdbasic:
#                        cval = colorv3
#                    else:
#                        cval = 'black'
#                        wval = linewdbasic
                    cval = colorv3
                    if wval < idle_weight_th:
                        wval = linewdbasic
                        cval = 'black'
                    else:
                        wval = 1.0*wval
                else:
                    cval = 'black'
                    wval = widthv_
                self.lines[sec][ii].set_color(cval)
                if (wval > 0):
                    self.lines[sec][ii].set_linewidth(wval)
                
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis='both',which='both',bottom='off',left='off',right='off',  top='off', labelbottom='off', labelleft='off')
                
    def show_segm_byName(self, fov=500, fig=None, ax=None, segmName='soma'):
        """Show a 2D plot of the morphology, highlight section with name 'segmName' """
        if fig is None:
            fig, ax = plt.subplots()
            # Set fov
            ax.set_xlim(-fov, fov)
            ax.set_ylim(-fov, fov)
            ax.set_aspect("equal")

        for sec in self.sections.values():
            self.lines = add_line(ax, self.lines, sec[0])
            
        for sec in self.sections:
            for ii in range(len(self.sections[sec][0].x)):
                if (str(sec) == segmName):
                    cval = 'red'
                    print (self.sections[sec][2]+ii)
                    self.lines[sec][ii].set_linewidth(3)
                else:
                    cval = 'black'
                self.lines[sec][ii].set_color(cval)
                
            
    def show_shape(self, fovx=200,fovy=200, fig=None, ax=None):
        """Show a 2D plot of the morphology"""
        if fig is None:
            fig, ax = plt.subplots()
            # Set fov
            ax.set_xlim(-fovx, fovx)
            ax.set_ylim(-fovy, fovy)
            ax.set_aspect("equal")

        for sec in self.sections.values():
            self.lines = add_line(ax, self.lines, sec[0])
            
#        cmap = plt.get_cmap('spectral')
#        print np.max(data)
#        cNorm = colorz.Normalize(vmin=0, vmax=1)
#        scalarMap = clrm.ScalarMappable(norm=cNorm, cmap=cmap)
#        for sec in self.sections:
#            for ii in range(len(self.lines[sec])):
##                print ii
#                vsec = self.sections[sec][1][ii]
#                cval = scalarMap.to_rgba(vsec)
#                self.lines[sec][ii].set_color(cval)

    def animate3d(self):
        """ Make an animation (3D) """
        try:
            nf = len(self.sections["soma"][1][0])
            print( 'nf: '+str(nf))
        except TypeError:
            print( "running simulation first")
            self.run()
            nf = len(self.sections["soma"][1][0])

        data = self.sections["soma"][1][0]
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        ax.set_xlim(-250, 250)
        ax.set_ylim(-250, 250)
        ax.set_zlim(-250, 250)
        ax.set_aspect('equal')
        self.show_shape3d(fig=fig, ax=ax)

        cmap = plt.get_cmap('afmhot')
        cNorm = colorz.Normalize(vmin=-75*mV, vmax=-0*mV)
        scalarMap = clrm.ScalarMappable(norm=cNorm, cmap=cmap)

        def anim3d(i):
            for sec in self.sections:
                for ii in range(len(self.lines[sec])):
#                    print ii
                    vsec = self.sections[sec][1][ii][i]
                    cval = scalarMap.to_rgba(vsec)
                    self.lines[sec][ii].set_color(cval)
            return self.lines,

        # call the animator.
        anim3d = animation.FuncAnimation(fig, anim3d, interval=20,
                                       frames=nf)
        mywriter = animation.FFMpegWriter()
        anim3d.save('basic_animation.avi', fps=30,writer=mywriter)
        
    def animate(self):
        """ Make an animation (2D) """
        try:
            nf = len(self.sections["soma"][1][0])
            print( 'nf: '+str(nf))
        except TypeError:
            print( "running simulation first")
            self.run()
            nf = len(self.sections["soma"][1][0])

        data = np.zeros([1,nf])
        for x in self.sections.values():
            data = np.append(data,np.array(x[1]),axis=0)
        data = data[1:,:]
        fig, ax = plt.subplots()
        ax.set_xlim(-500, 500)
        ax.set_ylim(-500, 500)
        ax.set_aspect('equal')
        self.show_shape(fig=fig, ax=ax)

        cmap = plt.get_cmap('afmhot')
        cNorm = colorz.Normalize(vmin=np.amin(data), vmax=np.amax(data))
        scalarMap = clrm.ScalarMappable(norm=cNorm, cmap=cmap)
        def anim(i):
            for sec in self.sections:
                for ii in range(len(self.lines[sec])):
#                    print ii
                    vsec = self.sections[sec][1][ii][i]
                    cval = scalarMap.to_rgba(vsec)
                    self.lines[sec][ii].set_color(cval)
            return self.lines,

        # call the animator.
        anim = animation.FuncAnimation(fig, anim, interval=20,
                                       frames=nf)
        mywriter = animation.FFMpegWriter()
        anim.save('basic_animation.avi', fps=30,writer=mywriter)
        
    def show_property(self, var_to_show,sspike=1,nspike=1,sspikemin=0,nspikemin=0):
        """Show a 2D plot of the morphology, highlight the distribution of a parameter 'var_to_show' """
        data = var_to_show
        
        fig, ax = plt.subplots()
        ax.set_xlim(-250, 75)
        ax.set_ylim(-200, 150)
        ax.set_aspect('equal')
#        ax.set_axis_bgcolor((0.93,0.93,0.93))
        self.show_shape(fig=fig, ax=ax)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis='both',which='both',bottom='off',left='off',right='off',  top='off', labelbottom='off', labelleft='off')

        minmin = np.amin(data)
        maxmax = np.amax(data)
                
#        cmap = plt.get_cmap('seismic') 
#        cNorm = colorz.Normalize(vmin=minmin, vmax= maxmax ) # np.amin(data)  np.amax(data)
#        scalarMap = clrm.ScalarMappable(norm=cNorm, cmap=cmap)
#        scalarMap.set_array([minmin,maxmax])
        
        orig_cmap = clrm.coolwarm #coolwarm, seismic,bwr,rainbow, jet
        shiftc = 1 - maxmax/(maxmax + abs(minmin))
        newcmap = shiftedColorMap(orig_cmap, start=0, midpoint=shiftc,
                                  stop=1.0, name='shiftedcmap',
                                  somspike=sspike,nmdaspike=nspike)
              
        cNorm = colorz.Normalize(vmin=minmin, vmax= maxmax )
        scalarMap = clrm.ScalarMappable(norm=cNorm,cmap=newcmap)
        scalarMap.set_array([minmin,maxmax])
        
        cbar = plt.colorbar(scalarMap,ticks = np.arange(minmin,maxmax+1))
        lblvar = list(range(sspikemin,sspikemin+sspike))+[' ']+list(range(nspikemin+nspike-1,nspikemin-1,-1))
        cbar.ax.set_yticklabels(lblvar)
        
        for sec in self.sections:
            sec_nr = self.sections[sec][2]
            for ii in range(len(self.lines[sec])):
                vsec = var_to_show[sec_nr+ii]
                cval = scalarMap.to_rgba(vsec)
                self.lines[sec][ii].set_color(cval)
        
#        title('Location-dependence evoking spikes',fontsize=25)
#        text(-350,150,'NMDA spike',color='r',fontsize=20)
#        text(-350,100,'Somatic spike',color='b',fontsize=20)
#        axis('off')  
        
    def show_segm_property(self,fig=None,ax=None,linewdbasic=1.5,cbaron=True,segm_range=[],var_to_show=[],widthv=0,segm_range1=[],widthv1=0,colorv_='k',widthv_=0,midv=100.,sspike=5,nspike=5):
        """ Show a 2D plot of the morphology, highlight the distribution of a parameter 'var_to_show' in range 'segm_range' by huanwen chen """
        if fig is None:
            fig, ax = plt.subplots()
            #ax.set_xlim(-240, 110)
            #ax.set_ylim(-200, 350)
            ax.set_aspect("equal")
		
		#add lines
        for sec in self.sections.values():
            self.lines = add_line(ax, self.lines, sec[0],linewdbasic=linewdbasic)

        data = np.array(var_to_show)
        var_list = list(var_to_show)
        segm_list = list(segm_range)
        if (len(segm_list) != len(var_list)):
            return			

        minmin = np.amin(data)
        maxmax = np.amax(data)
        if maxmax <= midv :
            orig_cmap = clrm.coolwarm   #winter,seismic,bwr,rainbow, jet,Blues
            maxmax = midv   
            startpoint = 0
            shiftc = 1
            stoppoint = 0.5
        elif minmin >= midv :
            orig_cmap = clrm.coolwarm   #winter,seismic,bwr,rainbow, jet,Blues
            minmin = midv
            startpoint = 0.5
            shiftc = 0
            stoppoint = 1            
        else:
            orig_cmap = clrm.coolwarm #coolwarm, seismic,bwr,rainbow, jet
            startpoint = 0
            stoppoint = 1   
            nrsn = sspike + nspike
            sval = (midv-minmin)/(maxmax-minmin)
            nval = 1.0 - sval
            shiftc = sval - 0.5/nrsn
            sspike = int(round(nrsn * sval))
            nspike = nrsn - sspike
            if sspike < 1:
                sspike = 1
                nspike = nrsn - sspike
            elif nspike < 1:
                nspike = 1
                sspike = nrsn - nspike

        #shiftc = 1 - maxmax/(maxmax + abs(minmin))

        newcmap = shiftedColorMap(orig_cmap, start=startpoint, midpoint=shiftc,
                                  stop=stoppoint, name='shiftedcmap',
                                  somspike=sspike,nmdaspike=nspike)		
        cNorm = colorz.Normalize(vmin=minmin, vmax= maxmax )
        scalarMap = clrm.ScalarMappable(norm=cNorm,cmap=newcmap)   #cmap=newcmap,cmap=orig_cmap
        scalarMap.set_array([minmin,maxmax])        
        cbar = plt.colorbar(scalarMap,ax=ax,shrink=0.6,pad=0.03)   #,ticks = np.arange(minmin,maxmax+1)
        plt.setp(cbar.ax.get_yticklabels(),fontsize=6,fontname='Times New Roman',color=(0.5,0.5,0.5))
        cbar.ax.tick_params(direction='in',length=2)			
        cbar.ax.set_visible(cbaron)

		# change colors
        for sec in self.sections:
            for ii in range(len(self.sections[sec][0].x)):
                segm_nr = self.sections[sec][2]+ii
                if (segm_nr in segm_list):
                    segm_index = segm_list.index(segm_nr)
                    vsegm = var_list[segm_index]                    
                    cval = scalarMap.to_rgba(vsegm)
                    if (segm_nr in segm_range1):
                        wval = widthv1
                    else:
                        wval = widthv
                else:
                    cval = colorv_
                    if (segm_nr in segm_range1):
                        wval = widthv1
                    else:
                        wval = widthv_
                self.lines[sec][ii].set_color(cval)
                if (wval > 0):
                    self.lines[sec][ii].set_linewidth(wval)
                
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis='both',which='both',bottom='off',left='off',right='off',  top='off', labelbottom='off', labelleft='off')

#        savefig('./'+'Neuron'+'.eps', format='eps', dpi=1000)
#        fig.suptitle(str(int(distMin))+'um to '+str(int(distMax))+' um')        

    def show_segm_property1(self,fig=None,ax=None,linewdbasic=1.5,cbaron=True,segm_range=[],var_to_show=[],widthv=0,segm_range1=[],widthv1=0,colorv_='k',widthv_=0,midv=100.):
        """ Show a 2D plot of the morphology, highlight the distribution of a parameter 'var_to_show' in range 'segm_range' by huanwen chen """
        if fig is None:
            fig, ax = plt.subplots()
            #ax.set_xlim(-240, 110)
            #ax.set_ylim(-200, 350)
            ax.set_aspect("equal")
		
		#add lines
        for sec in self.sections.values():
            self.lines = add_line(ax, self.lines, sec[0],linewdbasic=linewdbasic)

        data = np.array(var_to_show)
        var_list = list(var_to_show)
        segm_list = list(segm_range)
        if (len(segm_list) != len(var_list)):
            return			

        minmin = np.amin(data)
        maxmax = np.amax(data)
        oldcolors = clrm.get_cmap('coolwarm', 256)   #coolwarm,winter,seismic,bwr,rainbow, jet,Blues
        newcolors = oldcolors(np.linspace(0, 1, 256))
        black = np.array([0, 0, 0, 1])
        newcolors[128-10:128+10, :] = black
        orig_cmap = ListedColormap(newcolors)
        if maxmax <= midv :
            maxmax = midv   
            startpoint = 0
            shiftc = 1
            stoppoint = 0.5
        elif minmin >= midv :
            minmin = midv
            startpoint = 0.5
            shiftc = 0
            stoppoint = 1            
        else:
            startpoint = 0
            stoppoint = 1   
            shiftc = (midv-minmin)/(maxmax-minmin)

        newcmap = modifiedColorMap(orig_cmap, start=startpoint, midpoint=shiftc,
                                  stop=stoppoint, name='modifiededcmap')		
        cNorm = colorz.Normalize(vmin=minmin, vmax= maxmax )
        scalarMap = clrm.ScalarMappable(norm=cNorm,cmap=newcmap)   #cmap=newcmap,cmap=orig_cmap
        scalarMap.set_array([minmin,maxmax])        
        cbar = plt.colorbar(scalarMap,ax=ax,shrink=0.6,pad=0.03)   #,ticks = np.arange(minmin,maxmax+1)
        plt.setp(cbar.ax.get_yticklabels(),fontsize=6,fontname='Times New Roman',color=(0.5,0.5,0.5))
        cbar.ax.tick_params(direction='in',length=2)			
        cbar.ax.set_visible(cbaron)

		# change colors
        for sec in self.sections:
            for ii in range(len(self.sections[sec][0].x)):
                segm_nr = self.sections[sec][2]+ii
                if (segm_nr in segm_list):
                    segm_index = segm_list.index(segm_nr)
                    vsegm = var_list[segm_index]                    
                    cval = scalarMap.to_rgba(vsegm)
                    if (segm_nr in segm_range1):
                        wval = widthv1
                    else:
                        wval = widthv
                else:
                    cval = colorv_
                    if (segm_nr in segm_range1):
                        wval = widthv1
                    else:
                        wval = widthv_
                self.lines[sec][ii].set_color(cval)
                if (wval > 0):
                    self.lines[sec][ii].set_linewidth(wval)
                
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis='both',which='both',bottom='off',left='off',right='off',  top='off', labelbottom='off', labelleft='off')

    def property_showing(self,fig=None,ax=None,linewdbasic=1.5,cbaron=True,segm_range=[],var_to_show=[],widthv=0,segm_range1=[],widthv1=0,colorv_='k',widthv_=0,midv=100.,blackv=0.05):
        """ Show a 2D plot of the morphology, highlight the distribution of a parameter 'var_to_show' in range 'segm_range' by huanwen chen """
        if fig is None:
            fig, ax = plt.subplots()
            #ax.set_xlim(-240, 110)
            #ax.set_ylim(-200, 350)
            ax.set_aspect("equal")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis='both',which='both',bottom='off',left='off',right='off',  top='off', labelbottom='off', labelleft='off')		
		#add lines
        for sec in self.sections.values():
            self.lines = add_line(ax, self.lines, sec[0],linewdbasic=linewdbasic)

        data = np.array(var_to_show)
        var_list = list(var_to_show)
        segm_list = list(segm_range)
        if (len(segm_list) != len(var_list)):
            return			

        minmin = np.amin(data)
        maxmax = np.amax(data)
        if  minmin == maxmax:
            return
			
        #oldcolors = clrm.get_cmap('coolwarm', 256)   #coolwarm,winter,seismic,bwr,rainbow, jet,Blues
        #newcolors = oldcolors(np.linspace(0, 1, 256))
        black = np.array([0, 0, 0, 1])
        top = clrm.get_cmap('winter')
        bottom = clrm.get_cmap('Reds')
        ind1 = int(round(256*(midv*(1.0-blackv)-minmin)/(maxmax-minmin)))
        if ind1 < 0:
            ind1 = 0
        ind2 = int(round(256*(midv*(1.0+blackv)-minmin)/(maxmax-minmin)))
        if ind2 > 256:
            ind2 = 256    
        print(ind1,',',ind2)
        newcolors = np.vstack((top(np.linspace(0, 1, ind2)),
                   bottom(np.linspace(0, 1, 256-ind1))))
        newcolors[ind1:ind2, :] = black

        cMap = ListedColormap(newcolors)			
        cNorm = colorz.Normalize(vmin=minmin, vmax= maxmax ) # np.amin(data)  np.amax(data)
        scalarMap = clrm.ScalarMappable(norm=cNorm, cmap=cMap)
        scalarMap.set_array([minmin,maxmax])
        cbar = plt.colorbar(scalarMap,ax=ax,shrink=0.6,pad=0.03)   #,ticks = np.arange(minmin,maxmax+1)
        plt.setp(cbar.ax.get_yticklabels(),fontsize=6,fontname='Times New Roman',color=(0.5,0.5,0.5))
        cbar.ax.tick_params(direction='in',length=2)			
        cbar.ax.set_visible(cbaron)

		# change colors
        for sec in self.sections:
            for ii in range(len(self.sections[sec][0].x)):
                segm_nr = self.sections[sec][2]+ii
                if (segm_nr in segm_list):
                    segm_index = segm_list.index(segm_nr)
                    vsegm = var_list[segm_index]                    
                    cval = scalarMap.to_rgba(vsegm)
                    if (segm_nr in segm_range1):
                        wval = widthv1
                    else:
                        wval = widthv
                else:
                    cval = colorv_
                    if (segm_nr in segm_range1):
                        wval = widthv1
                    else:
                        wval = widthv_
                self.lines[sec][ii].set_color(cval)
                if (wval > 0):
                    self.lines[sec][ii].set_linewidth(wval)

    def property_showing_chw(self,fig=None,ax=None,linewdbasic=1.5,cbaron=True,segm_range=[],var_to_show=[],widthv=0,segm_range1=[],widthv1=0,colorv_='k',widthv_=0):
        """ Show a 2D plot of the morphology, highlight the distribution of a parameter 'var_to_show' in range 'segm_range' by huanwen chen """
        if fig is None:
            fig, ax = plt.subplots()
            #ax.set_xlim(-240, 110)
            #ax.set_ylim(-200, 350)
            ax.set_aspect("equal")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis='both',which='both',bottom='off',left='off',right='off',  top='off', labelbottom='off', labelleft='off')		
		#add lines
        for sec in self.sections.values():
            self.lines = add_line(ax, self.lines, sec[0],linewdbasic=linewdbasic)

        data = np.array(var_to_show)
        var_list = list(var_to_show)
        segm_list = list(segm_range)
        if (len(segm_list) != len(var_list)):
            return			

        minmin = np.amin(data)
        maxmax = np.amax(data)
        if  minmin == maxmax:
            return
			
        cmap = plt.get_cmap('seismic')   #coolwarm,winter,seismic,bwr,rainbow, jet,Blues
        cNorm = colorz.Normalize(vmin=minmin, vmax= maxmax ) # np.amin(data)  np.amax(data)
        scalarMap = clrm.ScalarMappable(norm=cNorm, cmap=cmap)
        scalarMap.set_array([minmin,maxmax])
        cbar = plt.colorbar(scalarMap,ax=ax,shrink=0.6,pad=0.03)   #,ticks = np.arange(minmin,maxmax+1)
        plt.setp(cbar.ax.get_yticklabels(),fontsize=6,fontname='Times New Roman',color=(0.5,0.5,0.5))
        cbar.ax.tick_params(direction='in',length=2)			
        cbar.ax.set_visible(cbaron)

		# change colors
        for sec in self.sections:
            for ii in range(len(self.sections[sec][0].x)):
                segm_nr = self.sections[sec][2]+ii
                if (segm_nr in segm_list):
                    segm_index = segm_list.index(segm_nr)
                    vsegm = var_list[segm_index]                    
                    cval = scalarMap.to_rgba(vsegm)
                    if (segm_nr in segm_range1):
                        wval = widthv1
                    else:
                        wval = widthv
                else:
                    cval = colorv_
                    if (segm_nr in segm_range1):
                        wval = widthv1
                    else:
                        wval = widthv_
                self.lines[sec][ii].set_color(cval)
                if (wval > 0):
                    self.lines[sec][ii].set_linewidth(wval)	

    def synapses_showing(self,fig=None,ax=None,linewdbasic=1.5,cbaron=True,segm_range=[],var_to_show=[],midv=100.,blackv=0.05,synR=2.0):
        """ Show a 2D plot of the morphology, highlight the distribution of a parameter 'var_to_show' in range 'segm_range' by huanwen chen """
        if fig is None:
            fig, ax = plt.subplots()
            #ax.set_xlim(-240, 110)
            #ax.set_ylim(-200, 350)
            ax.set_aspect("equal")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis='both',which='both',bottom='off',left='off',right='off',  top='off', labelbottom='off', labelleft='off')		
		#add lines
        for sec in self.sections.values():
            self.lines = add_line(ax, self.lines, sec[0],linewdbasic=linewdbasic)

        data = np.array(var_to_show)
        var_list = list(var_to_show)
        segm_list = list(segm_range)
        if (len(segm_list) != len(var_list)):
            return			

        minmin = np.amin(data)
        maxmax = np.amax(data)
        print(minmin,maxmax)
        if  minmin != maxmax:		
            #oldcolors = clrm.get_cmap('coolwarm', 256)   #coolwarm,winter,seismic,bwr,rainbow, jet,Blues
            #newcolors = oldcolors(np.linspace(0, 1, 256))
            black = np.array([0, 0, 0, 1])
            top = clrm.get_cmap('winter')
            bottom = clrm.get_cmap('Reds')
            if maxmax < midv*(1.0-blackv):
                cMap = clrm.get_cmap('winter')
            elif minmin > midv*(1.0+blackv):
                cMap = clrm.get_cmap('Reds')
            else:
                ind1 = int(round(256*(midv*(1.0-blackv)-minmin)/(maxmax-minmin)))
                if ind1 < 0:
                    ind1 = 0
                ind2 = int(round(256*(midv*(1.0+blackv)-minmin)/(maxmax-minmin)))
                if ind2 > 256:
                    ind2 = 256    
#                print(ind1,',',ind2)
                newcolors = np.vstack((top(np.linspace(0, 1, ind2)),bottom(np.linspace(0, 1, 256-ind1))))
                newcolors[ind1:ind2, :] = black
                cMap = ListedColormap(newcolors)			
            cNorm = colorz.Normalize(vmin=minmin, vmax= maxmax ) # np.amin(data)  np.amax(data)
            scalarMap = clrm.ScalarMappable(norm=cNorm, cmap=cMap)
            scalarMap.set_array([minmin,maxmax])
            cbar = plt.colorbar(scalarMap,ax=ax,shrink=0.6,pad=0.03)   #,ticks = np.arange(minmin,maxmax+1)
            plt.setp(cbar.ax.get_yticklabels(),fontsize=6,fontname='Times New Roman',color=(0.5,0.5,0.5))
            cbar.ax.tick_params(direction='in',length=2)
            cbar.ax.set_visible(cbaron)
        
        rScale = 1.0
        for ii in range(len(segm_list)):
            compNr = segm_list[ii]
            synV = var_list[ii]
            rScale = synV/100.0
            if minmin == maxmax:
#                rScale = 1.0
                cval = 'k'
            else:
#                rScale = synV/100.0
                cval = scalarMap.to_rgba(synV)
            startPt = self.compStart[compNr]
            endPt = self.compEnd[compNr]
            x0 = startPt[0]*1e-3
            y0 = startPt[1]*1e-3
            z0 = startPt[2]*1e-3
            x1 = endPt[0]*1e-3
            y1 = endPt[1]*1e-3
            z1 = endPt[2]*1e-3
            randC = np.random.rand()
            xc = x0 + randC*(x1-x0)
            yc = y0 + randC*(y1-y0)
            synapse = plt.Circle((xc,yc), rScale*synR, color=cval, alpha=1)
            ax.add_artist(synapse)			

    def marker_showing_chw(self, fov=500, fig=None, ax=None, linewdbasic=1.5,lw_=0,segm_range1 = [],mk1 = [],ms1 = 0,mfc1='k',segm_range2 = [],mk2 = [],ms2 = 0,mfc2='k'):
        """Show a 2D plot of the morphology, highlight sections in 3 ROIs. by huanwen chen 2020.04.02"""
        if fig is None:
            fig, ax = plt.subplots()
            # Set fov
            #ax.set_xlim(-240, 110)
            #ax.set_ylim(-200, 350)
            ax.set_aspect("equal")
        
        #add lines
        for sec in self.sections.values():
            self.lines = add_line(ax, self.lines, sec[0],linewdbasic=linewdbasic)
        
        # change colors
        for sec in self.sections:
            for ii in range(len(self.sections[sec][0].x)):
                if (self.sections[sec][2]+ii in segm_range1):
                    mkval = mk1
                    msval = ms1
                    mfcval = mfc1
                elif (self.sections[sec][2]+ii in segm_range2):
                    mkval = mk2
                    msval = ms2
                    mfcval = mfc2
                else:
                    mkval = []
					
                if mkval != []:
                    self.lines[sec][ii].set_marker(mkval)
                    self.lines[sec][ii].set_markersize(msval)
                    self.lines[sec][ii].set_markerfacecolor(mfcval)
                if (lw_ > 0):
                    self.lines[sec][ii].set_linewidth(lw_)
                
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis='both',which='both',bottom='off',left='off',right='off',  top='off', labelbottom='off', labelleft='off')
			
    def color_showing_chw(self,fig=None,ax=None,mappattern='seismic',linewdbasic=1.5,cbaron=True,segm_range=[],var_to_show=[],wInit=0.5,colorv_='k',widthv_=0):
        """ Show a 2D plot of the morphology, highlight the distribution of a parameter 'var_to_show' in range 'segm_range' by huanwen chen 2020.04.02"""
        if fig is None:
            fig, ax = plt.subplots()
            #ax.set_xlim(-240, 110)
            #ax.set_ylim(-200, 350)
            #ax.set_aspect("equal")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis='both',which='both',bottom='off',left='off',right='off',  top='off', labelbottom='off', labelleft='off')
		
        data = np.array(var_to_show)
        var_list = list(var_to_show)
        segm_list = list(segm_range)
        if (len(segm_list) != len(var_list)):
            return			
        minmin = np.amin(data)
        maxmax = np.amax(data)
        if  minmin == maxmax:
            if minmin > wInit:
                self.show_segm_3ROI(fig=fig,ax=ax,linewdbasic=0.5,segm_range1=segm_list,colorv1='r',widthv1=0.5,widthv_=0)
            else:
                self.show_segm_3ROI(fig=fig,ax=ax,linewdbasic=0.5,segm_range1=segm_list,colorv1='b',widthv1=0.5,widthv_=0)
            return	
		#add lines
        for sec in self.sections.values():
            self.lines = add_line(ax, self.lines, sec[0],linewdbasic=linewdbasic)			
        cmap = plt.get_cmap(mappattern)   #coolwarm,winter,seismic,bwr,rainbow, jet,Blues
        cNorm = colorz.Normalize(vmin=minmin, vmax= maxmax ) # np.amin(data)  np.amax(data)
        scalarMap = clrm.ScalarMappable(norm=cNorm, cmap=cmap)
        scalarMap.set_array([minmin,maxmax])
        cbar = plt.colorbar(scalarMap,ax=ax,shrink=0.6,pad=0.03)   #,ticks = np.arange(minmin,maxmax+1)  shrink=0.6,pad=0.03
        plt.setp(cbar.ax.get_yticklabels(),fontsize=6,fontname='Times New Roman',color=(0.5,0.5,0.5))
        cbar.ax.tick_params(direction='out',length=1)			
        cbar.ax.set_visible(cbaron)
		# change colors
        for sec in self.sections:
            for ii in range(len(self.sections[sec][0].x)):
                segm_nr = self.sections[sec][2]+ii
                if (segm_nr in segm_list):
                    segm_index = segm_list.index(segm_nr)
                    vsegm = var_list[segm_index]                    
                    cval = scalarMap.to_rgba(vsegm)
                else:
                    cval = colorv_
                self.lines[sec][ii].set_color(cval)
                if (widthv_ > 0):
                    self.lines[sec][ii].set_linewidth(widthv_)

    def show_multiarea_segm(self,fig=None,ax=None,linewdbasic=1.5,wTh=0.01,inputcomps=[],synapticweights=[],widthv=0,colorv_='k',widthv_=0):
        """Show a 2D plot of the morphology, highlight sections in multi-area. by huanwen chen 2020.04.25"""
        #colorvs = ['m','c','g','y','olive','brown','orange','skyblue','springgreen','darkblue']
        mpl.style.use('default')   # mpl.style.use('seaborn')
        incomps = [[] for _ in range(len(inputcomps))]
        for ii in range(len(inputcomps)):
            for jj in range(len(inputcomps[ii])):
                if synapticweights[ii][jj] > wTh:
                    incomps[ii].append(inputcomps[ii][jj])
        if fig is None:
            fig, ax = plt.subplots()
            #ax.set_xlim(-240, 110)
            #ax.set_ylim(-200, 350)
            #ax.set_aspect("equal")
		#add lines of pyramidal neuron
        for sec in self.sections.values():
            self.lines = add_line(ax, self.lines, sec[0],linewdbasic=linewdbasic)
        for sec in self.sections:	
            for ii in range(len(self.sections[sec][0].x)):
                connprop = []
                for jj in range(len(incomps)):
                    if (self.sections[sec][2]+ii in incomps[jj]):
                        connprop.append(jj)
                if len(connprop) == 1:
                    cval = 'C'+str(connprop[0]+1)
                    wval = widthv		
                elif len(connprop) > 1:
                    cval = 'C0'
                    wval = widthv		
                else:
                    cval = colorv_
                    wval = widthv_		
                self.lines[sec][ii].set_color(cval)
                if (wval > 0):
                    self.lines[sec][ii].set_linewidth(wval)						
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis='both',which='both',bottom='off',left='off',right='off',  top='off', labelbottom='off', labelleft='off')	
					
    def multitasks_synaptic_connection(self,fig=None,ax=None,xylimits=[],linewdbasic=0.5,wTh=0.01,wInit=0.5,inputcomps=[],synapticweights=[],widthv=0,colorv_='k',widthv_=0):
        """ 显示多个神经元群（多任务）与锥体神经元的突触连接强度。连接强度由突触权值表达。 by huanwen chen 2020.04.25"""	
        mpl.style.use('default')   # mpl.style.use('seaborn')		
        incomps = deepcopy(inputcomps)
        synweights = deepcopy(synapticweights)
        if incomps==[] or synweights==[]:
            print('inputcomps or synapticweights can not be empty!')
            return
        if fig is None:
            fig, ax = plt.subplots()
            #ax.set_xlim(-240, 110)
            #ax.set_ylim(-200, 350)
            #ax.set_aspect("equal")
		#add lines of pyramidal neuron
        self.show_multiarea_segm(fig=fig,ax=ax,linewdbasic=linewdbasic,wTh=wTh,inputcomps=incomps,synapticweights=synweights,widthv=widthv,colorv_=colorv_,widthv_=widthv_)
        #add lines of synaptic connection
        if xylimits == []:
            xlimits = ax.get_xlim()
            ylimits = ax.get_ylim()
        else:
            xlimits = xylimits[0]
            ylimits = xylimits[1]			
        inputx = xlimits[0]+0.1*abs(xlimits[1]-xlimits[0])
        inputys = np.linspace(ylimits[1]-0.1*abs(ylimits[1]-ylimits[0]),ylimits[0]+0.1*abs(ylimits[1]-ylimits[0]),len(incomps))
        idiam = 0.02*abs(xlimits[1]-xlimits[0])
        for ii in range(len(incomps)):
            segm_range = incomps[ii]
            weight_range = synweights[ii]
            ix = inputx
            iy = inputys[ii]
            circle = plt.Circle((ix, iy), idiam, color='C'+str(ii+1), alpha=0.5)
            ax.add_artist(circle)
            statetext = str(ii+1)
            plt.text(ix-1.2*idiam,iy,statetext,fontsize=6, va="center", ha="right", multialignment="left",color='black')   
            for sec in self.sections:
                for jj in range(len(self.sections[sec][0].x)):
                    if (self.sections[sec][2]+jj in segm_range):
                        segind = segm_range.index(self.sections[sec][2]+jj)
                        sx = np.asarray(self.sections[sec][0].x)*1e6
                        sy = np.asarray(self.sections[sec][0].y)*1e6
                        if weight_range[segind] >= wInit:
                            #ax.plot([ix,sx[jj]],[iy,sy[jj]],'-k',alpha=0.2,linewidth=0.5*weight_range[segind])
                            ax.plot([ix,sx[jj]],[iy,sy[jj]],'-',color='darkgray',alpha=0.2,linewidth=0.2)
                        elif (weight_range[segind] < wInit) and (weight_range[segind] > wTh):
                            #ax.plot([ix,sx[jj]],[iy,sy[jj]],'--k',alpha=0.2,linewidth=0.5*weight_range[segind])
                            ax.plot([ix,sx[jj]],[iy,sy[jj]],'--',color='darkgray',dashes=[idiam,0.5*idiam],alpha=0.2,linewidth=0.2)
                        del segm_range[segind]
                        del weight_range[segind]
        weightStat = np.zeros((len(synapticweights),2))
        countStat = np.zeros((len(synapticweights),2),dtype=int)
        for ii in range(len(synapticweights)):
            a = np.array(synapticweights[ii])
            a0 = a[a>=wInit]
            a1 = a[a<wInit]
            a1 = a1[a1>wTh]
            if len(a0) == 0:
                weightStat[ii,0] = 0
            else:
                weightStat[ii,0] = np.mean(a0)
            if len(a1) == 0:
                weightStat[ii,1] = 0
            else:
                weightStat[ii,1] = np.mean(a1)
            countStat[ii,0] = len(a0)
            countStat[ii,1] = len(a1)
        print(weightStat)
        print(countStat)
        data2show = ''
        for ii in range(len(synapticweights)):
            data2show = data2show+str(ii+1)+': '+'PW = '+str(np.round(weightStat[ii,0],decimals=4))+' , '+'PN = '+str(countStat[ii,0])+' , '  \
                         +'DW = '+str(np.round(weightStat[ii,1],decimals=4))+' , '+'DN = '+str(countStat[ii,1])+'.\n'
        ax.set_xlabel(data2show,fontsize=6,ha='right',ma='left',position=(0.9,0),color='gray')		
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis='both',which='both',bottom='off',left='off',right='off',  top='off', labelbottom='off', labelleft='off')			
		
    def v_record(self, neuron):
        """Set a monitor for the voltage of all segments"""
        return br2.StateMonitor(neuron, 'v', record=True)
        
    def run(self,morphodata):
        """run """
        
        # Set Biophysics
        neuron = self.biophysics(morphodata)
        # Record in every section
        monitor = self.v_record(neuron)
        
        # Initialize the model
        neuron.v = EL        
        neuron.I = 0.*nA
        neuron.I[0] = 0.2*nA
        
        #run
        br2.run(.5*ms)
        
        #store data in sec[1]
        for sec in self.sections.values():  
            kk = sec[2]
            sec[1] = monitor.v[kk:kk+len(sec[0].x)]
#            sec[1] = monitor.gKL[kk:kk+len(sec[0].x)]
        return monitor, neuron

if __name__ == "__main__":
#    test_MDL = '../0. Model/Branco2010_Morpho.swc'
#    morphodata = BrancoData    
#    distal_compartments = distal_compartments_Branco
#    proximal_compartments = proximal_compartments_Branco

    test_MDL = '../0. Model/Acker2008.swc'
    morphodata = AckerData
    distal_compartments = distal_compartments_Acker_eff
    proximal_compartments = proximal_compartments_Acker
    
    test_model = BRIANModel(test_MDL)
    
    
#    M, nrn = test_model.run(morphodata)
#    test_model.show_shape(fovx=200,fovy=500)
#    test_model.show_segm(segm_range=range(32,55,1))
    
#    dv,av = test_model.save_dist_vs_nr(1181)
#    tvect,div = test_model.calc_distCompartments(morphodata['basal'],dv)
    
#    test_model.show_property([5000]+list(morphodata['axon'])+list(morphodata['basal'])+list(morphodata['apical']))
    ii = 16
#    test_model.show_segm(segm_range=range(proximal_compartments_Acker[ii],
#        proximal_compartments_Acker[ii]+3))#
    test_model.show_segm(segm_range=distal_compartments,colorv='r',segm_range2=proximal_compartments,colorv2='b')
    
