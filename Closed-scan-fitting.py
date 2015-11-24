# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 08:42:33 2015

@author: Kaspar Martin Snashall

script designed to normalise and fit z-scan data using stochiastic method
includes two models for fitting v1 and v2

license = MIT
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
from os.path import isfile, join
from scipy import optimize
import sys, re
from matplotlib.widgets import Slider, Button, RadioButtons
from vertical import VertSlider

def ClosedScan(z,z0,delta_phi):
    """uses basic non geometrical equation"""
    z = np.array(z)
    x = (z/z0)
    return 1-4*delta_phi*x/((x**2+9)*(x**2+1))

     
def BeamWaist(Z0):
	#caluculate the beam waist
	wavelength = 1030*10**(-9) # m
	return np.sqrt(abs(Z0)*wavelength/np.pi)

def BeamWaistErr(w0,z0,deltaz0):
    return w0*np.sqrt(deltaz0/z0)

def Intensity(W0,energy,pulsetime):
	#calculate the intensity at focus
	power = energy/pulsetime
	i = power/(4*np.pi*W0**2)
	return i

def IntensityErr(i0,waistErr):
    return i0*waistErr

def L_eff(L,aplha):
	# calculate the effective length
	return (1-np.exp(-aplha*L))/aplha

def Beta(i0,l,q0):
	# calculate the non linear absorbtion coefficent
	return  q0/(i0*l)
 
def BetaErr(b,q0,q0Err,i0Err,i0):
    return b*np.sqrt((q0Err/q0)**2+(i0Err/i0)**2)
 
def ReadAllFiles(folder):
    return [ f for f in listdir(folder) if isfile(join(folder,f)) and '.csv' in f and 'correct' in f]
    
def GetParams(f):
    myl = f.split('_')[1]
    print myl
    if myl == '05nm':
        L = 0.5*10**(-9)
    elif myl == '1nm':
        L = 1*10**(-9)
    elif myl == '5nm':
        L = 5*10**(-9)
    else:
        myl = f.split('_')[2]
        if myl == '05nm':
            L = 0.5*10**(-9)
        elif myl == '1nm':
            L = 1*10**(-9)
        elif myl == '5nm':
            L = 5*10**(-9)
    if L == 1*10**(-9):
        alpha = 506174855.374563722
    elif L == 0.5*10**(-9):
        alpha = 272859143.560142626
    elif L == 5*10**(-9):
        alpha = 215938480.763119764
    else:
        print "unknow sample rewrite this code for new samples"
        sys.exit(0)
    e = f.split('_')[2]
    e2 = re.sub('nJ','',e)
    try:
        energy = float(e2)*10**-9
    except:
        e = f.split('_')[3]
        e2 = re.sub('nJ','',e)
        energy = float(e2)*10**-9
    return energy,L,alpha



    
def Main(i,fitted_data):
    saveDirectory = '/home/kaspar/Desktop/1nmlow'
    fileDirectory = '/home/kaspar/Desktop/1nmlow/'
    pulsetime = 340*10**(-15)     
    
    mybuttonslist = []
    f = ReadAllFiles(fileDirectory)[i]
    data = pd.read_csv(fileDirectory+f) # grab the data
    z = data.iloc[0:,0]/10**3 # z in meters
    ydata = data.iloc[0:,2] #the transmission
    guess = [0.0006,0.25]
    params, params_covariance = optimize.curve_fit(ClosedScan,z, ydata,guess)
    print f
    print params
    print params_covariance
    energy,L,alpah = GetParams(f)
    z0 = params[0]
    q0 = params[1]
    fig, ax = plt.subplots()
    plt.ylabel('Transmission')
    plt.xlabel('z(m)')
    mycolumn = "PtSe$_2$" +str(L*10**9) +"nm closed scan "+str(energy/10**-9)+"nJ"
    plt.title(mycolumn) # title 
    myscat = plt.scatter(z,ydata)
    improvedZ = np.linspace(-0.05,0.05,500)
    myfit, = plt.plot(improvedZ,ClosedScan(improvedZ,params[0],params[1]),color='red')
    center = plt.axvline(x=0, color='black')
    plt.subplots_adjust(right = 0.75, bottom=0.25)
    a0 = q0
    f0 = z0
    
    
    axcolor = 'lightgoldenrodyellow'
    axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
    axamp  = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
    xvals = plt.axes([0.92, 0.25, 0.03, 0.65], axisbg=axcolor)
    yvals = plt.axes([0.85, 0.25, 0.03, 0.65], axisbg=axcolor)  
    xslider = VertSlider(xvals, 'x offset', -0.01, 0.01, valinit=0,valfmt='%.3f')
    yslider = VertSlider(yvals, 'y offset', -0.3, 0.3, valinit=0,valfmt='%.4f')
    
    
    
    myZ0 = Slider(axfreq, 'Z0', 0.00001, 0.05, valinit=f0,valfmt='%.5f')
    myQ0 = Slider(axamp, 'q0', -1, 1, valinit=a0,valfmt='%.5f')
    
    def update(val):
        z0 = myZ0.val
        q0 = myQ0.val
        yshift = yslider.val
        xshift = xslider.val
        shiftedy = ClosedScan(improvedZ,z0,q0)+yshift
        myfit.set_ydata(shiftedy)
        myfit.set_xdata(improvedZ+xshift)
        fig.canvas.draw_idle()
        
    myZ0.on_changed(update)
    myQ0.on_changed(update)
    xslider.on_changed(update)
    yslider.on_changed(update)
    
    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
    def reset(event):
        myZ0.reset()
        myQ0.reset()
        xslider.reset()
        yslider.reset()
    button.on_clicked(reset)
    
    
    
    plt.show(block=False)

    answer = raw_input('Hit enter to continue\n')
    
    z0 = myZ0.val
    q0 = myQ0.val
    yshift = yslider.val
    xshift = xslider.val
    plt.close()
    
    
    plt.ylabel('Transmission')
    plt.xlabel('z(m)')
    mycolumn = "PtSe$_2$" +str(L*10**9) +"nm closed scan "+str(energy/10**-9)+"nJ"
    plt.title(mycolumn) # title 
    myscat = plt.scatter(z,ydata)
    improvedZ = np.linspace(-0.05,0.05,500)
    myfit, = plt.plot(improvedZ+xshift,ClosedScan(improvedZ,z0,q0)+yshift,color='red')
    z0 = z0*10**3
    plt.text(0.8,0.8,"$Z_0 (mm)$ ="+str('%.2f' % z0)+"\n"+" $\Delta \Phi_0$ = "+str('%.2f' % q0)+'\n' ,horizontalalignment='center',
				    verticalalignment='center',
				    transform = ax.transAxes)
    name = re.sub('correct.csv','',f)
    plt.savefig( saveDirectory+f+"closed.png")
    plt.close()
    
    return button
        
if __name__ == '__main__':
    fitted_data = pd.DataFrame(columns = ('Name','L','E(nJ)','Z0(mm)','Zerr','q0','qerr','Adj Z0','Adj q0','delta x','delta y','w0(mm)','werr','I0(Gw/cm^2)','Ierr','alpha','Beta','Beta err'),index=None)
    fileDirectory = '/home/kaspar/Desktop/1nmlow/'
    f = ReadAllFiles(fileDirectory)    
    for i in range(len(f)):
        button = Main(i,fitted_data)
    #fitted_data.to_csv('/home/kaspar/Desktop/fullfitv3.csv',index=False)
        
       
            
    
    
