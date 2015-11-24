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

def OpenScanv3(z,z0,beta):
	X = (z/z0)**2
	m=0 
	mysum = 0 # my sum
	myfactorials = [1,1,2,6,24,120,720,5040,40320,362880,3628800] #factorial numbers up to 10
	while m < 10:
		mysum += ((-beta*leff*Intensity(BeamWaist(z0),energy,pulsetime)*(1/(1+X)))**m)/myfactorials[m]
		m+=1
	return mysum

     
def BeamWaist(Z0):
	#caluculate the beam waist
	wavelength = 1030*10**(-9) # m
	return np.sqrt(abs(Z0)*wavelength/np.pi)

def BeamWaistErr(w0,z0,deltaz0):
    return w0*np.sqrt(deltaz0/z0)

def Intensity(W0,energy,pulsetime):
	#calculate the intensity at focus
	power = energy/pulsetime
	i = power/(np.pi*W0**2)
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
        alpha = 118501220
    elif L == 0.5*10**(-9):
        alpha = 219828946
    elif L == 5*10**(-9):
        alpha = 93780890
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
    saveDirectory = '/home/kaspar/Desktop/20151118-correct/'
    fileDirectory = '/home/kaspar/Desktop/20151118-correct/'
   
    mybuttonslist = []
    f = ReadAllFiles(fileDirectory)[i]
    data = pd.read_csv(fileDirectory+f) # grab the data
    z = data.iloc[0:,0]/10**3 # z in meters
    print f
    ydata = data.iloc[0:,1] #the transmission
    global energy,L,alpha,leff,pulsetime
    
    pulsetime = 340*10**(-15)     
    energy,L,alpha = GetParams(f)
    leff = L_eff(L,alpha)
    print leff
    
    guess= [0.009,0.5]
    params, params_covariance = optimize.curve_fit(OpenScanv3,z, ydata, guess)
    #print params
    #print params_covariance
    sigma1 = params_covariance[0]
    sigma2 = params_covariance[1]
    z0 = params[0]
    b = params[1] # the NL absorbtion coefficent
    if z0 < 0:
        z0 = -z0
    w0 = BeamWaist(z0) # the beam width
    i0 = Intensity(w0,energy,pulsetime) # the intensity at focus        
    
        
    # the one sigma error
    perr = np.sqrt(np.diag(params_covariance))        
    w0Err = BeamWaistErr(w0,z0,perr[0])
    i0Err = IntensityErr(i0,w0Err)
    bErr = perr[1]
    print params
    print perr
    print w0,w0Err
    print i0,i0Err
    print b,bErr
    
    
    fig, ax = plt.subplots()
    plt.ylabel('Transmission')
    plt.xlabel('z(m)')
    
    mycolumn = "PtSe$_2$" +str(L*10**9) +"nm open scan "+str(energy/10**-9)+"nJ"
    plt.title(mycolumn) # title 
    myscat = plt.scatter(z,ydata)
    plt.ylim( ax.get_ylim() )  
    improvedZ = np.linspace(-0.05,0.05,500)
    myfit, = plt.plot(improvedZ,OpenScanv3(improvedZ,params[0],params[1]),color='red')
    #center = plt.axvline(x=0, color='black')
    plt.subplots_adjust(right = 0.75, bottom=0.25)
    a0 = b
    f0 = z0
    
    
    axcolor = 'lightgoldenrodyellow'
    axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
    axamp  = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
    xvals = plt.axes([0.92, 0.25, 0.03, 0.65], axisbg=axcolor)
    yvals = plt.axes([0.85, 0.25, 0.03, 0.65], axisbg=axcolor)  
    xslider = VertSlider(xvals, 'x offset', -0.01, 0.01, valinit=0,valfmt='%.3f')
    yslider = VertSlider(yvals, 'y offset', -0.3, 0.3, valinit=0,valfmt='%.4f')
    
    myZ0 = Slider(axfreq, 'Z0', 0.00001, 0.05, valinit=f0,valfmt='%.5f')
    myQ0 = Slider(axamp, 'beta', -0.000005, 0.000005, valinit=a0,valfmt='%.10f')
    
    def update(val):
        z0 = myZ0.val
        q0 = myQ0.val
        yshift = yslider.val
        xshift = xslider.val
        shiftedy = OpenScanv3(improvedZ,z0,q0)+yshift
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
    b = myQ0.val
    yshift = yslider.val
    xshift = xslider.val
    plt.close()
    if z0 < 0:
        z0 = -z0
   
    w0 = BeamWaist(z0) # the beam width
    i0 = Intensity(w0,energy,pulsetime) # the intensity at focus        
    l = L_eff(L,alpha) # the effective length of the sample
    
   	
         
    # the one sigma error
    perr = np.sqrt(np.diag(params_covariance))        
    w0Err = BeamWaistErr(w0,z0,perr[0])
    i0Err = IntensityErr(i0,w0Err)
    bErr = perr[1]


    print "Z0 = ",z0,"beta = ",b
    print "err0r = ",perr
    print "w0 = ",w0,w0Err
    print "i0 = ",i0,i0Err
    print "b = ", b,bErr
    fig, ax = plt.subplots()
    plt.ylabel('Transmission')
    plt.xlabel('z(m)')
    mycolumn = "PtSe$_2$" +str(L*10**9) +"nm open scan "+str(energy/10**-9)+"nJ"
    plt.title(mycolumn) # title 
    myscat = plt.scatter(z,ydata)
    improvedZ = np.linspace(-0.05,0.05,500)
    myfit, = plt.plot(improvedZ+xshift,OpenScanv3(improvedZ,z0,b)+yshift,color='red')
    z0 = z0*10**3
    w0 = w0*10**6
    i0 = i0/(10**(9)*100*100)
    energy = energy/10**-9
    L = L*10**9
    plt.text(0.8,0.8,"$Z_0 (mm)$ ="+str('%.2f' % z0)+"\n"+'$w_0(\mu m)$ = '+str('%.2f' % w0),horizontalalignment='center',
				    verticalalignment='center',
				    transform = ax.transAxes)
    plt.text(0.2,0.8,"$I_0 (GW cm^{-2}) = $"+str('%.2f' % i0)+"\n"+r"$\beta = $"+ str(b),horizontalalignment='center',
				    verticalalignment='center',
				    transform = ax.transAxes)
            
    name = re.sub('correct.csv','v_4',f)
    plt.savefig(saveDirectory+name+'.png',bboxinches = 'tight') # save image
    plt.close()
    params[0] = params[0]*10**3
    perr[0] = perr[0]*10**3
    fitted_data.loc[i] = [name,L,energy,params[0],perr[0],params[1],perr[1],z0,b,xshift,yshift,w0,w0Err,i0,i0Err,alpha]
    return button
        
if __name__ == '__main__':
    fitted_data = pd.DataFrame(columns = ('Name','L','E(nJ)','Z0(mm)','Zerr','b','berr','Adj Z0','Adj b','delta x','delta y','w0(mm)','werr','I0(Gw/cm^2)','Ierr','alpha'),index=None)
    fileDirectory = '/home/kaspar/Desktop/20151118-correct/'
    f = ReadAllFiles(fileDirectory)    
    for i in range(len(f)):
        button = Main(i,fitted_data)
    fitted_data.to_csv('/home/kaspar/Desktop/2011118.csv',index=False)
        
       
            
    
    
