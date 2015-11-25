# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 08:42:33 2015

@author: Kaspar Martin Snashall

script designed to normalise and fit z-scan data using stochiastic method
has final versio  of algorithm v3

license = MIT

You will need python 2.7 and the following packages installed

Pandas
matplotlib
numpy
scipy

The data must be normalised in csv form with columns:

z(mm),open,closed

and name of file must be

x_thickness_energynJ_moreinfo.csv


Use the nomralising program on www.github/KasparSnashall/Z-scan-models
to normalise data from a z-scan.


Data is fitted using scipy's curve fit algorithm that is based on least squares.
To use this program change file directory and save directory to desired paths,
from terminal $ python Open-scan-fitting.py
If the data is in the correct form it should run, then present the user with a fitted plot and gui to change paramters is nessesary.
Once you are happy with the plot exit it and hit enter in the terminal, the program then continues to the next data point.


The program will have to be altered for anyones use however this is just meant to be a rough blue print for data fitting with a GUI.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
from os.path import isfile, join
from scipy import optimize
import sys, re
from matplotlib.widgets import Slider, Button
from vertical import VertSlider



def OpenScanv3(z,z0,beta):
	"""The plot used in curve fitting"""
	X = (z/z0)**2
	m=0 
	mysum = 0 # my sum
	myfactorials = [1,1,2,6,24,120,720,5040,40320,362880,3628800] #factorial numbers up to 10
	while m < 10:
		mysum += ((-beta*leff*Intensity(BeamWaist(z0),energy,pulsetime)*(1/(1+X)))**m)/myfactorials[m]
		m+=1
	return mysum

     
def BeamWaist(Z0):
	"""caluculate the beam waist"""
	wavelength = 1030*10**(-9) # m
	return np.sqrt(abs(Z0)*wavelength/np.pi)

def BeamWaistErr(w0,z0,deltaz0):

    return w0*np.sqrt(deltaz0/z0)

def Intensity(W0,energy,pulsetime):
	"""calculate the intensity at focus"""
	power = energy/pulsetime
	i = power/(np.pi*W0**2)
	return i

def IntensityErr(i0,waistErr):
    return i0*waistErr

def L_eff(L,aplha):
	"""calculate the effective length"""
	return (1-np.exp(-aplha*L))/aplha

def Beta(i0,l,q0):
	"""calculate the non linear absorbtion coefficent"""
	return  q0/(i0*l)
 
def BetaErr(b,q0,q0Err,i0Err,i0):
    return b*np.sqrt((q0Err/q0)**2+(i0Err/i0)**2)
 
def ReadAllFiles(folder):
    """reads all the files in a given directory"""
    return [ f for f in listdir(folder) if isfile(join(folder,f)) and '.csv' in f]
    
def GetParams(f):
    """Get the energy length and alpha from the file name
    Users MUST change this if you want to use this propgram"""
    myl = f.split('_')[1] # the file name
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
        print "unknown sample or error rewrite this code for new samples"
        sys.exit(0)
        
    # clean up the energy so its just a number
    energy1 = f.split('_')[2] 
    energy2 = re.sub('nJ','',energy1)
    #test to see if energy is a number 
    #this is because a mistake made in some file names where 2nd position was not the energy
    try:
        energy = float(energy2)*10**-9
    except:
        energy1 = f.split('_')[3]
        energy2 = re.sub('nJ','',energy1)
        energy = float(energy2)*10**-9
    # test to see if all the variables are numbers
    try:
        print "energy =",float(energy)
        print "length =",float(L)
        print "alpha =",float(alpha)
    except:
        print "energy, length or absorbtion coefficent is not a number"
    return energy,L,alpha



    
def Main(myfile,fitted_data,saveDirectory,fileDirectory):
    """The main function of the program"""
    f =  myfile# ith file name from list of files in directory
    data = pd.read_csv(fileDirectory+f) # grab the data
    z = data.iloc[0:,0]/10**3 # z in meters
    ydata = data.iloc[0:,1] #the transmission    
    
    global energy,L,alpha,leff,pulsetime # define some globals
    pulsetime = 340*10**(-15) # the pulse time
    energy,L,alpha = GetParams(f)
    leff = L_eff(L,alpha)
    
    
    #########################
    # fitting of data
    
    guess= [0.009,0.5] # initial guess
    params, params_covariance = optimize.curve_fit(OpenScanv3,z, ydata, guess) # the fitting
     
    z0 = params[0] # 1st param
    b = params[1] # the NL absorbtion coefficent -> 2nd param
    if z0 < 0:
        z0 = -z0 # make usre it hasnt made z0 negative



    w0 = BeamWaist(z0) # the beam width
    i0 = Intensity(w0,energy,pulsetime) # the intensity at focus        
    
        
    # the one sigma error for the paramters
    perr = np.sqrt(np.diag(params_covariance))        
    w0Err = BeamWaistErr(w0,z0,perr[0])
    i0Err = IntensityErr(i0,w0Err)
    bErr = perr[1]



    #print the initial values
    print params
    print perr
    print w0,w0Err
    print i0,i0Err
    print b,bErr
    
    #create a plot of the transmission and fitted curve
    fig, ax = plt.subplots()
    plt.ylabel('Transmission')
    plt.xlabel('z(m)')
    
    # example SampleX
    samplename =  f.split('_')[0]    
    
    mytitle =  str(samplename)+str(L*10**9) +"nm open scan "+str(energy/10**-9)+"nJ" # title
    plt.title(mytitle) # title 
    myscat = plt.scatter(z,ydata) # data scatter
    plt.ylim( ax.get_ylim() )  
    improvedZ = np.linspace(-0.05,0.05,500) # the x data for the fit
    myfit, = plt.plot(improvedZ,OpenScanv3(improvedZ,params[0],params[1]),color='red') # plot the fitted curve
    plt.subplots_adjust(right = 0.75, bottom=0.25) # alignment
    
    # intial slider values 
    a0 = b 
    f0 = z0
    
    # make sliders for interative plot
    axcolor = 'lightgoldenrodyellow'
    axz0 = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor) # z0 slider
    axb0  = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor) # beta slider
    xvals = plt.axes([0.92, 0.25, 0.03, 0.65], axisbg=axcolor) # x and y adjsutment sliders
    yvals = plt.axes([0.85, 0.25, 0.03, 0.65], axisbg=axcolor)  
    # slider implementation
    xslider = VertSlider(xvals, 'x offset', -0.01, 0.01, valinit=0,valfmt='%.3f')
    yslider = VertSlider(yvals, 'y offset', -0.3, 0.3, valinit=0,valfmt='%.4f')
    myZ0 = Slider(axz0, 'Z0', 0.00001, 0.05, valinit=f0, valfmt='%.5f')
    myB0 = Slider(axb0, 'beta', -0.000005, 0.000005, valinit=a0,valfmt='%.10f')
    
    
    # update function for matplotlibs sliders
    def update(val):
        z0 = myZ0.val
        q0 = myB0.val
        yshift = yslider.val
        xshift = xslider.val
        shiftedy = OpenScanv3(improvedZ,z0,q0)+yshift
        myfit.set_ydata(shiftedy)
        myfit.set_xdata(improvedZ+xshift)
        fig.canvas.draw_idle()
        
    myZ0.on_changed(update)
    myB0.on_changed(update)
    xslider.on_changed(update)
    yslider.on_changed(update)
    
    
    # reset button for sliders
    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
    def reset(event):
        myZ0.reset()
        myB0.reset()
        xslider.reset()
        yslider.reset()
    button.on_clicked(reset)
    
    
    #show the plot
    plt.show(block=False)
    
    # wait for user to finish fitting exit and enter input...
    answer = raw_input('Hit enter to continue\n')
    
    # retrive slider values
    z0 = myZ0.val
    b = myB0.val
    yshift = yslider.val
    xshift = xslider.val
    # close plot to prevent memory overload    
    plt.close()
    
    # ensure z0 is not less then 0 again
    if z0 < 0:
        z0 = -z0
   
    w0 = BeamWaist(z0) # the beam width
    i0 = Intensity(w0,energy,pulsetime) # the intensity at focus        

    
   	
         
    # the one sigma error
    perr = np.sqrt(np.diag(params_covariance))        
    w0Err = BeamWaistErr(w0,z0,perr[0])
    i0Err = IntensityErr(i0,w0Err)
    bErr = perr[1]

    # print the final versions
    print "Z0 = ",z0,"beta = ",b
    print "err0r = ",perr
    print "w0 = ",w0,w0Err
    print "i0 = ",i0,i0Err
    print "b = ", b,bErr
    
    #make the final plot
    fig, ax = plt.subplots()
    plt.ylabel('Transmission')
    plt.xlabel('z(m)')
    plt.title(mytitle) # title 
    plt.scatter(z,ydata)
    improvedZ = np.linspace(-0.05,0.05,500)
    myfit, = plt.plot(improvedZ+xshift,OpenScanv3(improvedZ,z0,b)+yshift,color='red')
    # format some of the data into the standard form
    z0 = z0*10**3
    w0 = w0*10**6
    i0 = i0/(10**(9)*100*100)
    energy = energy/10**-9
    L = L*10**9
    
    # add some text showing the values
    plt.text(0.8,0.8,"$Z_0 (mm)$ ="+str('%.2f' % z0)+"\n"+'$w_0(\mu m)$ = '+str('%.2f' % w0),horizontalalignment='center',
				    verticalalignment='center',
				    transform = ax.transAxes)
    plt.text(0.2,0.8,"$I_0 (GW cm^{-2}) = $"+str('%.2f' % i0)+"\n"+r"$\beta = $"+ str(b),horizontalalignment='center',
				    verticalalignment='center',
				    transform = ax.transAxes)
    
    # make a sname for the entry to be put into a csv
    name = str(re.sub('.csv','',f)) # remove .csv from the title
    plt.savefig(saveDirectory+name+'.png',bboxinches = 'tight') # save image
    plt.close()
    
    print name,L,energy,params[0],perr[0],params[1],perr[1],z0,b,xshift,yshift,w0,w0Err,i0,i0Err,alpha
    
    
    params[0] = params[0]*10**3 # ensure z0 is in meters for the csv
    perr[0] = perr[0]*10**3
    
    
    
    # add an entry to the csv
    fitted_data.loc[i] = [name,L,energy,params[0],perr[0],params[1],perr[1],z0,b,xshift,yshift,w0,w0Err,i0,i0Err,alpha]
    return button # return button this ensures the plot works (if removed the reset button will not work)
        
if __name__ == '__main__':
    # a blank dataframe to put data into    
    fitted_data = pd.DataFrame(columns = ('Name','L','E(nJ)','Z0(mm)','Zerr','b','berr','Adj Z0','Adj b','delta x','delta y','w0(mm)','werr','I0(Gw/cm^2)','Ierr','alpha'),index=None)
    
    fileDirectory = '/home/kaspar/Desktop/test/' # your file directory
    saveDirectory = '/home/kaspar/Desktop/test/' # your save directory
    f = ReadAllFiles(fileDirectory) # read all the fils in 
    for i in range(len(f)):
        # iterate over the files and fit
        button = Main(f[i],fitted_data,fileDirectory,saveDirectory)
    #once the fitting is complete create a csv        
    print fitted_data
    fitted_data.to_csv('/home/kaspar/Desktop/test.csv',index=False)
        
       
            
    
    
