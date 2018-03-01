#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#import matplotlib
import matplotlib.pyplot as plt
#import matplotlib.mlab as mlab
#import pmagplotlib as pmagplotlib
import numpy as np
import csv
import pylab
import os.path as path
import pmagplotlib as pmagpl
import pyscu_libs as scu
import sys
rad=np.pi/180.
deg=180./np.pi


def main():
    
        
    """
    NAME
        pyscu_draw.py

    DESCRIPTION
        plot the data calculated in pyscu_calc.py
    
    INPUT
        ouput files from pyscu_calc.py

    SYNTAX
        pyscu_draw.py [-h] [command line options]
    
    OPTIONS
        -h, plots help message and quits
        -f AFILE, specify file for input
            this is the main file output from the pyscu_calc.py program
            the other files have to been placed in the same folder
        -F RFILE, specify file for output
        -A, plot the A/N matrix
        -s [s/i/n], plot the SCI solutions, the intersections, or none.
        -fmt [svg, jpg, eps, pdf] format for output images.
        -i interactive entry of the remagnetization direction

    DEFAULTS
    -f, AFILE:  SCdata_main.txt
    -F, RFILE:  out
    DON'T plot the A/N matrix
    -s, DON'T plot the SCI solutions or the intersections
    -fmt, save as svg
    get de remagnetization direction from SCdata_Ref.txt       
    """

    
    print ('\nThis program uses the PmagPy and pySCu softwares utilities\n\tTauxe et al. 2016, G3, http://dx.doi.org/10.1002/2016GC006307\n\tCalvín et al. 2017, C&G, http://dx.doi.org/10.1016/j.cageo.2017.07.002')

    
   
    if '-h' in sys.argv:
        print(main.__doc__)
        sys.exit()

    if '-f' in sys.argv:
        ind = sys.argv.index('-f')
        infile = sys.argv[ind + 1]
    else: infile='SCdata_main.txt'
    
    if '-F' in sys.argv:
        ind=sys.argv.index('-F')
        outfile=sys.argv[ind+1]
    else: outfile=infile[:-9]
    
    if '-A' in sys.argv:
        preA='y'
    else: preA='n'
    
    if '-s' in sys.argv:
        ind=sys.argv.index('-s')
        preS=sys.argv[ind+1]
    else: preS='n'
    
    if '-fmt' in sys.argv:
        ind=sys.argv.index('-fmt')
        fmt='.'+sys.argv[ind+1]
    else: fmt='.svg'
    
    if '-i' in sys.argv:
        iRef='true'
        print('\nInput the Kent parameters (separated by spaces) of the remagnetization direction:')
        ref_input = input("\nDec Inc Eta Dec_Eta Inc_Eta Zeta Dec_Zeta Inc_Zeta: An example...\n329.9 39.5 10.5 155.1 50.4 4.9 62.3 2.6\n").split(' ')
        ref=[]
        for dato in ref_input:
            ref.append(float(dato))
    else: iRef='false'
    
    infile_m=infile[:-8]+'mat.txt'
    infile_ref=infile[:-8]+'Ref.txt'
    infile_inter=infile[:-8]+'inter.txt'
    infile_sci=infile[:-8]+'SCIs.txt'
    
    if path.exists(infile_ref):    Ref='true'
    else: Ref='false'
    
    if path.exists(infile_m):    matrix='true'
    else: matrix='false'
    
    if path.exists(infile_inter):    inter='true'
    else: inter='false'
    
    if path.exists(infile_sci):    SCIs='true'
    else: SCIs='false'
    
    if Ref=='false' and iRef=='false': 
        print("\nTake care, I don't found the file", infile_ref, " whit the reference direction")
    if matrix=='false' and preA=='y':
            print("\nTake care, I don't found the file", infile_m, ' whit the A/N matriz data')
    if inter=='false' and preS=='i':
            print("\nTake care, I don't found the file", infile_inter, ' whit the intersections directions')
    if SCIs=='false' and preS=='s':
            print("\nTake care, I don't found the file", infile_sci, ' whit the SCIs directions')
    
    
    out_name_bbc=outfile+'_bbc'+fmt
    out_name_bfd=outfile+'_bfd'+fmt
    out_name_atbc=outfile+'_atbc'+fmt
    out_name_mat=outfile+'_mat'+fmt
    
    
    print('\nPlease, wait a moment')
    print('\nPlots will be saved as', out_name_bbc, ', ', out_name_bfd, '...\n')
    
    
    #Saving the data in different list
    site,sc,geo,tilt,bfd=scu.getInFile_main(infile) #main file
    n=len(site)
    
    if Ref=='true' and iRef=='false': #reference direction
        reader=csv.reader(open(infile_ref), delimiter=' ')
        dat_Ref=list(reader)
        ref=[float(dat_Ref[1][1]),float(dat_Ref[1][2]),float(dat_Ref[1][3]),float(dat_Ref[1][5]),
        float(dat_Ref[1][6]),float(dat_Ref[1][4]),float(dat_Ref[1][7]),float(dat_Ref[1][8]),float(dat_Ref[1][11])]
    
    if inter=='true' and preS=='i': #intersections directions
        reader=csv.reader(open(infile_inter), delimiter=' ')
        dat_inter_h=list(reader)
        dat_inter=dat_inter_h[1:]
    
    if SCIs=='true' and preS=='s': #intersections directions
        reader=csv.reader(open(infile_sci), delimiter=' ')
        dat_SCIs_h=list(reader)
        dat_SCIs=dat_SCIs_h[1:]
    
    if matrix=='true' and preA=='y': #A/n values
        X,Y,Z,minA,maxA=scu.getInFile_mat(infile_m)  
    
    
    #Drawing...
    plt.figure(num=1,figsize=(6,7),facecolor='white')
    
    #Plotting the BBC directions, the SCs and the reference
    pmagpl.plotNET(1)
    pylab.figtext(.02, .045, 'pySCu v3.1')
    plt.text(0.85, 0.7, 'BBC', fontsize = 13)
    plt.scatter(0.8, 0.74, color='r',marker='s',s=30)
    plt.text(0.70, 0.85, 'n='+str(n), fontsize = 13)
    
    for dato in sc: #The SCs
        scu.smallcirc(dato,1)
    
    for dato in geo: #The BBC directions
        scu.plot_di_mean(dato[0],dato[1],dato[2],color='r',marker='s',markersize=8,label='Geo',legend='no',zorder=3)
        #You can change the marker (+, ., o, *, p, s, x, D, h, ^), the color (b, g, r, c, m, y, k, w) or the size as you prefere
    
    if Ref=='true': #The reference
        scu.plotCONF(ref)
        plt.text(0.51, -1.05, 'Reference', fontsize = 13)
        plt.scatter(0.45, -1, color='m',marker='*',s=100)
    plt.title('Before Bedding Correction',fontsize=15)
    plt.savefig(out_name_bbc)
    
    #Plotting the ATBC directions, the SCs and the reference
    plt.figure(num=2,figsize=(6,7),facecolor='white')
    pmagpl.plotNET(2)
    pylab.figtext(.02, .045, 'pySCu v3.1')
    plt.text(0.85, 0.7, 'ATBC', fontsize = 13)
    plt.scatter(0.8, 0.745, color='g',marker='^',s=40)
    plt.text(0.70, 0.85, 'n='+str(n), fontsize = 13)
    plt.title('After total bedding correction',fontsize=15)
    
    for dato in sc:
        scu.smallcirc(dato,1)
    if Ref=='true':
        scu.plotCONF(ref)
        plt.text(0.51, -1.05, 'Reference', fontsize = 13)
        plt.scatter(0.45, -1, color='m',marker='*',s=100)
    
    for dato in tilt:
        scu.plot_di_mean(dato[0],dato[1],dato[2],color='g',marker='^',markersize=9,label='Tilt',legend='no',zorder=3)
    
    plt.savefig(out_name_atbc)
    
    #Plotting the BFD directions, the SCs and the reference
    plt.figure(num=3,figsize=(6,7),facecolor='white')
    pmagpl.plotNET(3)
    pylab.figtext(.02, .045, 'pySCu v3.1')
    plt.text(0.85, 0.7, 'BFD', fontsize = 13)
    plt.scatter(0.8, 0.74, color='b',marker='o',s=30)
    plt.text(0.70, 0.85, 'n='+str(n), fontsize = 13)
    plt.title('After partial bedding correction',fontsize=15)
    
    for dato in sc:
        scu.smallcirc(dato,1)
    
    for dato in bfd:
        scu.plot_di_mean(dato[0],dato[1],dato[2],color='b',marker='o',markersize=5,label='BFD',legend='no',zorder=3)
    
    if Ref=='true': #Ploting the reference and the leyend
        scu.plotCONF(ref)
        plt.text(0.51, -1.05, 'Reference', fontsize = 13)
        plt.scatter(0.45, -1, color='m',marker='*',s=100)
    plt.savefig(out_name_bfd)
    
    #Plotting the A/n contour plot and/or the intersections
    if (preA=='y' and matrix=='true') or (preS=='i' and inter=='true') or (preS=='s' and SCIs=='true'):
        plt.figure(num=4,figsize=(9.5,9.5),facecolor='white')
        pmagpl.plotNET(4)
        pylab.figtext(.02, .045, 'pySCu v3.1')
        fig4='true'
    else: fig4='false'
            
    if preA=='y' and matrix=='true': #plotting the A/n contour plot
        max_z=max(Z)
        max_z_s=max_z+(5-max_z%5)+0.1
    
        min_z=min(Z)
        min_z_s=min_z-(min_z%5)
    
        levels5 = np.arange(min_z_s,max_z_s, 5)
        levels1 = np.arange(min_z_s,max_z_s, 1)
    
        CS=plt.tricontourf(X, Y, Z, vmin=min_z,vmax=max_z, cmap = 'Blues', levels=levels1) #Other colormaps (as 'rainbow') are possibles. Change 'Blues' for the choosed colormap
        cbar=plt.colorbar(CS, orientation='horizontal',pad=0.05)
        CS2=plt.tricontour(X,Y,Z, colors='k',linewidths = .5, levels=levels5)
    
        #plt.clabel(CS2,levels=levels5, inline=1, fmt='%1.0f', fontsize=10)
        cbar.ax.set_xlabel('A/n value'+' ('+str(round(minA,1))+'-'+str(round(maxA,1))+')')
        #cbar.add_lines(CS2)
        plt.axis((-1.35,1.35,-1.35,1.35))
    else:
        for dato in sc:
            scu.smallcirc(dato,1)
    
    
    if preS=='i' and inter=='true': #plotting the intersections
        text_i='SCs intersec. (n='+str(len(dat_inter))+')'
        plt.text(-0.3, -1.2, text_i, fontsize = 12)
        plt.scatter(-0.38, -1.12, color='k',marker='.',s=50)
        for dato in dat_inter:
            scu.plot_di_mean(float(dato[0]),float(dato[1]),0.,color='k',marker='.',markersize=1,label='Intersections',legend='no')
    
    if preS=='s' and SCIs=='true': #plotting the SCIs
        text_s='SCIs solutions (n='+str(len(dat_SCIs))+')'
        plt.text(-0.3, -1.2, text_s, fontsize = 12)
        plt.scatter(-0.38, -1.12, color='k',marker='.',s=50)
        for dato in dat_SCIs:
            scu.plot_di_mean(float(dato[0]),float(dato[1]),0.,color='k',marker='.',markersize=1,label='SCIs',legend='no')
            
    if fig4=='true' and Ref=='true': #Plotting the reference and the leyend
        scu.plotCONF(ref)
        plt.text(0.6, 0.83, 'Reference', fontsize = 13)
        plt.scatter(0.93, 0.73, color='m',marker='*',s=100)
        text_rat='mr/mp='+str(ref[8])+';'
        plt.text(-1.37, -1.2, text_rat, fontsize = 12)
    plt.title('A/n matriz',fontsize=15)
    
    if fig4=='true':
        plt.savefig(out_name_mat)
    
    
    
    plt.show()
    
main()
