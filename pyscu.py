#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import sys
import pyscu_libs as scu
rad=np.pi/180
deg=180/np.pi
from time import time

def main():
    
    
    """
    NAME
        pyscu_calc.py

    DESCRIPTION
        perform the routine explained in Calvín et al. (2017) 
                  https://doi.org/10.1016/j.cageo.2017.07.002
        calculates the Small Circle Intersection from a dataset
        calculates the paleodip of the beds
    
    INPUT
        spaced separated text file with header, with the next columns
        site dec inc kappa dipDir dip kappa

    SYNTAX
        pyscu_calc.py [-h] [command line options]
    
    OPTIONS
        -h plots help message and quits
        -f AFILE, specify file for input
        -F RFILE, specify file for output
        -i DON'T calculate the remagentization direction. Interactive entry of that
        -A calculate the A/N matrix
        -neg negative inclination of the remagnetization dir. is expected
        -nb N; specifies the number of bootstraps - default is 500

    DEFAULTS
        AFILE:  SCdata.txt
        RFILE:  SCdata_results.txt, SCdata_inter.txt, SCdata_mat.txt, SCdata_SCIs.txt
        Looking for a positive paleomagnetic direction
        Don't calculate the A/n matrix
        500 bootstraps
        """

    if '-h' in sys.argv:
        print(main.__doc__)
        sys.exit()

    if '-f' in sys.argv:
        ind = sys.argv.index('-f')
        files = sys.argv[ind + 1]
    else: files='SCdata.txt'
    
    if '-F' in sys.argv:
        ind=sys.argv.index('-F')
        out_file=sys.argv[ind+1]
    else: out_file=files[:-4]
    
    if '-nb' in sys.argv:
        ind=sys.argv.index('-nb')
        nb=int(sys.argv[ind+1])
    else: nb=500
    
    if '-neg' in sys.argv:
        pregunta_inc='neg'
    else: pregunta_inc='pos'
    
    if '-A' in sys.argv:
        pregunta_matriz='y'
    else: pregunta_matriz='n'


    data,site,geo,bed_d,bed_s,bed_pole,N_sites=scu.saveInputFile(files)   #Saving in diferent list the data
    tilt=scu.tilt_rot(geo,bed_s)                   #Calculating TILT directions


    time_ini=time()

    #Calculating BFD and A/N for a manual input remagnetization direction.
    if '-i' in sys.argv:
        pregunta_remag='n'
        DecInp = float(input("Declination of the remagnetization direction: "))
        IncInp = float(input("Inclination of the remagnetization direction: "))
        ref=[DecInp, IncInp]
        BFD,A_min=scu.calAQQ2(geo,bed_s,ref)
        Qmean_min=scu.fisher_mean(BFD)
        N=len(BFD)
        An=A_min/N
        Dir_remag=[N,ref[0],ref[1],round(A_min,3),round(An,3)]
        print('Used direction: ', 'Dec / Inc ', '(',"%.1f" % Dir_remag[1],'/',"%.1f" % Dir_remag[2],')', 'A/n ',"%.3f" % An)
    else: #Calculating SCI_solutions, their mean (the remagnetization direction) and the rest of parameters
        pregunta_remag='y'
        SCIs=[]
        pgeo=scu.para_dir(geo)
        pbed_pole=scu.para_dir(bed_pole)
        w,f,cont=0,0,nb/10
        print('\nplease, be patient... calculating',nb,'SCI_solutions')
        for l in range (nb): 
            w+=1
            dec=np.random.randint(0,359,1)
            inc=np.random.randint(1,89,1)
            if pregunta_inc=='neg': inc=inc*(-1)
            point=[dec[0],inc[0]]
            pgeo_u=scu.selec_para_geo(pgeo)
            pbed_u=scu.selec_para_pole(pbed_pole)
            pQ=scu.p_calAQQ2(pgeo_u, pbed_u, point)
            Qmean_min=scu.p_minA(pgeo_u, pbed_u,pQ)
            SCIs.append(Qmean_min)
            if w==cont:
                f+=1
                if f!=8:
                    print("%.0f"%cont)
                    cont+=nb/10.
                else: 
                    print('almost finished...')
                    cont+=nb/10.
            
        Qend=scu.dokent(SCIs,1.)
        ref=[Qend['dec'], Qend['inc']]
        BFD,A_min=[],[]
        BFD,A_min=scu.calAQQ2(geo,bed_s,ref)
        An=A_min/N_sites
        




    distance=[] #This is the angle between the BFD and the remagnetization direction, for each site
    for dato in BFD:
        distance_site=scu.ang2point(dato,ref)
        distance.append(distance_site)


    paleobed=scu.paleo_dip(tilt,bed_s,BFD)


    api=scu.cal_api(geo,bed_s)
    api2=scu.cal_api2(geo,bed_s)

    out_inter=scu.inter(api2,site)
    
    mr=len(out_inter)
    mp=N_sites*(N_sites-1.)/2.
    mr_mp=round(mr/mp,2)
    
    

    #joining the data in a unique list
    out_main=[]
    
    if pregunta_remag=='y':
        Dir_remag=[N_sites,round(Qend['dec'],1),round(Qend['inc'],1),round(Qend['Eta'],1),round(Qend['Zeta'],1),
        round(Qend['Edec'],1),round(Qend['Einc'],1),round(Qend['Zdec'],1),round(Qend['Zinc'],1),round(An,3),round(A_min,3),mr_mp,nb]
        print('\nKent mean remagnetization direction (Kent, 1982; Tauxe et al., 1991)', '\nDec / Inc: ', "%.1f" % Dir_remag[1],"/","%.1f" % Dir_remag[2])
        print('A/n: ',"%.3f" % An,'mr/mp: ',"%.2f" % mr_mp)
        print('Eta_95, dec, inc:', "%.1f" % Dir_remag[3], "%.1f" % Dir_remag[5], "%.1f" % Dir_remag[6])
        print('Zeta_95, dec, inc', "%.1f" % Dir_remag[4], "%.1f" % Dir_remag[7], "%.1f" % Dir_remag[8])
    
    for i in range(len(site)):
        site_main=[site[i][0],data[i][1],data[i][2],data[i][3],data[i][4],data[i][5],data[i][6],data[i][7],"%.1f" %tilt[i][0],"%.1f" %tilt[i][1],"%.1f" %BFD[i][0],"%.1f" %BFD[i][1],"%.0f" %bed_s[i][0],"%.2f" %api[i],"%.0f" %((paleobed[i][0]+90)%360),"%.0f" %paleobed[i][1],"%.1f" %distance[i]]
        out_main.append(site_main) 
        

    #saving the files
    header_main=['Site','Dec_BBC','Inc_BBC','a95', 'k','DipDir', 'Dip', 'k_bed','Dec_ATBC','Inc_ATBC','Dec_BFD','Inc_BFD','SC_Strike','SC_Api_angle','Paleo_DipDir','Paleo_dip','Ref_BFD_angle']
    name_main=out_file+'_main'
    scu.save_out_file(header_main,out_main,name_main)

    header_inter=['Dec','Inc','Site_i','Site_j']
    name_inter=out_file+'_inter'
    scu.save_out_file(header_inter,out_inter,name_inter)

    if pregunta_remag=='y':
        header_Ref=['N','Dec','Inc','Eta','Zeta','Dec_E','Inc_E','Dec_Z','Inc_Z','A/n','Asum','mr_mp','number of SCI_solutions']
        name_Ref=out_file+'_Ref'
        scu.save_out_file(header_Ref,[Dir_remag],name_Ref)
        header_SCIs=['Dec','Inc']
        name_SCIs=out_file+'_SCIs'
        scu.save_out_file(header_SCIs,SCIs,name_SCIs)

    if pregunta_matriz=='y':
        print('\n\nTo be patient! To calculate the matrix spend some min')
        mat=scu.cal_matriz(geo,bed_s)
        header_matriz=['Dec','Inc','x_eqarea','y_eqarea','x_eq_normalized','y_eq_normalized','A','A/n']
        name_matriz=out_file+'_mat'
        scu.save_out_file(header_matriz,mat,name_matriz)
        

    time_fin=time()
    time_tot=time_fin-time_ini
    print('\nExecution time: ',round(time_tot,1),'seg')
    
    print('\nYou can use pyscu_draw.py and pyscu_draw_labels.py to print your results')
    print('Good bye')
    

main()

