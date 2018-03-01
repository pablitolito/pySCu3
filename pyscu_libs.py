import pmag as pmag
import csv
import numpy as np
import sys
rad=np.pi/180
deg=180/np.pi
import matplotlib.pyplot as plt
import os.path as path
from past.utils import old_div

def selec_para_geo(pdir):
	'''
	Get one data from a input list (plist) which is a list of list with 30 objets per site of bedding or pmag. dir
	Out a list of list with Dec, Inc and alpha95 or DipDir, Strike, dip and bedding_error, depending the input data
	Generically, get randomly one of the 30 sublist which forms the N_sites lists
	'''
	
	plist_u=[]
	#index=np.random.randint(0,99,size=N_sites)
	#i=0
	
	for dato in pdir:
		i=np.random.randint(0,50)
		plist_u.append(dato[i])
		#i+=1
	return plist_u
	
def selec_para_pole(ppoles):
	'''
	Get one data from a input list (plist) which is a list of list with 30 pairs of dec and inc of the pole bedding and converts it to strike and dip
	Out a list of list with Strike and dip
	Generically, get randomly one of the 30 sublist which forms the N_sites lists
	'''
	
	plist_u=[]
	#index=np.random.randint(0,99,size=N_sites)
	#i=0
	plist_u_site=[]
	
	for dato in ppoles:
		#print dato[index[i]]
		i=np.random.randint(0,50)
		strike=(dato[i][0]+90)%360
		dip=90-dato[i][1]
		plist_u_site=[strike,dip]
		plist_u.append(plist_u_site)
		#i+=1
		
	return plist_u
	
def para_dir(list):
	'''
	Given a list with dec, inc, and kappa, calculates throught pmag.fishrot subrutine, 
	calculates para-datasets a list of list with 50 pairs of Dec, Inc of directions with the same fisherian parameters than input data
	'''
	
	plist=[]
	plist_site=[]
	for dato in list:
		plist_site=fishrot(dato[2],50,dato[0],dato[1])
		plist.append(plist_site)
	return plist

def fishrot(kappa,N,D,I): #from Pmagpy
    """
    Description: generates set of Fisher distributed data from specified distribution 
	Input: kappa (fisher distribution concentration parameter), number of desired subsamples, Dec and Inc
	Output: list with N pairs of Dec, Inc.
    """
    out_d=[]
    out=[]

    for k in range(N): 
        dec,inc= pmag.fshdev(kappa)  # send kappa to fshdev
        drot,irot=pmag.dodirot(dec,inc,D,I)
        out_d=[drot,irot]
        out.append(out_d)
    return out
	
def getk(alpha95,N):
	'''
	Calculater K value for a input data alpha95 and N
	'''
	
	P=((1./0.05)**(1./(N-1.)))-1.
	k=((N-1.)/N)*(1.+(P/(1.-np.cos(rad*alpha95))))
	return k

def getkaprox(alpha95,N):
	'''
	calculates K value using an approximate expresion
	'''
	k=(140.**2.)/(((alpha95)**2.)*N)
	
	return k

def saveInputFile(files):
	'''
	Get data from a txt file separated by spaces
	Input: string with the name of the file
		   file: (Site, Dec, Inc, alfa95, dip direction, dip, bedding error) In in situ coordinates
	Output: (i) data: list of list with Site, Dec, Inc, alpha95, dipdir, strike, dip (without header)
			(ii) geo: list of list with values of magnetization in situ; Dec, Inc, alpha95 (without header)
			(iii) bed: list of list with bedding data; dip direction, Strike (RHR), Dip, error
	
	'''
	
	#beig 'files' a text file separated with spaces, with 6 columns and header (Site, Dec, Inc, alfa95, , dip direction, dip, bedding error)
	#Return a list with the same but without header, and with the numbers as float
	#Saving the input file/data 
	#output have the same colums, but whithou header.
	#file = open(files)
	reader=csv.reader(open(files, 'rU'), delimiter=' ')
	dat=list(reader)
	#file.close()
	data=dat[1:]    #guardamos los datos sin el encabezado
	data_float=data[:]

	#Converting the input data to float
	for i, lista in enumerate(data):
		for j, val in enumerate(lista):
			if j>0:   #condicional es para dejar el nombre del site tranquilo
				data_float[i][j]=float(val)
	
	geo=[]
	bed_d=[]
	bed_s=[]
	bed_pole=[]
	site=[]
	#print data_float

	for dato in data_float:
		site_site=[dato[0]]
		site.append(site_site)
		geo_site=[dato[1],dato[2], dato[4]]
		geo.append(geo_site)
		bed_d_site=[dato[5],dato[6],dato[7]]
		bed_d.append(bed_d_site)
		bed_s_site=[(dato[5]-90)%360,dato[6],dato[7]]
		bed_s.append(bed_s_site)
		bed_p_site=[(dato[5]+180)%360,90-dato[6], dato[7]]
		bed_pole.append(bed_p_site)

	N=len(geo)
	
	return data_float, site, geo, bed_d, bed_s, bed_pole, N

def dir2car(dir): #being 'dir' a list of list, with pairs of Dec,Inc in degrees
	#changes from geographic to cartesians coord
	cart=[]   #output is a list of list, with x,y,z
	rad=np.pi/180

	for site in dir:
		x=np.cos(rad*site[0])*np.cos(rad*site[1]) 
		y=np.sin(rad*site[0])*np.cos(rad*site[1]) 
		z=np.sin(rad*site[1])                 
		cartSite=[x,y,z]
		cart.append(cartSite)
	return cart

def car2dir(cart): #being 'cart' a list of list, with x, y, z
	##changes from geographic to cartesians coord
	dir=[]    #output is a list of list, with Dec and Inc (in degree)
	deg=180/np.pi
	
	for site in cart:
		Dec=round(deg*np.arctan2(site[1],site[0]),1)
		Inc=round(deg*np.arcsin(site[2]), 1)
		dirSite=[Dec,Inc]
		dir.append(dirSite)
	return dir

def fisher_mean(dir):
	'''
	It does fisher mean
		Input: list of list, with pairs of Dec, Inc (degrees)
		Output: a list with [Dec_mean, Inc_mean, kappa, alfa95]
	'''
	
	fisher_mean=[]
	
	#change to cartesian coordinats using dir2car defined funcion
	car=dir2car(dir)
	N=len(dir)
	X,Y,Z=0,0,0
	
	for site in car:
		X+=site[0]
		Y+=site[1]
		Z+=site[2]

	R=np.sqrt(X**2+Y**2+Z**2)
	k=(N-1.)/(N-R)
	
	cosAlfa=1.-((N-R)/R)*(20.**(1./(N-1.))-1.)
	if cosAlfa<-1:cosAlfa=-1
	alpha95=deg*np.arccos(cosAlfa)

	#calculating mean direction
	DecMean=deg*np.arctan(Y/X)
	if X<0.:
		DecMean=DecMean+180.
	elif Y<0.:
		DecMean=DecMean+360.
	IncMean=deg*np.arcsin(Z/R)
	fisher_mean=[DecMean,IncMean,k,alpha95]
	return fisher_mean

def tilt_rot(m,b):
	"""
	This function allow to apply the bedding correction to a vector
	Input data: list of list with [m] the magnetic vector ( Dec and Inc in degrees) 
								  [b] the corresponding bed strike (RHR) and the dip of the bed
								  len(m)=len(b)
	Ouput data: list of list with the rotated vectors
	"""
		
	rotated=[]
	
	
	for datom, datob in zip(m,b):
		D=rad*datom[0]   #declination
		I=rad*datom[1]    #inclination
		s=rad*datob[0]    #strike
		d=rad*datob[1]    #dip
		dif=D-s
				  
		x=np.sin(s)*np.cos(dif)*np.cos(I)+np.cos(s)*np.cos(d)*np.sin(dif)*np.cos(I)+np.cos(s)*np.sin(d)*np.sin(I)
		y=np.cos(s)*np.cos(dif)*np.cos(I)-np.sin(s)*np.cos(d)*np.sin(dif)*np.cos(I)-np.sin(s)*np.sin(d)*np.sin(I)
		z=-np.sin(d)*np.sin(dif)*np.cos(I)+np.cos(d)*np.sin(I)
		
		Dec=deg*np.arctan(x/y)
		Inc=deg*np.arcsin(z)
		
		if y<0:
			Dec=Dec+180.
		elif x<0:
			Dec=Dec+360.
		
		outSite=[Dec,Inc]
		rotated.append(outSite)

	return rotated

def paleo_dip(tilt,bed,q):
	"""
	This function allow to calculate the paledip
	For this, it rotates the TITL magnetization and the BFD with a axis perpendicular to the strike.
	In this way, the strike will be upright, the unfolding angle is the difference between TILT and BFD and the sense of rotation is clear
	
	Input data: list of list with [tilt] Dec and Inc (in degrees) of the tilt magnetic vector
								  [bed] strike (RHR) and dip of bedding
								  [q] Dec and Inc (in degrees) of the BFD for each site
	Ouput data: list of list with the paleobuz (Strike and dip)
	"""
		
	paleobed=[]    
	rot=[]
	
	#calculating the rotation axis that allow rotate the strike to the upright position (looking down)
	for dato in bed:
		s=(dato[0]+90.)%360.   #declination
		sit=[s,90.]
		rot.append(sit)
	
	#rotating both tilt and q data.
	tilt_r=tilt_rot(tilt,rot)
	q_r=tilt_rot(q,rot)
	
	#calculation BFD for each site
	for dato_tilt, dato_q, dato_bed in zip(tilt_r,q_r,bed):
		paleo_dip=(dato_q[0]-dato_tilt[0])%360.
		if paleo_dip>=180:
			paleo_st=(dato_bed[0]-180.)%360.
			paleo_dip=360.-paleo_dip
		else:
			paleo_st=dato_bed[0]
		
		sit=[paleo_st,paleo_dip]
		paleobed.append(sit)  
		
	return paleobed

def calAQQ2(geo,bed,point):
	"""
	This function calculates (i)the directions over each small circle closest to a one point: Q
							 (ii)the sum of all minimum angles between one point (P) and all SCs: A
	Input data: is a list of list (without headets) with in situ Dec, Inc and alpha95: geo
				list of list with bedding data (strike -RHR-, dip and error): bed
				list with dec, inc of the reference direction: point
	Output: (i) Ai, the sum of the angles.
			(ii)out, a list of list with the closest directions: Q
	"""
		
	q_rot=[]
	rot=[]
	unrot=[]
	point_l=[]
	A=0
	
	for dato in bed:
		s=(rad*dato[0]-rad*90.)%(2.*np.pi)   #rotation parameters
		sit=[deg*s,90.]
		unsit=[deg*s,-90.]
		rot.append(sit)
		unrot.append(unsit)
		point_l.append(point)
	
	geo_rot=tilt_rot(geo,rot)
	p_rot=tilt_rot(point_l,rot)
	
	for dato_geo, dato_p in zip(geo_rot,p_rot):       
		Dqi=dato_p[0]
		Iqi=dato_geo[1]
		
		alfai=abs(dato_geo[1]-dato_p[1])      
		A += alfai
		
		site=[Dqi,Iqi]
		q_rot.append(site)
		
	q=tilt_rot(q_rot,unrot)
		
	return q,A

def calA(geo,bed,point):
	"""
	This function calculates (i)the sum of all minimum angles between one point (P) and all SCs
	Input: geo: list of list with Dec and Inc (in degrees)
		   bed: list of list with Strike and dip
		   point: list with the Dec and Inc of one direction
	Output: (i) A, the sum of the angles.

	"""
		
	q_rot=[]
	rot=[]
	point_l=[]
	A=0
	
	for dato in bed:
		s=(rad*dato[0]-rad*90.)%(2.*np.pi)   #rotation parameters
		sit=[deg*s,90.]
		rot.append(sit)
		point_l.append(point)
	
	geo_rot=tilt_rot(geo,rot)
	p_rot=tilt_rot(point_l,rot)
	
	for dato_geo, dato_p in zip(geo_rot,p_rot):       
		alfai=abs(dato_geo[1]-dato_p[1])
		A += alfai
		
	return A

def cal_matriz(geo,bed):
	
	matriz_pos=[]
	d_mat=0.
	n=len(geo)
	i=0
	j=0
	count=3240
	
	while d_mat<360:
		i_mat=0.000001
		while i_mat<90:
			point=[d_mat,i_mat]
			x=2.*np.sin(rad*(90.-np.absolute(i_mat))/2.)*np.sin(rad*d_mat)
			xx=x/np.sqrt(2)
			y=2.*np.sin(rad*(90.-np.absolute(i_mat))/2.)*np.cos(rad*d_mat)
			yy=y/np.sqrt(2)
			A=round(calA(geo,bed,point),3)
			An=round(A/n,3)
			out=[d_mat,round(i_mat,0),x,y,xx,yy,A,An]
			matriz_pos.append(out)
			i_mat+=1
			i+=1
			if i==count:
				j+=10
				count+=3240
				print(j,'%')
		d_mat+=1
	d_mat=0
	
	return matriz_pos

def save_out_file(header,data,name):
	name_csv=name+'.txt'
	data_out=[header]
	data_out.extend(data)
	file_out=open(name_csv, 'w')
	writer=csv.writer(file_out, delimiter=' ', lineterminator='\n')
	for row in data_out:
		writer.writerow(row)
	del writer
	file_out.close()

def ang2point(a,b):
	"""
	It calculates the angle between two points along a great circle.
	The coordinates of the points are introduced as two list; Dec and Inc (degrees).
	Te output is a two decimal float with the angle in degrees
	"""

	xa=np.cos(rad*a[0])*np.cos(rad*a[1]) 
	ya=np.sin(rad*a[0])*np.cos(rad*a[1]) 
	za=np.sin(rad*a[1])                 
	a=[xa,ya,za]
	
	xb=np.cos(rad*b[0])*np.cos(rad*b[1]) 
	yb=np.sin(rad*b[0])*np.cos(rad*b[1]) 
	zb=np.sin(rad*b[1])
	b=[xb,yb,zb]
	c=deg*np.arccos(np.dot(a,b))
	c=round(c,2)
	return c

def cal_api(geo,bed):
	"""
	The apical angle of the small circle is calculated for each site.
	input: 
		geo: list of list with Dec and Inc (in degrees)
		strike: list of list with strike an dip (only strike is used)
	"""
	
	api=[]
	
	for dato_geo, dato_bed in zip(geo,bed):
		api_site=deg*np.arccos(np.cos(rad*(dato_geo[0]-dato_bed[0]))*np.cos(rad*dato_geo[1]))
		api.append(api_site)
		
	return api

def cal_api2(geo,bed):
	"""
	The apical angle of the small circle is calculated for each site.
	input: 
		geo: list of list with Dec and Inc (in degrees)
		strike: list of list with strike an dip (only strike is used)
	output:
		At: a list with the apical angle and the strike (in degrees)
	"""
	
	At=[]
	for dato_geo, dato_bed in zip(geo,bed):
		api_site=deg*np.arccos(np.cos(rad*(dato_geo[0]-dato_bed[0]))*np.cos(rad*dato_geo[1]))
		trend=dato_bed[0]
		at_site=[api_site,trend]
		At.append(at_site)
	return At
	
def intersec(sc1,sc2):
	"""
	This function calculates the posible intersection between two give SCs
	input:
		SC1, SC2, a couple of SC with [apical angle, trend]
	output:
		[Dec, Inc]
	"""
	d1=np.cos(rad*sc1[0])
	d2=np.cos(rad*sc2[0])
	rt1=rad*sc1[1]
	rt2=rad*sc2[1]
	a=d1*np.cos(rt2)-d2*np.cos(rt1)
	b=d1*np.sin(rt2)-d2*np.sin(rt1)
	D=np.arctan(-a/b)
	I=deg*np.arccos(d1/np.cos(D-rt1))
	D=(deg*D+360.)%360.
	out=[round(D,2),round(I,2)]
	
	return out

def inter(api2,site):
	"""
	This function calculates the posible intersection between the SCs
	First, it test if two SCs intersect, and if true, intersec(SC1,SC2) calculates de intersection
	input:
		api2: a list of list with the apical angle and the trend (in degrees) of each site
		
	"""
	i=0
	j=1
	l=len(api2)
	
	in_out=[]

	while i<l:
		while j<l:
			difSt=deg*np.arccos(np.cos(rad*(api2[i][1]-api2[j][1]))) #calculando la diferencia de strike
			sumAp=api2[i][0]%90+api2[j][0]%90
			resAP=np.abs(api2[i][0]-api2[j][0])
			if sumAp>=difSt:
				if resAP<=difSt:
					inter_ij=intersec(api2[i],api2[j])
					inter_ij.extend(site[i])
					inter_ij.extend(site[j])
					in_out.append(inter_ij)
			j+=1
		i+=1
		j=i+1
	return in_out
		
def minA(geo,bed,Q,A):
	'''
	This function found the remagnetization direction following the iterative method
	'''
	angle=2. #inicialiting angle
	Qmean=fisher_mean(Q)
	while angle>0.05:
		Qnew,A=calAQQ2(geo,bed,Qmean)
		Qnew_mean=fisher_mean(Qnew)
		angle=ang2point(Qnew_mean,Qmean)
		Qmean=Qnew_mean[:]
		
	return Qnew_mean,Qnew,A
		   
def p_minA(geo,bed,Q):
	'''
	This function found the pseudo_remagnetization direction following the iterative method
	'''
	angle=2. #inicialiting angle
	Qmean=fisher_mean(Q)
	while angle>0.01:
		Qnew=p_calAQQ2(geo,bed,Qmean)
		Qnew_mean=fisher_mean(Qnew)
		#print Qnew_mean
		#print Qmean
		angle=ang2point(Qnew_mean,Qmean)
		Qmean=Qnew_mean[:2]
		
	return Qmean

def p_calAQQ2(geo,bed,point):
	"""
	This function calculates (i)the directions over each small circle closest to a one point: Q
							 (ii)the sum of all minimum angles between one point (P) and all SCs: A
	Input data: is a list of list (without headets) with in situ Dec, Inc and alpha95: geo
				list of list with bedding data (strike -RHR-, dip and error): bed
				list with dec, inc of the reference direction: point
	Output: (i) Ai, the sum of the angles.
			(ii)out, a list of list with the closest directions: Q
	"""
	
	q_rot=[]
	rot=[]
	unrot=[]
	point_l=[]
	
	for dato in bed:
		s=(rad*dato[0]-rad*90.)%(2.*np.pi)   #rotation parameters
		sit=[deg*s,90.]
		unsit=[deg*s,-90.]
		rot.append(sit)
		unrot.append(unsit)
		point_l.append(point)
	
	geo_rot=tilt_rot(geo,rot)
	p_rot=tilt_rot(point_l,rot)
	
	for dato_geo, dato_p in zip(geo_rot,p_rot):       
		Dqi=dato_p[0]
		Iqi=dato_geo[1]       
		site=[Dqi,Iqi]
		q_rot.append(site)
		
	q=tilt_rot(q_rot,unrot)
		
	return q

def getInFile_main(files):
	'''
	Save a *txt file (files) in a list of list without its header.
	Input: string with the name of the file
			This file is the *_main.txt file, which is the main output of pySCu_calc.py module	
	'''
	
	
	#Saving the input file/data 
	reader=csv.reader(open(files, 'rU'), delimiter=' ')
	dat=list(reader)
	data=dat[1:]    #removing the header
	data_float=data[:]

	#Converting the input data to float
	for i, lista in enumerate(data):
		for j, val in enumerate(lista):
			if j>0:   #condicional es para dejar el nombre del site tranquilo
				data_float[i][j]=float(val)
				
	geo=[]
	tilt=[]
	site=[]
	sc=[]
	bfd=[]

	for dato in data_float:
		site_site=[dato[0]]
		site.append(site_site)
		dif=(dato[12]-dato[1])%360
		if dif>90 and dif<270: 
			st=(dato[12]+180)%360
		else: st=dato[12]
		sc_site=[st,0.,dato[13]]
		sc.append(sc_site)
		geo_site=[dato[1],dato[2], dato[3]]
		geo.append(geo_site)
		tilt_site=[dato[8],dato[9], dato[3]]
		tilt.append(tilt_site)
		bfd_site=[dato[10],dato[11], dato[3]]
		bfd.append(bfd_site)
	
	return site,sc,geo,tilt,bfd
 
def getInFile_mat(files):
	'''
	Save a *txt file (files) in a list of list without its header.
	Input: string with the name of the file
		   file: (Site, Dec, Inc, alfa95, dip direction, strike, dip)
	Output: (i) data: list of list with Site, Dec, Inc, alfa95, strike, dip (without header)
			(ii) geo: list of list with values of magnetization GEO; Dec, Inc, alfa95 (without header)
			(iii) bed: list of list with bedding data; Strike (RHR) and Dip
	
	'''
	
	#beig 'files' a csv comma separated with 5 columns and header (Site, Dec, Inc, alfa95, bed strike)   #Return a list with the same but without header, and with the numbers as float
	#Saving the input file/data 
	#output have the same colums, but whithou header.
	#file = open(files)
	reader=csv.reader(open(files, 'rU'), delimiter=' ')
	dat=list(reader)
	data=dat[1:]    #removing the header
	data_float=data[:]

	#Converting the input data to float
	for i, lista in enumerate(data):
		for j, val in enumerate(lista):
			data_float[i][j]=float(val)
	
	X=[]
	Y=[]
	Z=[]
	
	for dato in data_float:
		x=dato[4]
		X.append(x)
		y=dato[5]
		Y.append(y)
		z=dato[7]
		Z.append(z)
	minA=min(Z)
	maxA=max(Z)
	return X,Y,Z,minA,maxA

def smallcirc(a,zorder=1):#Modified from PmagPy (Tauxe et al., 2016)
	Ds,Is=pmag.circ(a[0],a[1],a[2])
	Xcirc,Ycirc=[],[]
	for k in range(len(Ds)):
		XY=pmag.dimap(Ds[k],Is[k])
		Xcirc.append(XY[0])
		Ycirc.append(XY[1])
	plt.plot(Xcirc,Ycirc,'0.75',linewidth = 0.5,zorder=zorder)

def plot_di_mean(dec,inc,a95,color='k',marker='o',markersize=20,label='',legend='no',zorder=3):#Modified from PmagPy (Tauxe et al., 2016)
	"""
	Plot a mean direction (declination, inclination) with alpha_95 ellipse on
	an equal area plot.

	Before this function is called, a plot needs to be initialized with code
	that looks something like:
	>fignum = 1
	>plt.figure(num=fignum,figsize=(10,10),dpi=160)
	>ipmag.plot_net(fignum)

	Required Arguments
	-----------
	dec : declination of mean being plotted
	inc : inclination of mean being plotted
	a95 : a95 confidence ellipse of mean being plotted

	Optional Keywords
	-----------
	color : the default color is black. Other colors can be chosen (e.g. 'r').
	marker : the default is a circle. Other symbols can be chosen (e.g. 's').
	markersize : the default is 20. Other sizes can be chosen.
	label : the default is no label. Labels can be assigned.
	legend : the default is no legend ('no'). Putting 'yes' will plot a legend.
	"""
	DI_dimap=pmag.dimap(dec,inc)
	if inc < 0:
		plt.scatter(DI_dimap[0],DI_dimap[1],
		edgecolors=color ,facecolors='white',
		marker=marker,s=markersize,label=label,zorder=4)
	if inc >= 0:
		plt.scatter(DI_dimap[0],DI_dimap[1],
		edgecolors=color,facecolors=color,
		marker=marker,s=markersize,label=label,zorder=zorder)
	Xcirc,Ycirc=[],[]
	Da95,Ia95=pmag.circ(dec,inc,a95)
	if legend=='yes':
		plt.legend(loc=2)
	for k in range(len(Da95)):
		XY=pmag.dimap(Da95[k],Ia95[k])
		Xcirc.append(XY[0])
		Ycirc.append(XY[1])
	plt.plot(Xcirc,Ycirc,c=color,linewidth=0.5,zorder=3)
	plt.tight_layout()

def plotELL(pars,col,lower,plot): #Modified from PmagPy (Tauxe et al., 2016)
    """
    function to calculate points on an ellipse about Pdec,Pdip with angle beta,gamma
    """
    #pylab.figure(num=fignum)
    Pdec,Pinc,beta,Bdec,Binc,gamma,Gdec,Ginc=pars[0],pars[1],pars[2],pars[3],pars[4],pars[5],pars[6],pars[7]
    if beta > 90. or gamma>90:
        beta=180.-beta
        gamma=180.-beta
        Pdec=Pdec-180.
        Pinc=-Pinc
    beta,gamma=beta*rad,gamma*rad # convert to radians
    X_ell,Y_ell,X_up,Y_up,PTS=[],[],[],[],[]
    nums=201
    xnum=float(nums-1.)/2.
# set up t matrix
    t=[[0,0,0],[0,0,0],[0,0,0]]
    X=pmag.dir2cart((Pdec,Pinc,1.0)) # convert to cartesian coordintes
    if lower==1 and X[2]<0:
       for i in range(3):
           X[i]=-X[i]
# set up rotation matrix t
    t[0][2]=X[0]
    t[1][2]=X[1]
    t[2][2]=X[2]
    X=pmag.dir2cart((Bdec,Binc,1.0))
    if lower==1 and X[2]<0:
       for i in range(3):
           X[i]=-X[i]
    t[0][0]=X[0]
    t[1][0]=X[1]
    t[2][0]=X[2]
    X=pmag.dir2cart((Gdec,Ginc,1.0))
    if lower==1 and X[2]<0:
       for i in range(3):
           X[i]=-X[i]
    t[0][1]=X[0]
    t[1][1]=X[1]
    t[2][1]=X[2]
# set up v matrix
    v=[0,0,0]
    for i in range(nums):  # incremental point along ellipse
        psi=float(i)*np.pi/xnum
        v[0]=np.sin(beta)*np.cos(psi) 
        v[1]=np.sin(gamma)*np.sin(psi) 
        v[2]=np.sqrt(1.-v[0]**2 - v[1]**2)
        elli=[0,0,0]
# calculate points on the ellipse
        for j in range(3):
            for k in range(3):
                elli[j]=elli[j] + t[j][k]*v[k]  # cartesian coordinate j of ellipse
        PTS.append(pmag.cart2dir(elli))
        R=np.sqrt( 1.-abs(elli[2]))/(np.sqrt(elli[0]**2+elli[1]**2)) # put on an equal area projection
        if elli[2]<0:
#            for i in range(3): elli[i]=-elli[i]
            X_up.append(elli[1]*R)
            Y_up.append(elli[0]*R)
        else:
            X_ell.append(elli[1]*R)
            Y_ell.append(elli[0]*R)
    if plot==1:
        if X_ell!=[]:plt.plot(X_ell,Y_ell,col,linewidth=2,zorder=5)#pylab.plot(X_ell,Y_ell,col,zorder=3)
        if X_up!=[]:plt.plot(X_up,Y_up,col,linewidth=2,zorder=5)#pylab.plot(X_up,Y_up,'g-',zorder=3)
        #pylab.draw()
    else: 
        return PTS
	
def plotCONF(pars): #Modified from PmagPy (Tauxe et al., 2016)
	"""
	plots directions and confidence ellipses 
	"""


	#
	# put on the mean direction
	#
	x,y=[],[]
	XY=pmag.dimap(float(pars[0]),float(pars[1]))
	x.append(XY[0])
	y.append(XY[1])
	#pylab.figure(num=1)
	
	if pars[1]<1:plt.scatter(x, y, edgecolors='m',facecolors='w',marker='*',s=100,zorder=5)
	else:plt.scatter(x, y, edgecolors='w',facecolors='m',marker='*',s=100,zorder=5)
	#pylab.scatter(x,y,marker='^',s=80,c='g',zorder=2)
	#pylab.title(s)
	#
	# plot the ellipse
	#
	plotELL(pars,'m',0,1)

def dokent(data, NN):	#From PmagPy (Tauxe et al., 2016)
    """
    gets Kent  parameters for data ([D,I],N)
    """
    X, kpars = [], {}
    N = len(data)
    if N < 2:
        return kpars
#
#  get fisher mean and convert to co-inclination (theta)/dec (phi) in radians
#
    fpars = pmag.fisher_mean(data)
    pbar = fpars["dec"] * np.pi / 180.
    tbar = (90. - fpars["inc"]) * np.pi / 180.
#
#   initialize matrices
#
    H = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]]
    w = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]]
    b = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]]
    gam = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]]
    xg = []
#
#  set up rotation matrix H
#
    H = [[np.cos(tbar) * np.cos(pbar), -np.sin(pbar), np.sin(tbar) * np.cos(pbar)], [np.cos(tbar)
                                                                                     * np.sin(pbar), np.cos(pbar), np.sin(pbar) * np.sin(tbar)], [-np.sin(tbar), 0., np.cos(tbar)]]
#
#  get cartesian coordinates of data
#
    for rec in data:
        X.append(pmag.dir2cart([rec[0], rec[1], 1.]))
#
#   put in T matrix
#
    T = pmag.Tmatrix(X)
    for i in range(3):
        for j in range(3):
            T[i][j] = old_div(T[i][j], float(NN))
#
# compute B=H'TH
#
    for i in range(3):
        for j in range(3):
            for k in range(3):
                w[i][j] += T[i][k] * H[k][j]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                b[i][j] += H[k][i] * w[k][j]
#
# choose a rotation w about North pole to diagonalize upper part of B
#
    psi = 0.5 * np.arctan(2. * b[0][1] / (b[0][0] - b[1][1]))
    w = [[np.cos(psi), -np.sin(psi), 0],
         [np.sin(psi), np.cos(psi), 0], [0., 0., 1.]]
    for i in range(3):
        for j in range(3):
            gamtmp = 0.
            for k in range(3):
                gamtmp += H[i][k] * w[k][j]
            gam[i][j] = gamtmp
    for i in range(N):
        xg.append([0., 0., 0.])
        for k in range(3):
            xgtmp = 0.
            for j in range(3):
                xgtmp += gam[j][k] * X[i][j]
            xg[i][k] = xgtmp
# compute asymptotic ellipse parameters
#
    xmu, sigma1, sigma2 = 0., 0., 0.
    for i in range(N):
        xmu += xg[i][2]
        sigma1 = sigma1 + xg[i][0]**2
        sigma2 = sigma2 + xg[i][1]**2
    xmu = old_div(xmu, float(N))
    sigma1 = old_div(sigma1, float(N))
    sigma2 = old_div(sigma2, float(N))
    g = -2.0 * np.log(0.05) / (float(NN) * xmu**2)
    if np.sqrt(sigma1 * g) < 1:
        eta = np.arcsin(np.sqrt(sigma1 * g))
    if np.sqrt(sigma2 * g) < 1:
        zeta = np.arcsin(np.sqrt(sigma2 * g))
    if np.sqrt(sigma1 * g) >= 1.:
        eta = old_div(np.pi, 2.)
    if np.sqrt(sigma2 * g) >= 1.:
        zeta = old_div(np.pi, 2.)
#
#  convert Kent parameters to directions,angles
#
    kpars["dec"] = fpars["dec"]
    kpars["inc"] = fpars["inc"]
    kpars["n"] = NN
    ZDir = pmag.cart2dir([gam[0][1], gam[1][1], gam[2][1]])
    EDir = pmag.cart2dir([gam[0][0], gam[1][0], gam[2][0]])
    kpars["Zdec"] = ZDir[0]
    kpars["Zinc"] = ZDir[1]
    kpars["Edec"] = EDir[0]
    kpars["Einc"] = EDir[1]
    if kpars["Zinc"] < 0:
        kpars["Zinc"] = -kpars["Zinc"]
        kpars["Zdec"] = (kpars["Zdec"] + 180.) % 360.
    if kpars["Einc"] < 0:
        kpars["Einc"] = -kpars["Einc"]
        kpars["Edec"] = (kpars["Edec"] + 180.) % 360.
    kpars["Zeta"] = zeta * 180. / np.pi
    kpars["Eta"] = eta * 180. / np.pi
    return kpars
