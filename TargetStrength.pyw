#! /usr/bin/env python

'''
Created: 9 Oct 2014 - CBassett
Last Edited: 13 Oct 2014 - RLevine (QtGui and MPL canvas within window)

This codes is a version of the standard target code.
This version is only suitable for backscatter (theta = 180 degrees)
The main difference between the standard target code "wInputs" is the GUI


'''
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4 import QtSql
from ui_TargetStrengthfin import Ui_mainWindow
import sys

#For the math things:
import scipy.special as sps
import numpy as np
import math as m
import matplotlib.pyplot as plt


class TargetStrength(QMainWindow,Ui_mainWindow):

	def __init__(self, parent=None):        
		super(TargetStrength, self).__init__(parent)
		self.setupUi(self)

		
		self.connect(self.Calculate_but, SIGNAL("clicked()"),self.run)
		self.connect(self.save_fig_but, SIGNAL("clicked()"),self.save_fig)
		self.connect(self.target_material,SIGNAL("currentIndexChanged(QString)"),self.target_mat)



	def target_mat(self):
		cur_mat = self.target_material.currentText()
		print cur_mat
		if cur_mat == 'Tungsten-Carbide':
			self.target_density.setText('14900')
			self.target_ts_long.setText('6853')
			self.target_ts_trans.setText('4171')
		elif cur_mat == 'Copper':
			self.target_density.setText('8947')
			self.target_ts_long.setText('4770')
			self.target_ts_trans.setText('2289')
	

	def run(self):
		rhow= None
		cw= None
		densitytarget= None
		cl= None
		ct= None
		DeprecationWarning= None
		freqmin= None
		freqmax = None
	
		rhow = self.water_density.text()
		rhow = float(rhow)
		cw = self.water_soundspeed.text()
		cw = float(cw)
		densitytarget = self.target_density.text()	# density of target
		densitytarget = float(densitytarget)
		print densitytarget
		cl = self.target_ts_long.text()   				# longitudinal wave speed of target
		cl = float(cl)		
		ct = self.target_ts_trans.text()				# transverse wave speed of target
		ct = float(ct)
		D = self.target_diameter.text()					# diameters in mm of target
		D = float(D) 					
		freqmin = self.min_freq.value()				# minimum frequency of interest
		freqmin = freqmin
		freqmax = self.max_freq.value()				# maximum frequency of interest
		freqmax = freqmax


 	
		# This loop checks the minimum and maximum frequencies
		# If they are switched (f_min > f_max), this switches them
		if freqmin > freqmax:
			fm1 = freqmin
			freqmin = freqmax
			freqmax = fm1

		def wavenumber(freq,c):
			lamb = np.divide(c,freq)
			k = np.divide( np.multiply(2, np.pi), lamb   )
 			return (k)

		R = D / (1000 * 2.0)  	# converts diameters in mm to radius in m
		Nka = 10000 				# total number of points in the curve
		densitytarget = 14900.0 	# target density in kg/m^3
		cl = 6848.0					# longitudinal wave speeds
		ct = 4161.0    			    # transverse wave speed
		hl = float(cl / cw)   		# ratio of longitudinal wave speeds
		ht = float(ct / cw)			# ratio of transverse wave speeds
		g = densitytarget/rhow   	# density ratio
	 
		angle_inc = 180 # incidence angle, 180 for backscatter
		theta = np.deg2rad(angle_inc)  # converted to rad
		x = np.cos(theta)   #

		ka = np.linspace(0.1,50,Nka)

		Nmax = max(ka) + 10;
		n = range(0, int(Nmax)+1)


		kal= np.divide(ka , hl) # longitudinal wavenumbers
		kat= np.divide(ka , ht) # transverse wavenumbers

		##########
		#Legendre
		##########
		pn, dpnit = sps.lpn(max(n),x) 

		jn = np.zeros((Nka,Nmax+1))
		djn = np.zeros((Nka,Nmax+1))
		yn = np.zeros((Nka,Nmax+1))
		dyn = np.zeros((Nka,Nmax+1))
		jnl = np.zeros((Nka,Nmax+1))
		djnl = np.zeros((Nka,Nmax+1))
		jnt = np.zeros((Nka,Nmax+1))
		djnt = np.zeros((Nka,Nmax+1))

		for j in range(len(ka)):
	
			xj, dxj = sps.sph_jn(Nmax,ka[j])
			xy, dxy = sps.sph_yn(Nmax,ka[j])
			xl, dxl = sps.sph_jn(Nmax,kal[j])
			xt, dxt = sps.sph_jn(Nmax,kat[j])
	
			jn[j,:] = xj
			djn[j,:] = dxj
			yn[j,:] = xy
			dyn[j,:] = dxy
			jnl[j,:] = xl
			djnl[j,:] = dxl
			jnt[j,:] = xt
			djnt[j,:] = dxt


		#next loop 
		nl = np.multiply(2.0, n)+ 1 # summation term
		f = np.zeros(len(ka))
		nn =np.array(n) * np.array(n) + n

		for j in range(len(ka)):
	
			tan1 = -np.divide(np.multiply(kal[j], djnl[j,:]), jnl[j,:])
			tan2 = -np.divide(np.multiply(kat[j], djnt[j,:]), jnt[j,:])
			tan3 = -np.divide(np.multiply(ka[j], djn[j,:]), jn[j,:])
			tan_beta = -np.divide(np.multiply(ka[j], dyn[j,:]), yn[j,:])
			tan_del =  -np.divide(jn[j,:], yn[j,:])
			d1 = tan1 + 1.0
			d2 = nn - 1.0 - kat[j] ** 2 / 2.0 + tan2;
			term1a = np.divide(tan1, d1)
			term1b = np.divide(nn, d2)
			term2a = np.divide( nn - kat[j] ** 2 / 2.0 + 2.0*tan1, d1)
			term2b = np.divide( np.multiply(nn , tan2 + 1), d2)
			td = np.multiply(-0.5 * kat[j] ** 2, np.divide ( (term1a - term1b), (term2a - term2b)))
			tan_phi = np.divide(-td, g)
			tan_eta = np.multiply(tan_del, np.divide( tan_phi + tan3 , tan_phi + tan_beta ))
			cos_eta = np.divide( 1.0 , (1 + np.multiply(tan_eta, tan_eta)) ** 0.5 )
			sin_eta = np.multiply(tan_eta, cos_eta)
			bn = np.multiply(sin_eta, 1j * cos_eta) -  np.multiply(sin_eta, sin_eta)
			s = np.multiply(np.multiply(nl, pn),bn) 
			esum =   sum(s)  
			f[j] = abs(esum)

		theta = np.deg2rad(angle_inc)
		# this is the solution for f/ka
	
		fka = np.divide(f,ka)
		#F = cw.* np.divide(theta,2*np.pi*R);
		TS = 10.0 * np.log10( np.divide(fka, 2) * R **2)

		# find ka_min
		ka_min = wavenumber(freqmin,cw)* R * 1000
		ka_max = wavenumber(freqmax,cw)* R * 1000

# print "The ka range is %.1f to %.1f" %  (ka_min, ka_max)

		# find range 
		plotinds = np.where((ka >= ka_min) & (ka <= ka_max))
		ka_out = ka[plotinds]
		f_out = np.divide(np.multiply(ka_out, cw) , (2 * np.pi * R)) / 1000

		TS_out = TS[plotinds]
		
		
		plt = self.widget.canvas.ax
		plt.clear()
		
		# plt.ion()
		#plt.figure(1, figsize=(10,10))
		plt.plot(f_out, TS_out,'k')
		plt.set_ylabel('TS (dB)')
		plt.set_xlabel('f (kHz)')
		# #plt.savefig('StandardTarget_fromGUI.png', bbox_inches='tight')
		# plt.draw()
		self.widget.canvas.draw()
		
	def save_fig(self):
		freqmin = self.min_freq.value()	
		freqmax = self.max_freq.value()	
		name = self.target_material.currentText() + '_fmin' + str(freqmin) + '_fmax' + str(freqmax) + '.png'
		print name
		self.widget.canvas.print_figure(str(name))
		
		
		
# main, runs all from command line	
if __name__ == "__main__":	
    import sys
    app = QApplication(sys.argv)
    form = TargetStrength()
    form.show()
    sys.exit(app.exec_())