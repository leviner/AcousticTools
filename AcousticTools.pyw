#! /usr/bin/env python
# -*- coding: utf-8 -*-
'''
Acoustic Tools v0.8.5
Things left to do:
- check against Demer


'''

# Python module dependencies:
# General:
import sys
# For Gui:
from PyQt4.QtCore import *
from PyQt4.QtGui import *
#from PyQt4 import QtSql
from ui_AcousticToolsUpdated2_bandalt import Ui_Form
# For saving data files
import csv
from itertools import izip
#For math and figures:
import scipy.special as sps
import numpy as np
#import math as m
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

class TargetStrength(QDialog,Ui_Form):

    def __init__(self, parent=None):
		super(TargetStrength, self).__init__(parent)
		self.setupUi(self)
		# Define button actions
		self.connect(self.calculateTSButton, SIGNAL("clicked()"),self.calcTS)
		self.connect(self.saveFigButton, SIGNAL("clicked()"),self.save_TS_fig)
		self.connect(self.exportDataButton, SIGNAL("clicked()"),self.save_TS_data)
		self.connect(self.targetMaterial,SIGNAL("currentIndexChanged(QString)"),self.target_mat)
		self.connect(self.beamButton,SIGNAL("clicked()"),self.footprint)
		self.connect(self.soundSpeedButton, SIGNAL("clicked()"), self.go2SoundSpeed)
		self.connect(self.soundSpeedButtonBeam, SIGNAL("clicked()"), self.go2SoundSpeed)
		self.connect(self.soundSpeedCalcButton, SIGNAL("clicked()"), self.calculateSoundSpeed)
		self.connect(self.applySSButton, SIGNAL("clicked()"), self.applySS)
		self.connect(self.alphaCalcButton, SIGNAL("clicked()"),self.calculateAlpha)
		self.connect(self.centerFR, SIGNAL("toggled(bool)"), self.TSfrType)
		self.connect(self.bandFR, SIGNAL("toggled(bool)"), self.TSfrType)
		# Set default field values when gui is initialized
		self.targetDensity.setText('14900')
		self.longSoundSpeed.setText('6853')
		self.tranSoundSpeed.setText('4171')
		self.targetDiameter.setText('38.1')
		# Set of necessary class variables
		self.ssplot = []
		self.centerFR.setChecked(True)
		self.TSfrType()

    def TSfrType(self):
        if self.centerFR.isChecked():
            self.minFreq.setEnabled(False)
            self.label_5.setEnabled(False)
            self.maxFreq.setEnabled(False)
            self.label_10.setEnabled(False)
            self.cFreq.setEnabled(True)
            self.label_55.setEnabled(True)
            self.label_6.setEnabled(True)
            self.cFreqTau.setEnabled(True)
            self.label_9.setEnabled(True)
            self.bandTypeTS = 0
        if self.bandFR.isChecked():
            self.cFreq.setEnabled(False)
            self.label_55.setEnabled(False)
            self.label_6.setEnabled(False)
            self.cFreqTau.setEnabled(False)
            self.label_9.setEnabled(False)
            self.minFreq.setEnabled(True)
            self.label_5.setEnabled(True)
            self.maxFreq.setEnabled(True)
            self.label_10.setEnabled(True)
            self.bandTypeTS = 1



    def target_mat(self): # Set the fields for the default material properties for copper or WC targets
        cur_mat = self.targetMaterial.currentText()
        if cur_mat == 'Tungsten-Carbide':
            self.targetDensity.setText('14900')
            self.longSoundSpeed.setText('6853')
            self.tranSoundSpeed.setText('4171')
        elif cur_mat == 'Copper':
            self.targetDensity.setText('8947')
            self.longSoundSpeed.setText('4770')
            self.tranSoundSpeed.setText('2289')

    def calcTS(self): # Calculate the target strength at a given frequency based on given environmental and target material properties
        # preallocate variables
        rhow= None
        cw= None
        densitytarget= None
        cl= None
        ct= None
        DeprecationWarning= None
        freqmin= None
        freqmax = None
        # Get values from GUI
        rhow = self.waterDensity.text()
        rhow = float(rhow)
        cw = self.soundSpeed.text()
        cw = float(cw)
        densitytarget = self.targetDensity.text()   # density of target
        densitytarget = float(densitytarget)
        cl = self.longSoundSpeed.text()                 # longitudinal wave speed of target
        cl = float(cl)
        ct = self.tranSoundSpeed.text()             # transverse wave speed of target
        ct = float(ct)
        D = self.targetDiameter.text()                  # diameters in mm of target
        D = float(D)
        pulse = self.cFreqTau.value()* (10**-6)
        if self.bandTypeTS == 0:
            fnom = self.cFreq.value()
            freqmin = fnom - (1/(pulse * cw))/2
            freqmax = fnom + (1/(pulse * cw))/2
        else:
            freqmin = self.minFreq.value()              # minimum frequency of interest
            freqmin = freqmin
            freqmax = self.maxFreq.value()              # maximum frequency of interest
            freqmax = freqmax

        def wavenumber(freq,c): # determine k, requires nominal frequency and sound speed
            lamb = np.divide(c,freq)
            k = np.divide( np.multiply(2, np.pi), lamb   )
            return (k)

        R = D / (1000 * 2.0)    # converts diameters in mm to radius in m
        Nka = 10000                 # total number of points in the curve
        #densitytarget = 14900.0    # target density in kg/m^3
        #cl = 6848.0                    # longitudinal wave speeds
        #ct = 4161.0                    # transverse wave speed
        hl = float(cl / cw)         # ratio of longitudinal wave speeds
        ht = float(ct / cw)         # ratio of transverse wave speeds
        g = densitytarget/rhow      # density ratio
        angle_inc = 180 # incidence angle, 180 for backscatter
        theta = np.deg2rad(angle_inc)  # converted to rad
        x = np.cos(theta)   #
        ka = np.linspace(0.1,50,Nka)
        Nmax = max(ka) + 10;
        n = range(0, int(Nmax)+1)
        kal= np.divide(ka , hl) # longitudinal wavenumbers
        kat= np.divide(ka , ht) # transverse wavenumbers

        #Legendre
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

        # find range for plot
        plotinds = np.where((ka >= ka_min) & (ka <= ka_max))
        ka_out = ka[plotinds]
        self.f_out = np.divide(np.multiply(ka_out, cw) , (2 * np.pi * R)) / 1000
        self.TS_out = TS[plotinds]
        #ave_TS = sum(self.TS_out)/len(self.TS_out)
        # plot frequency response
        plt = self.TSCurveFigure.canvas.ax
        self.TSCurveFigure.canvas.fig.subplots_adjust(left =.17, bottom =.1)
        plt.clear()
        plt.plot(self.f_out, self.TS_out,'k')
        plt.set_ylabel('TS (dB)')
        plt.set_xlabel('f (kHz)')
        plt.text(.1,  .9 ,  'Mean TS = ' + str(round(np.mean(self.TS_out), 2)) + ' dB', transform = plt.transAxes)
        plt.yaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
        self.TSCurveFigure.canvas.draw()

    def save_TS_fig(self): # Save the TS vs. Frequency plot as a .png
        freqmin = self.minFreq.value()
        freqmax = self.maxFreq.value()
        name = self.targetMaterial.currentText() + '_fmin' + str(freqmin) + '_fmax' + str(freqmax) + '.png'
        self.TSCurveFigure.canvas.print_figure(str(name))

    def save_TS_data(self): # Save the data for the TS vs. Frequency plot as text file
        freqmin = self.minFreq.value()
        freqmax = self.maxFreq.value()
        name = self.targetMaterial.currentText() + '_fmin' + str(freqmin) + '_fmax' + str(freqmax) + '.txt'
        with open(str(name), "wb") as f:
            writer = csv.writer(f)
            writer.writerows(izip(self.f_out, self.TS_out))

    def go2SoundSpeed(self): # Jump to sound speed calculator tab
        self.tabWidget.setCurrentIndex(2)

    def applySS(self): # Set sound speed fields in all tabs equal to result of calculation
        self.soundSpeed.setText(self.soundSpeedResults.text()) # set text in TS tab
        self.soundSpeedBeam.setText(self.soundSpeedResults.text())# set text in beam pattern tab
        self.soundSpeedAlpha.setText(self.soundSpeedResults.text()) # set text in alpha calculation box

    def calculateSoundSpeed(self): # Calculate sound speed from gui fields
        T = float(self.tempSS.text()) # temperature
        S = float(self.salinitySS.text()) # salinity
        D = float(self.depthSS.text()) # depth
        c = self.unesco_eq(S,T,D/100) # calculate sound speed, return c
        self.soundSpeedResults.setText(str(round(c, 3)))
        # calculate a range of sound speed values across a depth range from 0 to 2*D to plot
        D_all = range(int(D)*2)
        c_all = np.zeros(len(D_all))
        for i in D_all:
            c_all[i] = self.unesco_eq(S,T,(float(D_all[i])/100))
        plt = self.soundSpeedFig.canvas.ax
        plt.clear()
        # If this is the first attempt to plot sound speed, invert the y-axis, else, leave inverted
        if self.ssplot != 1:
            plt.invert_yaxis()
            self.ssplot = 1
        plt.plot(c_all,D_all,'k')
        plt.xaxis.set_label_position('top')
        plt.xaxis.set_ticks_position('top')
        plt.set_ylabel('Depth (m)')
        plt.set_xlabel('Sound Speed')
        self.soundSpeedFig.canvas.fig.subplots_adjust(left =.15)
        plt.xaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
        self.soundSpeedFig.canvas.draw()

    def unesco_eq(self, S,T,P): # requires salinity (ppt), temperature (C), and pressure (bars)
        # The unesco equation for calculate sound speed based on
        C_00 = 1402.388;C_01 = 5.03830;C_02 = -5.81090*(10**-2);C_03 = 3.3432*(10**-4);C_04 = -1.47797*(10**-6);C_05 = 3.1419*(10**-9)
        C_10 = .153563;C_11 = 6.8999*(10**-4);C_12 = -8.1829*(10**-6);C_13 = 1.3632*(10**-7);C_14 = -6.126*(10**-10)
        C_20 = 3.126*(10**-5);C_21 = -1.7111*(10**-6);C_22 = 2.5986*(10**-8);C_23 = -2.5353*(10**-10);C_24 = 1.0415*(10**-12)
        C_30  = -9.7729*(10**-9);C_31 = 3.8513*(10**-10);C_32 = -2.3654*(10**-12)
        A_00 = 1.389;A_01 = -1.262*(10**-2);A_02 = 7.166*(10**-5);A_03 = 2.008*(10**-6);A_04 = -3.21*(10**-8)
        A_10 = 9.4742*(10**-5);A_11 = 1.2583*(10**-5);A_12 = -6.4928*(10**-8);A_13 = 1.0515*(10**-8);A_14 = -2.0142*(10**-10)
        A_20 = -3.9064*(10**-7);A_21 = 9.1061*(10**-9);A_22 = -1.6009*(10**-10);A_23 = 7.994*(10**-12)
        A_30 = 1.1*(10**-10);A_31 = 6.651*(10**-12);A_32 = -3.391*(10**-13)
        B_00 = -1.922*(10**-2);B_01 = -4.42*(10**-5)
        B_10 = 7.3637*(10**-5);B_11 = 1.795*(10**-7)
        D_00 = 1.727*(10**-3)
        D_10 = -7.9836*(10**-6)

        D_tp = D_00 + (D_10*P)
        B_tp = B_00 + (B_01*T) + ((B_10+(B_11*T))*P)
        A_tp = (A_00 + (A_01*T) + (A_02 * (T**2)) + (A_03 * (T**3)) + (A_04 * (T**4))) +\
                ((A_10 + (A_11*T) + (A_12*(T**2)) + (A_13*(T**3)) + (A_14*(T**4)))*P) +\
                ((A_20 + (A_21*T) + (A_22*(T**2)) + (A_23*(T**3))) * (P**2)) +\
                ((A_30 + (A_31*T) + (A_32*(T**2)))*(P**3))
        Cw_tp = (C_00 + (C_01*T) + (C_02*(T**2)) + (C_03*(T**3)) + (C_04*(T**4)) + (C_05*(T**5))) +\
                ((C_10 + (C_11*T) + (C_12*(T**2)) + (C_13*(T**3)) + (C_14*(T**4)))*P) + \
                ((C_20 + (C_21*T) + (C_22*(T**2)) + (C_23*(T**3)) + (C_24*(T**4)))*(P**2)) + \
                ((C_30 + (C_31*T) + (C_32*(T**2)))*(P**3))
        c = Cw_tp + (A_tp*S) + (B_tp*(S**(3/2))) + (D_tp*(S**2))
        return c

    def calculateAlpha(self): # Calculate alpha based on values in GUI fields
        c = float(self.soundSpeedAlpha.text()) # sound speed
        pH = float(self.phAlpha.text()) # acidity
        T = float(self.tempAlpha.text()) # temperature
        D = float(self.depthAlpha.text()) # depth
        S = float(self.salinityAlpha.text()) # salinity
        f = float(self.freqAlpha.text()) # nominal frequency

        alpha = self.FGequation(c, pH, T, D, S,f)/1000 # calculate and return alpha, in dB/meter
        self.alphaResults.setText(str(round(alpha, 5)))

        # calculate a range of slpha values across a frequency bandwidth from 0 to 2*f to plot
        F_all = range(int(f)*2)
        alpha_all = np.zeros(len(F_all))
        for i in F_all:
            alpha_all[i] = self.FGequation(c, pH, T, D, S,F_all[i])/1000
        plt = self.alphaFig.canvas.ax
        plt.clear()
        plt.plot(F_all,alpha_all,'k')
        plt.set_ylabel('alpha (dB/m)')
        plt.set_xlabel('Frequency (kHz)')
        self.alphaFig.canvas.fig.subplots_adjust(left =.15)
        self.alphaFig.canvas.draw()

    def FGequation(self, c, pH, T, D, S,f): # requires sound speed (m/s), pH, temp(C), depth(m), salinity(ppt), and nominal frequency(kHz)
        # Attenuation Coefficient is based on Francois and Garrison, 1982 - "Sound absorption based on ocean measurements.
        # Boric Acid Contribution, P1 = 1
        A1=((8.86/c)*(10**(0.78*pH-5)))
        f1=((2.8*((S/35)**0.5))*(10**(4-(1245/(T+273)))))
        # MgSO4 Contribution
        A2=((21.44*(S/c))*(1+(0.025*T)))
        P2=(1-(1.37*(10**-4)*D)+(6.2*(10**-9)*(D**2)))
        f2=((8.17*(10**(8-(1990/(T+273)))))/(1+.0018*(S-35)))
        # Pure water contribution, where A3 is temperature dependent
        if T > 20:
            A3=((3.964*(10**-4))-(1.146*(10**-5)*T)+(1.45*(10**-7)*(T**2))-(6.5*(10**-10)*(T**3)))
        else:
            A3=((4.937*(10**-4))-(2.59*(10**-5)*T)+(9.11*(10**-7)*(T**2))-(1.5*(10**-8)*(T**3)))
        P3=((1-(3.83*(10**-5)*D)) + (4.9*(10**-10)*(D**2)))
        # Calculate and return Alpha
        alpha = (((f**2)*A1*f1)/(((f1**2)) + (f**2)))+ ((A2*P2*f2*(f**2))/((f2**2) + (f**2))) + (A3*P3*(f**2))
        return alpha

    def footprint(self):
		fmin = self.minFreqBeam.value()
		fmax = self.maxFreqBeam.value()
		fc = self.footprintFreq.value()
		c = float(self.soundSpeedBeam.text())
		ran = self.targetRangeBeam.value()
		din = self.xducerSizeBeam.value()
		unit  = self.xducerSizeUnitBeam.currentText()

		if unit == 'cm':
			dia = din/100
		elif unit == 'in':
			dia = din*2.54/100

		if self.beamType3.isChecked() == True:
			thresh = .5
		elif self.beamType43.isChecked() == True:
			thresh = 1 / np.exp(1)

		num_pts = fmax - fmin      # calculate number of points to include both endpoints
		f = np.linspace(fmin,fmax,int(num_pts)+1,endpoint=False)  # frequency at 1 kHz spacing
		f *= 1000                   # Convert frequency to Hz for calculations
		fc *= 1000                  # Convert center frequency to Hz for calculations
		wavlen = c / f              # wavelength
		k = 2 * np.pi / wavlen      # wavenumber
		r = dia/2                   # divide into radius for calculations
		wavlenfc = c / fc
		kfc = 2 * np.pi / wavlenfc
		theta = np.linspace(0.1,70,3500,endpoint=False) # create array with 0.05 deg. spacing
		bw = []
		for i in range(len(f)):
			[D, bout] = self.beamwidth(thresh,theta, k[i], r)
			bw.append(bout)
		fresnel_r = np.sqrt(wavlenfc * ran / 2) * 100     # only first Fresnel radius
		[dtmp, bwfc] = self.beamwidth(thresh, theta, kfc, r)
		fpfc = r * 100 + ran * 100 * np.tan(np.pi * np.array(bwfc) / 180)
		fp = r * 100 + ran * 100 * np.tan(np.pi * np.array(bw) / 180)
		ff = r ** 2 / wavlenfc              # range to farfield

		plt = self.beamFootprintFig.canvas.ax
		plt.clear()
		plt.plot(f / 1000, bw,'k')
		plt.set_ylabel('Beamnwidth (deg)')
		plt.set_xlabel('Frequency (kHz)')
		self.beamFootprintFig.canvas.draw()

		plt = self.beamFootprintFig_2.canvas.ax
		plt.clear()
		plt.plot(f / 1000, fp,'k')
		plt.set_xlabel('Frequency (kHz)')
		plt.set_ylabel('Footprint Rad. (cm)')
		self.beamFootprintFig_2.canvas.draw()

		amp = 10 * np.log10( D )
		amp = np.hstack([amp[::-1], amp])
		trad = theta * np.pi / 180
		trad = np.hstack([-trad[::-1], trad])
		plt = self.beamPatternFig.canvas.ax
		self.beamPatternFig.canvas.figure.clf()
		plt.cla()
		ax1 = plt.figure.gca(polar=True)
		ax1.clear()
		ax1 = plt.figure.gca(polar=True)
		ax1.set_theta_zero_location('N')
		#ax1.set_theta_direction('clockwise')
		ax1.set_ylim(-80, 0)
		ax1.set_yticks(np.arange(-80,0,20))

		thetaticks = [0,45, 90, 135,180,225,270,315]
		ax1.set_thetagrids(thetaticks)

		#ax1.set_xticks(np.pi/180. * np.linspace(0,360, 8, endpoint=False))
		#labels = ['1','2','3','4','5','6','7','8']
		#ax1.set_xticklabels(labels)
		ax1.plot(trad, amp, 'k')
		self.beamPatternFig.canvas.draw()


    def beamwidth(self,threshold,angle,wavenumber, radius):
        bw_sintheta = np.sin( np.array(angle) * np.pi / 180)
        bw_krsintheta = wavenumber * radius * np.array(bw_sintheta)
        D = np.power(( 2 * sps.jn(1,bw_krsintheta) / bw_krsintheta ),2)
        D[0] = 1 # D has NAN in first indices, replace with 1 because D = 1
        bwind=np.argmin(abs(D-threshold)) # find index for closest value to beamwidth
        bwout = angle[bwind]        # theta for beamwidth
        return (D,bwout)


# Main, runs if called from command line
if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    form = TargetStrength()
    form.show()
    sys.exit(app.exec_())
