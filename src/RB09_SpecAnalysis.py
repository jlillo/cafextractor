import GLOBALutils
import scipy
import pyfits
import numpy as np
import os
import matplotlib.pyplot as plt
import CAFEutilities
import CAFEx_SetupFile as CS
from astropy.io import fits
from scipy.signal import find_peaks_cwt
from scipy.optimize import curve_fit
from sklearn.metrics.pairwise import euclidean_distances#pairwise_distances
from astroML.stats import sigmaG
import sys
from astropy import constants as c
from astropy.table import Table, Column
from astropy.io import ascii
import sys
import rvtoolbox as rvtbx
import matplotlib.gridspec as gridspec # GRIDSPEC !


def spec_norm(w_frames, cv, frame_names):
	
	NORMdicts = []
	for frame in w_frames:
		wave  = frame[3,:,:]
		flux  = frame[1,:,:]
		eflux = frame[2,:,:]
		fnorm  = flux*0.0 +1.
		efnorm = flux*0.0 +1.
		norm  = flux*0.0 +1.
		myflux_all  = flux*0.0 +1.
		snr_oo = np.zeros(CS.Nominal_Nord)
		for oo in range(CS.Nominal_Nord):
		
			if np.count_nonzero(~np.isnan(flux[oo,:])) != 0:
				# ===== SNR from der_snr
				snr_oo[oo] = CAFEutilities.der_snr(flux[oo,:])
			
				# ===== Flux normalization
				exclude_CR = np.where(flux[oo,:] > np.nanmedian(flux[oo,:])+3.*sigmaG(flux[oo,~np.isnan(flux[oo,:])]))[0]
				exclude_AL = np.where(flux[oo,:] < np.nanmedian(flux[oo,:])-3.*sigmaG(flux[oo,~np.isnan(flux[oo,:])]))[0]
				myflux = flux[oo,:]*1.
				myflux[exclude_CR] = np.nan
				myflux[exclude_AL] = np.nan
				myflux_all[oo,:] = myflux
			
				if np.count_nonzero(~np.isnan(myflux)) != 0:
					try:
						coeff = np.polyfit(wave[oo,~np.isnan(myflux)],myflux[~np.isnan(myflux)],2)
						p = np.poly1d(coeff)
						norm[oo,:] = p(wave[oo,:])
						fnorm[oo,:] = flux[oo,:]/norm[oo,:]
						efnorm[oo,:] = eflux[oo,:]/norm[oo,:]
					except:
						fnorm[oo,:] = flux[oo,:]/np.nanmean(flux[oo,:])
						efnorm[oo,:] = eflux[oo,:]/np.nanmean(flux[oo,:])
						
		
		snr = snr_oo[CS.ordID_5500]
		
# 		for i in range(CS.Nominal_Nord): 
# 			plt.plot(wave[i,200:-200],flux[i,200:-200],c='green')
# 			plt.plot(wave[i,200:-200],myflux_all[i,200:-200],c='red')
# 			plt.plot(wave[i,200:-200],norm[i,200:-200],c='k')			
# 		plt.show()
# 		plt.close()
		
# 		for i in range(CS.Nominal_Nord): 
# 			plt.plot(wave[i,200:-200],fnorm[i,200:-200])
# 		plt.show()
# 		plt.close()

		NORMdict =  { 'fnorm':fnorm,
					'efnorm':efnorm,
					'SNR':snr,
					'SNR_o':snr_oo
					}
		NORMdicts.append(NORMdict)

	return NORMdicts

def merge1D(w_frames, NORMdicts, cv, frame_names):
	
	MERGEdicts = []
	
	for ii,frame in enumerate(w_frames):
		Ndict	= NORMdicts[ii]
		wave  	= frame[3,:,:]
		flux  	= frame[1,:,:]
		eflux 	= frame[2,:,:]
		fnorm  	= Ndict['fnorm']
		efnorm 	= Ndict['efnorm']
		before = 0
		
		wmerge = []
		fmerge = []

		for oo in np.flip(np.arange(CS.Nominal_Nord),0):
			#| Check if order matches next order
			if np.max(wave[oo,200:-200]) > np.min(wave[oo-1,200:-200]):
				disp = np.mean(wave[oo,1:]-wave[oo,0:-1])
				rang = np.max(wave[oo,200:-200]) -  np.min(wave[oo,200:-200])
				wnew = np.linspace(np.min(wave[oo,200:-200]) + before ,np.max(wave[oo,200:-200]), 2048-400)
					
				f0tmp = np.interp(wnew, wave[oo,200:-200], fnorm[oo,200:-200]) 
				f1tmp = np.interp(wnew, wave[oo-1,200:-200], fnorm[oo-1,200:-200], left = np.nan, right=np.nan) 
				wmerge.append(wnew)
				fmerge.append(np.nanmean([f0tmp,f1tmp], axis=0))
				before = np.max(wave[oo,200:-200]) - np.min(wave[oo-1,200:-200])
# 				if oo < 83:
# 					plt.plot(wave[oo,200:-200], flux[oo,200:-200],c='red')
# 					plt.plot(wave[oo-1,200:-200], flux[oo-1,200:-200],c='b')
# 					plt.plot(wnew,np.nanmean([f0tmp,f1tmp], axis=0), c='k' )
# 					plt.show()
# 					sys.exit()
		
		
			else:
				wmerge.append(wave[oo,200:-200])
				fmerge.append(fnorm[oo,200:-200])


		wmerge, fmerge = np.array(wmerge), np.array(fmerge)

		wmerge = np.concatenate(wmerge, axis=0)
		fmerge = np.concatenate(fmerge, axis=0)
				
		MERGEdict =  { 'fmerge':fmerge,
					'wmerge':wmerge,
					}
		MERGEdicts.append(MERGEdict)


	return MERGEdicts






