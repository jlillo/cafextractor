import numpy as np
from astropy.io import fits
import os
import sys
import glob



class store_results:
  def __init__(self, obj, MYPATH=None):
	hjd, snr, rv, erv, rvcorr, berv, texp = [], [], [], [], [], [], []
	if MYPATH is not None:
		wpath = MYPATH
	else:
		wpath = os.getcwd()
	for file in glob.glob(wpath+'/*.fits'):
		a = fits.open(file)
		hdr = a[0].header
		try:
			if hdr["OBJECT"].split(' ')[0] == obj:
				hjd.append(hdr["HIERARCH CAFEX HJD"])
				rv.append(hdr["HIERARCH CAFEX RV"])
				erv.append(hdr["HIERARCH CAFEX ERV"])
				snr.append(hdr["HIERARCH CAFEX SNR"])
				berv.append(np.float(hdr["HIERARCH CAFEX BERV"]))
				texp.append(np.float(hdr["EXPTIME"]))
				rvcorr.append(hdr["HIERARCH CAFEX RVCORR"])
		except:
			tmp=2
		
	hjd,rv,erv,snr, berv, texp = np.array(hjd), np.array(rv), np.array(erv), np.array(snr)  , np.array(berv), np.array(texp) 
	
	self.hjd 	= hjd
	self.rv		= rv		
	self.erv	= erv	
	self.snr	= snr	
	self.rvcorr	= rvcorr	
	self.berv	= berv	
	self.texp	= texp	

class spectrum(file):
  def __init__(self, file):
	wpath = os.getcwd()
	a = fits.open(wpath+'/'+file)
		
	self.wave	= a[2].data
	self.flux	= a[1].data	
	self.eflux	= a[3].data
	self.fnorm	= a[4].data
	self.efnorm	= a[5].data	
	self.wmerge	= a[6].data
	self.fmerge	= a[7].data
	self.dvel	= a[8].data
	self.ccf	= a[9].data
	self.eccf	= a[10].data

def read_obj(obj,MYPATH=None):
	tmp = store_results(obj, MYPATH=MYPATH)
	return tmp

def read_spec(file):
	tmp = spectrum(file)
	return tmp
