import numpy as np
from astropy.io import fits
import os
import sys
import glob



class store_results:
  def __init__(self, obj, MYPATH=None):
	hjd, snr, rv, erv, rvcorr, berv, texp = [], [], [], [], [], [], []
	telfocus, telfocus_scale = [],[]
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
				telfocus.append(np.float(hdr["HIERARCH CAHA TEL FOCU VALUE"]))
				telfocus_scale.append(hdr["HIERARCH CAHA TEL FOCU SCALE"])
		except:
			tmp=2
		
	hjd,rv,erv,snr, berv, texp = np.array(hjd), np.array(rv), np.array(erv), np.array(snr)  , np.array(berv), np.array(texp) 
	telfocus, telfocus_scale = np.array(telfocus), np.array(telfocus_scale)
	self.hjd 	= hjd
	self.rv		= rv		
	self.erv	= erv	
	self.snr	= snr	
	self.rvcorr	= rvcorr	
	self.berv	= berv	
	self.texp	= texp	
	self.telfocus	= telfocus	
	self.telfocus_scale	= telfocus_scale	

class spectrum(file):
  def __init__(self, file, FULL_PATH=False):
	if FULL_PATH:
		mypath = file
	else:
		wpath = os.getcwd()
		mypath = wpath+'/'+file
	a = fits.open(mypath)
		
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
	self.head	= a[0].header

def read_obj(obj,MYPATH=None):
	tmp = store_results(obj, MYPATH=MYPATH)
	return tmp

def read_spec(file,FULL_PATH=False):
	tmp = spectrum(file,FULL_PATH=FULL_PATH)
	return tmp
