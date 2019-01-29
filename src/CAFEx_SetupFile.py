
# ==============================
# Basic paths:
# root		=	'/Users/lillo_box/00_Instrumentation/CAFE/CAFExtractor/test_data/'
# raw			=	'/Users/lillo_box/00_Instrumentation/CAFE/CAFExtractor/test_data/00_RAW/'
# redfolder	=	'/Users/lillo_box/00_Instrumentation/CAFE/CAFExtractor/test_data/12_REDUCED/'
root		=	'/Volumes/willyfog/gang5/jlillo/22_RUNS/2019_01_CAHA_2.2_CAFE_Recommisioning_Run2/'
raw			=	'/Volumes/willyfog/gang5/jlillo/22_RUNS/2019_01_CAHA_2.2_CAFE_Recommisioning_Run2/00_RAW/'
redfolder	=	'/Volumes/willyfog/gang5/jlillo/22_RUNS/2019_01_CAHA_2.2_CAFE_Recommisioning_Run2/11_REDUCED/'
RefFrames	=	'/Users/lillo_box/00_Instrumentation/CAFE/CAFExtractor/cafextractor/ReferenceFrames/'

RefArc		=	'arc__180718_0031.fits'
RefFlat		=	'flat__180718_0011.fits'
# ==============================

# Range (in days) to group the different calibration frames 
# to perfom the MasterFrame (MasterBias, MasterFlat, MasterArc).
# RB03.
 
bandwidth=	{	'bias': 1.0/24.,
		      	'flat': 1.0/24.,
		      	'arcs': 1.0/24.,
		      	'sci':	1.
		     }

MinMembersBias = 4
MinMembersFlat = 4
MinMembersArcs = 9

# Sigma clipping value. Any frame with a mean above sigclip times the std dev   
# will be rejected to perform the MasterFrame (MasterBias, MasterFlat, MasterArc)
# RB03.

sigclip=	{	'bias': 5.,
		   		'flat': 0.01, # 1% of the mean, just to avoid wrongly classified flats 
		   		'arcs': 5.,
		    }

# ==============================
# RB04: Tracing orders

order_aperture_ampl = 5
order_trace_degree  = 4

y0Nominal_first = 115.433 #135.9
Nominal_Nord = 84		# Number of orders to be exrtacted
ordID_5500 = 43			# Order corresponding to 5500A
order0 = 60				# Order corresponding to first extracted order



# ==============================
# RB05: Wavelength calibraiton
nx	= 3				# X-order of the polynomial for wavelength calibration
nm  = 8				# Y-order of the polynomial for wavelength calibration
#Selected_MasterARC = 'evening' 	# Default MasterArc to use for wavelength calibration
Selected_MasterARC = 0 	# Default MasterArc to use for wavelength calibration


# ==============================

Marsh_alg       = 0
ext_aperture    = 5
NSigma_Marsh    = 5
NCosmic_Marsh   = 5
S_Marsh         = 0.4
N_Marsh         = 3      # grado polinomio 
min_extract_col = 64
max_extract_col = 1984
RO_fl 			= 3.3/20.0
GA_fl 			= 1.0
npools 			= 1
