# CAFExtractor: a pipeline to make a good CAFE

CAFExtractor (hereafter cafex) is the observatory pipeline of the upgraded CAFE instrument (called CAFE2). This pipeline is meant to be used by the observatory every night. The purpose is that CAFE2 users can have their data reduced by the morning after the observations. CAFExtractor is partly bsaed on the CERES pipeline (Brahm et al., 2017). For additional information on CERES, please visit https://github.com/rabrahm/ceres.

Version: v0.8

##Installation

Just clone this directory to a local path like:

```bash
mkdir /path_to_folder/cafextractor
cd /path_to_folder/cafextractor
git clone https://github.com/jlillo/cafextractor.git
```
Then make sure you tell you computer where the pipeline is:

```tcsh
setenv PYTHONPATH ${PYTHONPATH}:/path_to_folder/cafextractor/src
```

Before starting, make sure you have the following python packages installed:
- python 2.7, numpy, scipy, pyfits, astropy, termcolor, ntpath, argparse


##Usage

Once installed, CAFExtractor is relatively simple to use. A simple example is:

1. Place your data from night YYMMDD on a folder called:

```bash
mkdir /path_to_run_folder/00_RAW/YYMMDD
```

2. Go to CAFExtrcator source folder:

```bash
cd /path_to_cafextractor/cafextractor/src/
```

3. Run CAFExtractor:
```bash
python cafex.py YYMMDD --root /path_to_run_folder/
```

##Example

The CAFExtractor package contains an "example" folder where you can test if the instalation was successful. To run the example just go to the cafextractor/src/ folder and  do:

```bash
python cafex.py 190709 --root /path_to_cafextractor/examples
```

This will create the folder /path_to_cafextractor/examples/11_REDUCED/190709 with the reduction of the science frames in the folder located in /path_to_cafextractor/examples/00_REDUCED/190709

##Version control (summary)
- v0.1	  01/2019	First release
- v0.2	  01/2019	ThAr RV measurements changed
- v0.3	  02/2019	New mask for science RV + fancy CCF plots
- v0.4	  05/2019	New order geometry, new RV mask, RVCORR for SNR, new ref frames
- v0.5	  05/2019	New RV corrections, new --root option, including CCF_FWHM in header
- v0.6	  08/2019	Modified to reduced data from 2018
- v0.7	  -		-
- v0.8	  11/2019	First version publicly released  

##Future improvements

- Order tracing: At high S/N there is a tilt of the extracted orders. To be investigated
- RV dependency with S/N: there is a strong dependency of RV with the S/N of the spectrum. Under investigation.
- Order normalization: Orders #X and #Y have a very bad normalization due to the strong atmospheric absorption at those wavelengths.
- Order merging: to be optimized
- Create IRAF-readable reduced spectra files (also MIDAS if possible)

##Additional resources: 

**RV re-computation**

In some cases, you would like to recompute your RVs because of different reasons. The script RV_measure_onredspec.py on the src folder allows you to do this with different options like including the Simbad-searchable name of your object for a better computation of the BERV in case the OBJECT name put in the header is not resolved by Simbad. Alternatively, you can provide the coordinates of your target. Some examples follows:

```bash
python RV_measure_onredspec.py /path_to_cafextractor/examples/11_REDUCED/190709/reduced/HD109358__190709_0052_red.fits --COORD 12:33:44.54 +41:21:26.92 --RVguess 6.2 --RVampl 100. --UPDATERV
```

**From fits to ascii**

The script fit2ascii_cafex.py also allows to convert the spactra into ascii files:

```bash
python fits2ascii_cafex.py /path_to_cafextractor/examples/11_REDUCED/190709/reduced/
```

##Citation

If you make use of the products of the CAFExtractor pipeline please make sure to cinclude the following reperences:
- Lillo-Box et al., 2019, MNRAS (https://ui.adsabs.harvard.edu/abs/2019arXiv190604060L/abstract)
- Brahm et al., 2017, PASP, 129, 973 (https://ui.adsabs.harvard.edu/abs/2017PASP..129c4002B/abstract)
A suggested sentence to include both references is as follows: "The data were reduced using the CAFExtractor pipeline \citep{lillo-box2019}, partly based on the CERES algorithms \citep{brahm2016}"

##Ownership

These tools have been developed by Dr. Jorge Lillo-Box. While the routines are not public and in the absence of any paper to be cited, any data reduced with the pipeline must contact the pipeline owner (jlillo@cab.inta-csic.es or jlillobox@gmail.com) to ask for permission and discuss about usage conditions. 
