# CAFExtractor: a pipeline to make a good CAFE

CAFExtractor (hereafter cafex) is the observatory pipeline of the upgraded CAFE instrument (called CAFE2). This pipeline is meant to be used by the observatory every night. The purpose is that CAFE2 users can leave with their data already reduced.

Version: v0.2

**Version control**
- v0.1	  01/2019	First release
- v0.2	  01/2019	ThAr RV measurements changed
- v0.3	  02/2019	New mask for science RV + fancy CCF plots

**Future improvements**
- Radial velocity uncertainty estimation. Perform a CCF-dependent estimation of the RV uncertainty.
- Order tracing: At high S/N there is a tilt of the extracted orders. To be investigated
- RV dependency with S/N: there is a strong dependency of RV with the S/N of the spectrum. Under investigation.
- Order normalization: Orders #X and #Y have a very bad normalization due to the strong atmospheric absorption at those wavelengths.
- Order merging: to be optimized
- Create IRAF-readable reduced spectra files (also MIDAS if possible)

**Usage**



**Ownership**
These tools have been developed by Dr. Jorge Lillo-Box. While the routines are not public and in the absence of any paper to be cited, any data reduced with the pipeline must contact the pipeline owner (jlillobox@eso.org or jlillobox@gmail.com) to ask for permission and discuss about usage conditions. 
