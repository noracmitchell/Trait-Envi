Guide to Trait-Envi Files:
7/16/2014: cleaned up files and made additional notes to be ready for posting on GitHub


PHYLOGENY:

-FINAL_2011_2012_ultrametric.nex: ultrametricized Pelargonium phylogeny derived from Rigelberg et al. (in prep)

-valenteetree.nex: ultrametricized Protea tree derived from DNA portion of Valente et al. (2010) Protea phylogeny

-schnitzler.nex: ultrametricized Protea tree derived from Schnitzler et al. (2011)

DATA

Pel_CFR_2011_2012: 
-Pelargonium trait data collected from 2011/2012 (sites within the CFR)
-Contains all 4 of our target traits (LMA, FWC, Canopy_area, LWR)

Protea_CFR_DATA:
-contains trait data for all Dimensions Protea collected in the CFR for 2011, 2012, 2013.
-Species_name represents the actual species name we have listed
-Schnitzler_name matches the Schnitzler et al. (2011) phylogeny
-Valente_name matches re-analyzed Valente et al. (2010) phylogeny
-Contains all 4 of our target traits (LMA, FWC, Canopy_area, LWR)

Protea_Pellie_Climate.csv: 
-contains climatic variables for all Protea and Pelargonium sites used in the analysis
-data is for individual population GPS localities
-climate variables are 20 year averages derived from Wilson & Silander (2013)
-elevation and insolation are derived from the ASTER DEM( NASA Land Processes Distributed Active Archive Center)


MODEL FILES AND CODE

traits-environment-with-phylo.R:
-code to run model in R2jags
-settings can be changed to run EITHER Protea or Pelargonium data
-debug mode also an option

traits-enviroment-with-phylo.txt:
-code for model, defining priors, etc.

posterior-comarisons.Rtex, .tex, .pdf, .bib
-code and output for running posterior comparisons to check for differences between Protea and Pelargonium regression coefficients
-LaTeX run through R


RESULTS

results-2014-03-08-12-23-06.txt:
-results from the Pelargonium analysis

results-2014-03-08-13-25-09.Rsave:
-results from the Pelargonium analysis as an Rsave file

results-2014-03-06-08-47-40.txt:
-results from the Protea analysis

results-2014-03-06-09-01-02.Rsave:
-results from the Protea analysis as an Rsave file



