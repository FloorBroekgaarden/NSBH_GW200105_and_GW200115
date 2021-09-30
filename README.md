# Code to reproduce: "Formation of the First Two Black Hole–Neutron Star Mergers (GW200115 and GW200105) from Isolated Binary Evolution" <https://ui.adsabs.harvard.edu/abs/2021arXiv210805763B/abstract> 
This is a publicly available Github repostitory to recreate all the results and figures from the manuscript "Formation of the First Two Black Hole – Neutron Star Mergers (GW200115 and GW200105) from Isolated Binary Evolution" by Floor S. Broekgaarden and Edo Berger (accepted. to ApJ Letters, <https://ui.adsabs.harvard.edu/abs/2021arXiv210805763B/abstract>). 


For any questions or any errors/bugs/questions with the code or data, please don't hesitate to reach out to Floor Broekgaarden: floor.broekgaarden@cfa.harvard.edu


The repostitory is structured as follows 

-----------------------------------------

## Reproducing Figure 1:

#### Fig 1. Data:

To recreate Figure 1 you will need the following .csv files
 - rates_MSSFR_Models_BHBH_AllDCOsimulation.csv # 22.9 kB
 - rates_MSSFR_Models_BHNS_AllDCOsimulation.csv # 22.9 kB
 - rates_MSSFR_Models_NSNS_AllDCOsimulation.csv # 22.9 kB


that are available in the directory 
https://github.com/FloorBroekgaarden/NSBH_GW200105_and_GW200115/tree/main/dataFiles/summary_data_Fig_1

The same files are also available in the file "csvFilesForFigure2_and_3_DCOpaper.zip" from the Zenodo entry  https://zenodo.org/record/5178777

#### Fig 1. Code:
The code to reproduce Figure 1 is available in the jupyter notebook:
https://github.com/FloorBroekgaarden/NSBH_GW200105_and_GW200115/tree/main/plottingCode/Fig_1/make_Figure_1.ipynb  

-----------------------------------------

## Reproducing Figure 2:

#### *Fig 2. Data:*

To recreate Figure 2 you will need to download all the BHNS data hdf5 files from https://zenodo.org/record/5178777  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5178777.svg)](https://doi.org/10.5281/zenodo.5178777)

This includes all the files below:

 * fiducial.zip,  the Fiducial model (A) and Optimistic CE model (H)
 * massTransferEfficiencyFixed_0_25.zip, the β = 0.25 model (B) 
 * massTransferEfficiencyFixed_0_5.zip, the β = 0.5 model (C) 
 * massTransferEfficiencyFixed_0_75.zip, the β = 0.75 model (D)
 * unstableCaseBB.zip, the unstable case BB mass transfer model (E) and unstable case BB & optimistic CE model (F) 
 * alpha0_1 zip, the α=0.1 model (G) 
 * alpha0_5.zip, the α=0.5 model (I) 
 * alpha2_0.zip, the α=2.0 model (J) 
 * alpha10_0.zip, the α=10.0 model (K) 
 * rapid.zip, the rapid SN model (L) 
 * maxNSmass2_0.zip, the max mNS=2M⊙ model (M) 
 * maxNSmass3_0.zip, the max mNS=3M⊙ model (N)
 * noPISN.zip, the no PISN model (O) 
 * ccSNkick_100km_s.zip, the σcc= 100 km/s model (P) 
 * ccSNkick_30km_s.zip, the σcc= 30 km/s model (Q)
 * noBHkick.zip, the no BH SN kick model (R)
 * wolf_rayet_multiplier_0_1.zip, the model with Wolf-Rayet wind factor fWR=0.1 (S)
 * wolf_rayet_multiplier_5.zip, the model with Wolf-Rayet wind factor fWR=5 (T)



#### *Fig 2. Code:*
The code to reproduce Figure 2 is available in the jupyter notebook:
https://github.com/FloorBroekgaarden/NSBH_GW200105_and_GW200115/tree/main/plottingCode/Fig_2/Chirp_Mass_Matching_Models.ipynb 

-----------------------------------------

## Reproducing Figure 3 and Figure 4:

#### *Fig 3 & 4. Data:*

To recreate Figure 3 & 4 you will need to download the same data as for Fig 2: all the BHNS data hdf5 files from https://zenodo.org/record/5178777  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5178777.svg)](https://doi.org/10.5281/zenodo.5178777)

For recreating Fig. 3 you will also need to download the posterior samples for GW200105 and GW200115, which are publicly available through https://www.gw-openscience.org/about/ (the high spin posteriors from the combined analysis)



#### *Fig 3 & 4. Code:*
The code to reproduce Figure 3 & 4 is available in the jupyter notebook:
https://github.com/FloorBroekgaarden/NSBH_GW200105_and_GW200115/tree/main/plottingCode/Fig_3_and_4/Triangle_plots.ipynb 

-- 





