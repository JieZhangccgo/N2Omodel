To run the model, experimental information should be firstly created as a matlab file 'expInfo.mat'.
The 'expInfo.mat' should include the following information:

***The matrices mentioned below have 12 rows indicating 12 variables in the sequence of: (1) CO2, (2) DOC (mmolC/L), (3) O2, (4) b_AER, (5) NH4, (6) NO2, (7) NO3, (8) N2O, (9) b_AOB, (10) N2, (11) b_DEN, and (12) b_NOB.

*msInfo.icMesh: an array with the size m [0,...,0.1], the depths of soil samples in a 10 cm soil core. [m] 
*msInfo.icData: a 12*m matrix, the initial concentration profiles along the soil depth for the 12 variables. [mmol/L] 
*msInfo.thetaMat: a 12*m matrix, the water-filled porosity or air-filled porosity profiles along the soil depth for the 12 variables. A value of 1 should be assigned to the rows of microbes. [-]  
*msInfo.DeffMat: a 12*m matrix, the effective diffusion coefficient profiles along the soil depth for the 12 variables. A value of 0 should be assigned to the rows of microbes. [m2/d] 
*msInfo.theta_tot: a number, total porosity in the soil core. [-]  
*msInfo.N_NH4: a number, the dimensionless Freundlich isotherm exponent 0.74 in the paper  
*msInfo.KF_NH4: a number, the dimensionless Freundlich isotherm exponent 4.89 in the paper
*msInfo.rhob: a number, bulk density [kg/L soil]
*msInfo.restSOC_t0: a number, SOC content in the soil. [gC/g soil]
*msInfo.POC_t0: a number, POC content in the manure. [gC/g soil]
*msInfo.biomass_base: an array with the size of four, the base values of four microbial biomass: AER, AOB, NOB and DEN in order. [g biomass/g soil]

***See the BG paper (Zhang et al., 2023) for more details. 