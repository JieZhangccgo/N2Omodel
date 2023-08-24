To run the model, experimental information should be firstly created as a matlab file 'expInfo.mat'.
The 'expInfo.mat' should include the following information:

***The matrices mentioned below have 12 rows indicating 12 variables in the sequence of: (1) CO2, (2) DOC (mmolC/L), (3) O2, (4) b_AER, (5) NH4, (6) NO2, (7) NO3, (8) N2O, (9) b_AOB, (10) N2, (11) b_DEN, and (12) b_NOB.

*expInfo.icMesh: an array with the size m [0,...,0.1], the depths of soil samples in a 10 cm soil core. [m] 
*expInfo.icData: a 12*m matrix, the initial concentration profiles along the soil depth for the 12 variables. [mmol/L] 
*expInfo.thetaMat: a 12*m matrix, the water-filled porosity or air-filled porosity profiles along the soil depth for the 12 variables. A value of 1 should be assigned to the rows of microbes. [-]  
*expInfo.DeffMat: a 12*m matrix, the effective diffusion coefficient profiles along the soil depth for the 12 variables. A value of 0 should be assigned to the rows of microbes. [m2/d] 
*expInfo.theta_tot: a number, total porosity in the soil core. [-]  
*expInfo.N_NH4: a number, the dimensionless Freundlich isotherm exponent 0.74 in the paper  
*expInfo.KF_NH4: a number, the dimensionless Freundlich isotherm exponent 4.89 in the paper
*expInfo.rhob: a number, bulk density [kg/L soil]
*expInfo.restSOC_t0: a number, SOC content in the soil. [gC/g soil]
*expInfo.POC_t0: a number, POC content in the manure. [gC/g soil]
*expInfo.biomass_base: an array with the size of four, the base values of four microbial biomass: AER, AOB, NOB and DEN in order. [g biomass/g soil]

***See the BG paper doi.org/10.5194/bg-2023-98, Zhang et al. (2023) for more details. 
