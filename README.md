# Visualization of APEX model outputs
Visualizations of biomass outputs from the APEX model and comparison to field data and remote sensing results.

## Data Types
* .SAD = biomass and soil water output file from APEX
  * Biomass should be the sum of STL + STD + GZSLkg/ha + GZSD/ha
  * STL = standing live (Mg/ha)
  * STD = standing dead (Mg/ha)
  * GZSLkg/ha = standing live grazed (kg/ha)
  * GZSDkg/ha = standing dead grazed (kg/ha)
  * A_DDM = plant daily growth after accounting for plant stress (Mg/ha)
  * BIOM = total above and belowground biomass, doesn't include dead material (Mg/ha)
  * X-XCM = volumetric water content between the specified depths

* .CMP = simulated biomass output from APEX
  * Only in PEST version of APEX, each run had a different set of parameter values derived from Montecarlo sampling
  * CPNM2 = plant functional group
  * CAGE_B = cage biomass (kg/ha), ungrazed clipped dry biomass
  * VOR_B = visual obstruction reading, grazed biomass (kg/ha)
  * Biomass = Field collected biomass (kg/ha) from CPER
  * S_BIOMAS = Biomass simulated by APEX
  * WTG = livestock weight gain
     * ADWG = average daily weight gain (kg/animal)
  * Summary statistics are provided for each run and data type
  
* CPER = Central Plains Experimental Range (CPER), Nunn, CO, USA

* CARM Biomass
  * Biomass clipped from grazing exclosure cages (kg/ha) for each plant functional group
  * Cages are moved at the end of each growing season
  
* PastID_ecosite
  * Identifying information for each pasture subarea
  
 * Remote Sensing Data
   * CPER NDVI
   * CPER biomass mean: biomass estimates derived from NDVI by Sean Kearney (lbs/acre)
   
 * Soil Moisture
   * Soil moisture data collected at CPER (volumetric water content [%])
  
