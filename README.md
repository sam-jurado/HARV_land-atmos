# Increasing large precipitation events and low available water holding capacity create the conditions for dry land-atmosphere feedbacks in the Northeastern United States

## PIs
- Sam Jurado, PhD. Student, Yale University
- Jackie Matthes, Senior Scientist, Harvard University

## Project Summary
As a warmer climate enables an increase in atmospheric humidity, extreme precipitation events have become more frequent in the Northeastern United States. Understanding the impact of evolving precipitation patterns is critical to understanding water cycling in temperate forests and moisture coupling between the atmosphere and land surface. Although the role of soil moisture in evapotranspiration has been extensively studied, few have analyzed the role of soil texture in determining ecosystem-atmosphere feedbacks. In this study, we utilized long term data associated with ecosystem water fluxes to deduce the strength of land-atmosphere coupling at Harvard Forest, Petersham, MA, USA. We found a 1.5% increase in heavy precipitation contribution per decade where high-intensity events compose upwards of 42% of total yearly precipitation in 2023. Intensifying precipitation trends were found in conjunction with a long-term soil drying at the Harvard Forest despite no significant increase in evapotranspiration over 32 years. This suggests that soil water holding capacity is a key mediating variable controlling the supply of water to ecosystems and the atmosphere. We found that these land surface changes directly impacted the lifted condensation level (LCL) height over Harvard Forest which was found to be increasing at a rate of 6.62 meters per year while atmospheric boundary layer (ABL) heights have fallen at a modest rate of 1.76 meters per year. This has amplified dry feedbacks between the land surface and the atmosphere such that 80% of observed summers ending in a water deficit also had an anomalously low soil water content in the spring.

## Data
- Hava a google drive of data in the same way in which I have it on my computer so it is easier to run directly from this code?
- Dataframe structure and units information?


## Code
The following is the reccomended order in which to run the code.
1. **Soil_Moisture_Final.R** This code creates Figure 5a. Additionally, the dataframes df_anom and df_anom_daily are created, which are necessary to run Evap_Perc.R and ABL_LCL_Final.R.
2. **Evap_Perc.R** This code creates Figure 6.
3. **ABL_LCL_Final.R** This code creates Figures 7,8 and 9.

The following code can be run in any order, and is not dependent on the files created in Soil_Moisture_Final.R
- **NEUS_pcp_sm.R** This code creates Figure 2.
- **Precip_HF.R** This code creates Figure 10.
- **SL_T_Final.R** This code creates Figure 4.
- **HF-soilmoisture.R** This code creates Figure 3.
- **HF-streamflow.R** This code creates Figure 5b.


## Deprecated Folder
- This folder contains old code no longer needed to recreate the figures found in the manuscript.
