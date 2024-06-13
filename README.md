Code used to generate results, figures and tables for Proceedings B submission:* Novel ecological change drives new compositional trajectories*

The data used to generate this publication were obtained from the Neotoma Paleoecological Database: details are located in dataProcessing.R.

Running scripts will require a copy of the entire repository, including subfolders, which contain raw data and supplementary functions. Please add correct working directory entries in toplevel .R scripts and relative file paths will do the rest.

R scripts make extensive use of code-folding functionality in RStudio.

REPOSITORY TREE:

**./dataProcessing.R**: commented code to access Neotoma data and harmonize taxonomy. Raw outputs have been included in repository so this script does not need to be re-run to reproduce results.

**./nove\_plants\_legacy.R**: commented code that reproduces analyses, figures and tables found in the main text.

**./functions**: Small ease-of-use functions read into top-levels scripts.

**./rawdata**: Input data objects.
* 	**./rawdata/Neotoma vascular plant records.rds**: storage of raw Neotoma download. See Neotoma docs for metadata.
* 	**./rawdata/processedRecords.rds**: Neotoma data with harmonized taxonomy. See Neotoma docs for metadata.
* 	**./rawdata/analysisRecords.rds**: Neotoma data subset and transformed, ready for analysis. Notably, these data include pollen counts corrected via Relative Pollen Production values (relative to Poaceae).
* 	**./rawdata/WFOTaxaTable.csv**: Reference table for taxa (unique "variablenames" in Neotoma records), with abundance ("count") and World Flora harmonized taxonomy for family, genus and species.
* 	**./rawdata/WFOTaxaTableEDITED.csv**: Version with manual adjustment of some taxa where World Flora could not parse records.
*  **./rawdata/RPP\_Dataset\_v2\_Table\_5\_6.xlsx**: Relative pollen production values from *M. Wieczorek, U. Herzschuh, Compilation of relative pollen productivity (RPP) estimates and taxonomically  harmonised RPP datasets for single continents and Northern Hemisphere extratropics. Earth Syst Sci Data 12, 3515–3528 (2020).*
* 	**./rawdata/climate/**: Climatic data used to generate temperature change variables. FIles are raw downloads from references or data repositories,
*	**/rawdata/climate/Shakun2012\_retreat\_temp.csv:** 6500 - 22000 ybp, global temperature estimates from *J. D. Shakun, et al., Global warming preceded by increasing carbon dioxide concentrations during the last deglaciation. Nature 484, 49–54 (2012).*
*	**/rawdata/climate/temp12k\_allmethods\_percentiles.csv:** 0 - 12000 ybp, global temperature estimates from from *D. Kaufman, et al., Holocene global mean surface temperature, a multi-method reconstruction approach. Sci Data 7, 201 (2020).*
*	**/rawdata/climate/trace.01-36.22000BP.clm2.TSA.22000BP\_decavg\_400BCE.nc**: 0 - 22000 ybp temperature estimates, gridded at a "regional" grid resolution of ~3.6&deg;. From the TRacE project: *F. He, P. U. Clark, Freshwater forcing of the Atlantic Meridional Overturning Circulation re-visited. Nat Clim Chang 12, 449–454 (2022).*

**./plots**: Figure outputs from novel_plants_legacy.R, used in main text and online supplement.

**./outputs**: File outputs from novel_plants_legacy.R. Includes intermediate analysis steps, including novel community detection framework ("all neotoma novelty.rds"), analysis data-frame ("post200.csv"), as well as model summary tables.
**./outputs/post200.csv**: Dataset of values used for model selection process described in methods. Separate metadata is available as "post200Metadata.md".

