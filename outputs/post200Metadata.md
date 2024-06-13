Metadata for "./outputs/post200.csv", the analysis dataset used to fit statistical models.
Columns include variable name, type, use in statistical models (Response, fixed or random) and description.

Additional columns in post200.csv with an "S" subscript are scaled variables for modelling. These are identical to the non-S variant of the column, but scaled to a mean of 0 and a standard deviation of 1. These are not described in the metadata.

| Variable name  | Type     | Modelled? | Description |
| ------------- | ------------- | ------ | --------- | 
| site | factor |   | Neotoma site ID |
| novelID| factor | Random | Neotome site ID and sampling bin of detected novel community. Acts as a unique ID for each novel state. |
| bin | numeric |  | Center of 200-year sampling bin used to average pollen sample data. |
| dP | numeric |  | Bray-Curtis dissimilarity from the *pre-novel* state (see methods for details).|
| dN | numeric |  | Bray-Curtis dissimilarity from the *novel* state (see methods for details).| 
| delta.dP | numeric |  | Change in dP from previous post-novel state. Not modelled or used in analysis.|
| delta.dN | numeric |  | Change in dN from previous post-novel state. Not modelled or used in analysis.|
| dP1 | numeric |  | Change in dP from first post-novel state. Not modelled or used in analysis.|
| dN1 | numeric |  | Change in dN from first post-novel state. Not modelled or used in analysis.|
| novel.time | numeric | | Center of 200-year sampling bin where the novel state was detected. Forms part of novelID.
| time.since.novel | numeric | Fixed | Time between the center of post-novel sampling bin and center of novel state sampling bin.|
| bin.n | numeric | | Count of the position of novel state along the Neotoma time series, as as the position of 200-year sampling bins with sample data (e.g., "4" = fourth sampling bin along the time series ordered oldest to youngest).  |
| H0 | numeric | | Hill number order 0 of post-novel pollen composition data. Equates to plant family count. |
| H1 | numeric | | Hill number order 1 of post-novel pollen composition data. Equivalent to the exponent of Shannon's Entropy. |
| H2 | numeric | | Hill number order 2 of post-novel pollen composition data. Equivalent to the inverse of Simpson's Index. |
| rareDivN | numeric | | Sample N for "rareDiv", set to the minimum of 500 or total pollen grains in the sample. |
| rareType | numeric | | Rarefaction adjusted diversity of post-novel pollen composition data. Rarefaction set to 10 pollen grains.
| rareDiv | numeric | | Rarefaction adjusted diversity of post-novel pollen composition data. Rarefaction set to rareDivN. |
| novelH0 | numeric | | Hill number order 0 of novel state pollen composition data. Equates to plant family count. |
| novelH1 | numeric | | Hill number order 1 of novel state pollen composition data. Equivalent to the exponent of Shannon's Entropy. The natural-log transformation of this variable, "novelH1L", below was used as a fixed effect in models. |
| novelH2 | numeric | | Hill number order 2 of novel state pollen composition data. Equivalent to the inverse of Simpson's Index. |
| novelRareType | numeric | | Rarefaction adjusted diversity of novel state pollen composition data. Rarefaction set to 10 pollen grains. |
| novelRareDiv | numeric | | Rarefaction adjusted diversity of novel state pollen composition data. Rarefaction set to rareDivN. |
| preNovH0 | numeric | | Hill number order 0 of pre-novel state pollen composition data. Equates to plant family count. |
| preNovH1 | numeric | | Hill number order 1 of pre-novel state pollen composition data. Equivalent to the exponent of Shannon's Entropy. |
| preNovH2 | numeric | | Hill number order 2 of pre-novel state pollen composition data. Equivalent to the inverse of Simpson's Index. |
| preNovRareType | numeric | | Rarefaction adjusted diversity of pre-novel state pollen composition data. Rarefaction set to 10 pollen grains. |
| preNovRareDiv | numeric | | Rarefaction adjusted diversity of pre-novel state pollen composition data. Rarefaction set to rareDivN. |
| novAbund | numeric | | Total pollen grain count across all samples averaged into the 200-year sampling bin of the novel state. |
| novAbundProp | numeric | | novAbund expressed as a proportion of the average across all 200-year sampling bins in the time series. Values > 1 indicate more pollen than average in the novel sampling bin, and vice versa.|
downerN | numeric | | Count of taxa that decreased in abundance from the preceding state (the novel state for the first post-novel observation, or the preceding post-novel state for each successive observation.|
upperN | numeric | | Count of taxa that increased in abundance from the preceding state (the novel state for the first post-novel observation, or the preceding post-novel state for each successive observation.|
| relAbundChange | numeric | | Total % of relative abudance change from the preceding state |
| gamma | numeric | Fixed | Time series gamma diversity. The total number of plant families detected across the entire study period.|
| tslength | numeric | Fixed | The total number of 200-year sampling bins across the time series.|
| cat.bef | factor | |  Excess column from novel community detection framework that is not used here.
| seq.dist | numeric | | Measures of "instantaneous dissimilarity" as per *Pandolfi, J. M., Staples, T. L., & Kiessling, W. (2020). Increased extinction in the emergence of novel ecological communities. Science, 370(6513), 220–222. https://doi.org/10.1126/science.abb3996* |
| raw.min.dist | numeric | | Measures of "cumulative dissimilarity" as per *Pandolfi, J. M., Staples, T. L., & Kiessling, W. (2020). Increased extinction in the emergence of novel ecological communities. Science, 370(6513), 220–222. https://doi.org/10.1126/science.abb3996* |
| seq.exp | numeric | | Modelled expectation of instantaneous dissimilarity|
| min.exp | numeric | | Modelled expectation of cumulative dissimilarity|
min.p | numeric | Fixed | Probability of cumulative dissimilarity falling within modelled expectations. Lower = greater magnitude of compositional novelty. 1 - min.p was used as a covariate to describe the "magnitude of novelty". |
| bin.lag | numeric | Fixed | Time difference between novel state and pre-novel state. Most are 200 years, but some time series had empty sampling bins.| 
| temp | numeric | | Global temperature estimate for the post-novel state, from interpolated data (see Methods). All post-novel observations in the same 200-year bin, regardless of site, have the same value. |
| temp.sd | numeric | | Standard deviation of "temp" global temperature estimate.|
| temp.se | numeric | | Standard error of "temp" global temperature estimate.|
| localTemp | numeric | | Regional temperature estimate from TraCE gridded data. Post-novel observations in the same 200-year bin and in the same TraCE grid cell will have the same value. |
globalLag[1, 5, 10, 15, 20, 25] | numeric | Fixed | Variants of global temperature change, using different past temperature estimates to calculate temperature change. Number following variable refers to the number of 200-year bins used as a lag, so 1 = 200 years, 5 = 1000 years, ..., 25 = 5000 years.|
localLag[1, 5, 10, 15, 20, 25] | numeric | Fixed | Variants of regional temperature change, using different past temperature estimates to calculate temperature change. Number following variable refers to the number of 200-year bins used as a lag, so 1 = 200 years, 5 = 1000 years, ..., 25 = 5000 years.|
| lat | coordinates | | Site latitude in decimal degrees (WGS84). |
| long | coordinates | | Site longitude in decimal degrees  (WGS84). |
| novelH1L | numeric | Fixed | Natural-log transformation of Hill order 1 of the novel state. |
| H1L | numeric | | Natural-log transformation of the Hill order 1 of the post-novel observation. |
| preNovH1L | numeric | | Natural-log transformation of the Hill order 1 of the pre-novel state. |
| deltaH1 | numeric | Fixed | novelH1 - preNovH1L, a measure of change in ln-transformed diversity between the pre-novel state and the novel state. Positive values mean diversity increased in the novel state and vice versa. |
| deltaH1abs | numeric | Fixed | Absolute version of deltaH1. Used as a model variant to discriminate predictive power between directionality and magnitude of diversity change. |
| deltaPost | numeric | Fixed | H1L - novelH1L, a measure of change in ln-transformed diversity between the novel state and the post-novel observation. Positive values mean diversity increased in the post-novel observation and vice versa. |
| deltaPostabs | numeric | Fixed | Absolute version of deltaPost. Used as a model variant to discriminate predictive power between directionality and magnitude of diversity change. |
| tsProp | numeric | Fixed | binN / tsLength, the proportion along the time series where the novel state occurred. Used as a covariate to test for sampling effects in early vs late novel states, correcting for different length of time series.|
| novAbundPropAbs | numeric | Fixed | Absolute variant of novAbundProp, where values below 1 were inverted. This was used as a covariate to test whether variation in reponse variables could be explained by small vs large pollen samples in the novel state, or simply greater deviation from average (directionality vs magnitude)|
| dMag | numeric | Response | "Turnover" measurement of rotated dN and dP values for each post-novel observation. Used as response variable for turnover models. See methods for details.|
| dRatio | numeric | Response | "Persistence" measurement of rotated dN and dP values for each post-novel observation. Used as response variable for persistence models. See methods for details.|