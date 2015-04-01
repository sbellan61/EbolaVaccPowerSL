# Organization of Repository

## Incidence Fits
1. data Folder

* Data Ebola (Public).xlsx - raw data downloaded from UNMEER's humanitarian data exchange.
* slData.R processes this file, cleans it, and produces cleanSLData.Rdata as a long-format data table of district-level Ebola virus disease (EVD) incidence in Sierra Leone.
* createHT.Rdata is produced by Plot_ExpFit.R (see below) and stores exponential model fits to this data.

2. ExpFit.R - Functions to fit exponential models with negative binomial likelihoods to district-level EVD incidence data.

3. Plot_ExpFits.R - Perform & plot Sierra Leone district-level exponential model fits along with incidence projections. Plot example cluster-level hazard trajectories and individual-level variation around cluster means.

## Trial Simulations

1. simFuns.R - Functions to set up a trial population, with time-varying hazards by cluster and individual hazard variation, to set up individuals' vaccination schedules for a specified study design, and then simulate infections.

2. PlotTrial.R - Plot simulated hazard trajectories (either from Sierra Leone projections, or phenomenological model).

3. AnalysisFuns.R - Functions to take as input a simulated trial object and return as output data objects of the form analyzable by survival analysis or by cluster-level Poisson regressions. Also includes functions to truncate/censor data based on trial design and analysis, plot trial designs, compile output statistics of interest, and finally to perform an entire trial simulation and analysis from start to finish.

4. CoxFxns.R - Functions to analyze simulated trial data with Cox proportional hazard gamma frailty models, Poisson regressions, and permutation or bootstrap approaches over either of those methods. Also includes functions to run GEEs and GLMMs, which exhibited convergence problems or caused segmentation faults and were subsequently excluded from the manuscript.

6. startSim.R - Initialize simulation runs with parameters fed in from R CMD BATCH.

7. Files starting with MK - make TXT files with each line an R CMD BATCH call to startSim.R, to be sent to a high performance computing (HPC) cluster.

8. SLSumm.R - Collects results from HPC cluster (too big to store in repository), and summarizes each parameter set's performance over 2040 simulated trials in a data table.

9. Files starting with Plot - Plot results from output summaries produced by SLSumm.R. Also makes tables, and extracts numbers for the manuscript text. Sources plotting functions in multiplot.R and ggplot themes in ggplotTheme.R.

10. trialDiagrams.R - plots trial diagrams (e.g. Figs 3, S1, S2)


