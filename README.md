# semicomp_dosefinding
Code (main branch) and data (ResultsData branch) associated with paper "A semi-competing risks framework to incorporate dose reductions and dose interruptions into oncology dose-finding"

Files beginning with ‘simulation_’ are the main simulation files that generate data, fit with given model, and write results files. These were run on UVA’s Rivanna HPC system. Most of these were split into 10 identical copies and run n=100 iterations, to be able to run in parallel and get n=1000 iterations for each simulation scenario.
The label attached to ‘gen’ in these files tells from what types of hazards the data is generated from. The options are constant (const, small increasing (weibinc), small decreasing (weibsd), and big decreasing (weibbd) as described in the paper. The last piece of the file name is what model (and prior) was used to analyze the data. ‘constanalysis’ fits with the Constant-Skeleton model. ‘weibanalysis’ uses the Weibull-Skeleton model. ‘titecrm’ uses the TITE-CRM. ‘TITECRMMC’ uses the TITE-CRMMC. The CRM-MC is created after the fact with results from the TITE-CRMMC. The ‘uninform’ prior is the uninformative normal prior, ‘inform’ is informative normal, and ‘uniformprior’ is uniform prior. All instances of Weibull-Skeleton analysis use ‘uninform’ prior. 

The respective model and prior combination are implemented in Stan via one of four files: ‘stancode_constant_uninform.txt’, ‘stancode_constant_inform.txt’, ‘stancode_constant_uniformprior.txt’, or ‘stancode_weibull_uninform.txt’. These files are called within the simulation files when analysis occurs. 

All of the simulation files that use one of the proposed models applies an equal-weighted Total Distance dosing decision rule during the trial. The Maximum Uniformly Tolerated Dose is calculated in subsequent analysis files from results produced. 

The files beginning with ‘cgen’, ‘wbdgen’, ‘wincgen’, or ‘wsdgen’ are the same simulation files as above, just using the 3:1 Total Distance weighting for dosing decisions (DLT to DR) instead of equal weighting. 

‘TITECRMMC_analysis.R’ uses the TITE-CRMMC analysis and creates results for both TITE-CRMMC and CRMMC. 

‘titecrm_analysis.R’ generates analysis from TITE-CRM results. 

‘weib_analysis.R’ generates analysis from the results of the Weibull-Skeleton simulations (all generating hazards, uninformative prior).

‘weibanalysis_3_1.R’ generates same analysis, but when 3:1 weighted Total Distance dosing decision rule is used. 

‘const_const_analysis.R’ generates analysis from the results of the Constant-Skeleton simulations (constant generating hazards, 3 priors).

‘weibbd_const_analysis.R’ generates same analysis (Big Decreasing hazards, 3 priors). 

‘weibinc_const_analysis.R’ generates same analysis (Small Increasing hazards, 3 priors). 

‘weibsd_const_analysis.R’ generates same analysis (Small Decreasing hazards, 3 priors). 

‘genericprobcalculations.R’ was used to manually determine dose-toxicity curves based on the data generating models. It plots dose-toxicity curves as well. 

‘resultsviz.R’ plots all PCS, Overdose, and Accuracy Indices for proposed models and comparator models (CRMs) as found in the paper. 

‘weib_accuracyanalysis_finitesample.R’ plots bias, variance, and MSE for Constant-Skeleton and Weibull-Skeleton model as found in the paper. 

The simulation results files needed for the analysis files are found in the ResultsData branch of this repository. The main set of analysis data gathered into one document is ‘SummaryResults_AllMethods.xlsx’
