# GOF_Test
Goodness-of-fit (GOF) test for cox model under monotone constraint implemented by R.

The `time_dep` folder contains (updating) programs for simulation study of test for time-dependent covariate:
- `BootstrapHypothesisTesting_SurvFixKM.R` is the main function, which defines the simulation settings. The underlying hazard is defined in this program.
- `Bootstrap_exp_SurvFixKM.R` is called by the main function to conduct the test. A sample is generated by the program using the given hazard. Survival time is generated by given hazard when baseline hazard follows exponential(1) distribution.
- `CenTimeFix_Td.R` is used to generate censoring time in bootstrap sample by Kaplan-Meier estimator conditional on new censoring time is larger or equal to the observed survival time. New censoring time is only generated for uncensored subjects. This function is called by `Bootstrap_exp_SurvFixKM`.
- `SurvTime_td.R` is used to generate survival time in bootstrap sample using the inverse method from the `Breslow-Type` estimator of the distribution function of the cumulative baseline hazard.
- `isoph_td.R` is used to estimate the monotone hazard using pseudo-iterative convex minorant algorithm proposed by Chung 2018.

The `time-indep` folder contains programs for simulation study for time-independent covariate:
- `BootstrapHypothesisTesting_SurvFixKM.R` is the main function, which defines the simulation settings. The underlying hazard is defined in this program. `Bootstrap_exp_SSurvFixKM.R` is called by the main function and calls `CenTimeFix.R`, `SurvTimeSmoothSurv.R` and `nlpl_Censor_ind.R` (calculating the negative log partial likelihood with monotone hazard estimator or Cox hazard). `nlpl_noCensor_ind.R` can be used to evaluate the negative log partial likelihood when there is no censoring.
-   `SimulationResult_BHT02SSurvFixKMmaxX1000s500b.rda` is the simulation result based on 500 Monte Carlo samples of size 1000 while the critical value is obtained by bootstrapping 500 times. `readSimulationResult.R` is used to read the output.
-   `Data Example.Rmd` contains the programs to apply the test to some data examples and information of the datasets.
-   `MartingaleResidualMethod` folder of the similar structure contains the result for residual-based GOF test. 

