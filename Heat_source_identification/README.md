## steady_3d_conduction_convection_turb_single_qsource.ipynb
load: A matrix (G specific) (from AG)

plot (general heat source shape): 2d T and q distributions, centerline T and q, vertices

example outputs for IRN calculation (including 'IT_conversion_A_div_I_conversion'):

 steady_3d_conduction_convection_turb_v1_single_qsource_rectangle_G4215_0852kW
 
 steady_3d_conduction_convection_turb_v1_single_qsource_halfcross_G4215_0852kW
 
 steady_3d_conduction_convection_turb_v1_single_qsource_arc1_G4215_0852kW
 
 steady_3d_conduction_convection_turb_v1_single_qsource_groove1_G4215_0852kW

## IRN_automu_sim_reinitialize.ipynb
single IRN solver (with reinitialization and adaptive mu)
 smoothed_holder_weights.m
 
load: (G, q, and shape specific)

plot: initial guess, converged result, iteration statistics, edges

## IRN_automu_sim_reinitialize_runs.ipynb
 IRN_automu_sim_reinitialize_f.m
 
load: (G, q, and shape specific)

plot: MAE, regularization, estimated noise, edge feature, area feature, original shape

## IRN_Bayes_hybrid.ipynb
 MRF_line_process.m
 
load: (G, q, and shape specific) & Tqbcs (G and q specific)

plot: IRN solution, chain history, Bayes solution, heat flux statistics

## Other codes
construction_all_pointsources_allqbcs

load: A matrix (G specific)

output for Bayes calculation: T_qbcs_allqbcs (temperature sensitivities to point heat sources)
