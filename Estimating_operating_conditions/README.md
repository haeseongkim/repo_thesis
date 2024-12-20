----- operating condition -----
lsq_estimation_multirun_QG_sim_interp_divload_centerline.ipynb
load: 3d data (G specific) (from 1.05 1.05), centerline data (G range)
plot: samples noise ellipse, analytical noise ellipse

display_QG_ellipse_allerrors_sim_centerline.ipynb
load: 3d data (G specific), centerline data (G range) for [ref,Prt,vt]
plot: total UQ ellipse, noise ellipse, physics modeling rhombus

display_QG_determinant_zsensors_large_centerline.ipynb
 suboptimum_sensor_location_f.m
load: centerline data (G range)
plot: normalized Jacobian, best cases n_sensor=2-101, two-sensor contour, histogram n_sensor=2-5

steady_3d_conduction_convection_turb_Adivload_Tqresult.ipynb
load: 3d data (G specific)
plot(fast calculation of rectangle shape): 2d T and q distributions, centerline T and q

steady_3d_conduction_convection_turb_MCMC_centerline.ipynb
 myLinearInterp.m
load: 3d data (G specific), centerline GH table for [ref,Prt,vt]
plot: Bayesian inference chain history, gaussian pdf (bivariate, G, Q)

construction_Amatrix_v1: whole term
construction_Aconv_v1: advection term
 velocity_profile_turb (*1.15)
 eddy_diffusivity_momentun_VanDriest
 Xsteam

construction_divisions
load: A_matrix, A_conv (1.05 1.05 does not exist: use AG)
output1: A_div_L_qbcs1, A_div_A_conv_A_div_L_qbcs1, A_div_A_conv_A_div_L_Tbcs
output2: centerline_A_div_L_qbcs1, centerline_A_div_A_conv_A_div_L_qbcs1, centerline_A_div_A_conv_A_div_L_Tbcs
