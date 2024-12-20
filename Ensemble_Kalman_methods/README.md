## transient_1d_conduction_heater_Rgap_insulation.ipynb
 Rgap_T.m
 
 k_T_heater.m
 
 rhocp_T_heater.m
 
 k_T_zro2.m
 
 rhocp_T_zro2.m
 
plot: q(t), R(t), T(t), T(r,t)

## EnKF_augmented_1d_conduction_heater_Rgap_insulation.ipynb
load: (q specific)

plot: q_reconst(t), T_peak_reconst(t), R_gap_reconst(t), T(r,t)

## EnKS_augmented_1d_conduction_heater_Rgap_insulation.ipynb
 Jacobian_expression_augmented_constproperty_f.m
  
load: (q specific)

plot: q_reconst(t), T_peak_reconst(t), R_gap_reconst(t), T(r,t)

## sensitivity_Kalman_methods_sensor_location.ipynb
 EnKF_augmented_1d_conduction_heater_Rgap_insulation_f.m
 
 EnKS_augmented_1d_conduction_heater_Rgap_insulation_f.m
 
load: (q specific)

plot: summary error q_reconst, summary error T_peak_reconst

