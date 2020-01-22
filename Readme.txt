Dears users,

Some problem have been identified on the previous version of the codes :
"Combined_Calibration_BleachingModel.m"
"Inversion_OSLsurf_time.m"
"Inversion_two_variables.m"
"Inversion_two_variables_t_unknown.m"

The main problem was on the fact that the statistics were made on the entire likelihood distribution, although the way the statistics are made assuming the same frequency of value over the spectrum of values tested, with is not.
Thee new version of the codes resamples the likelihood by keeping the the likelihood value which superior of a random value between 0 and 1. In this way the uncertainties of the likelihood is preserve and the frequency of value over the spectrum of values tested is uniform.

"Inversion_two_variables_t_known_resampling.m"
Invert the values of µ and sigma phi x time from the bleaching luminescence signal. The inversion is made assuming the exposure time in unknown to not propagate the uncertainties of t on the determination of µ. Sigma phi is determined later by dividing sigma phi x time by the time.

"Inversion_two_variables_t_unknown_resampling.m"
Invert the values of µ and sigma phi x time from the bleaching luminescence signal. 

"Inversion_OSLsurface_exposure_age_resampling.m"
Invert the exposure times from the bleaching luminescence signal knowing  µ and sigma phi from calibration. 