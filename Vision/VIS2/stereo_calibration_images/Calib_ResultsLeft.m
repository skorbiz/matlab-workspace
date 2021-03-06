% Intrinsic and Extrinsic Camera Parameters
%
% This script file can be directly excecuted under Matlab to recover the camera intrinsic and extrinsic parameters.
% IMPORTANT: This file contains neither the structure of the calibration objects nor the image coordinates of the calibration points.
%            All those complementary variables are saved in the complete matlab data file Calib_Results.mat.
% For more information regarding the calibration model visit http://www.vision.caltech.edu/bouguetj/calib_doc/


%-- Focal length:
fc = [ 832.511895249810891 ; 838.884609640151211 ];

%-- Principal point:
cc = [ 510.674221352745576 ; 413.892896292243620 ];

%-- Skew coefficient:
alpha_c = 0.000000000000000;

%-- Distortion coefficients:
kc = [ -0.366370249282006 ; 0.215078117967467 ; -0.009280324909704 ; -0.000268874992000 ; 0.000000000000000 ];

%-- Focal length uncertainty:
fc_error = [ 2.849551743702734 ; 2.098171191715270 ];

%-- Principal point uncertainty:
cc_error = [ 3.380205640200274 ; 5.325378719204483 ];

%-- Skew coefficient uncertainty:
alpha_c_error = 0.000000000000000;

%-- Distortion coefficients uncertainty:
kc_error = [ 0.007394974991379 ; 0.041193413372192 ; 0.001944288226494 ; 0.000388724769101 ; 0.000000000000000 ];

%-- Image size:
nx = 1024;
ny = 768;


%-- Various other variables (may be ignored if you do not use the Matlab Calibration Toolbox):
%-- Those variables are used to control which intrinsic parameters should be optimized

n_ima = 25;						% Number of calibration images
est_fc = [ 1 ; 1 ];					% Estimation indicator of the two focal variables
est_aspect_ratio = 1;				% Estimation indicator of the aspect ratio fc(2)/fc(1)
center_optim = 1;					% Estimation indicator of the principal point
est_alpha = 0;						% Estimation indicator of the skew coefficient
est_dist = [ 1 ; 1 ; 1 ; 1 ; 0 ];	% Estimation indicator of the distortion coefficients


%-- Extrinsic parameters:
%-- The rotation (omc_kk) and the translation (Tc_kk) vectors for every calibration image and their uncertainties

%-- Image #1:
omc_1 = [ 2.050959e+00 ; 9.049245e-01 ; -3.166605e-01 ];
Tc_1  = [ -1.337331e+02 ; -6.839603e+01 ; 6.751729e+02 ];
omc_error_1 = [ 4.683689e-03 ; 2.949744e-03 ; 4.622037e-03 ];
Tc_error_1  = [ 2.751388e+00 ; 4.147641e+00 ; 2.378311e+00 ];

%-- Image #2:
omc_2 = [ 1.843380e+00 ; 1.772321e+00 ; -6.226471e-01 ];
Tc_2  = [ -1.186995e+02 ; -1.736331e+02 ; 7.181988e+02 ];
omc_error_2 = [ 2.894319e-03 ; 3.805099e-03 ; 5.462680e-03 ];
Tc_error_2  = [ 2.956683e+00 ; 4.209280e+00 ; 2.559888e+00 ];

%-- Image #3:
omc_3 = [ 1.961552e+00 ; 1.903528e+00 ; -4.577139e-01 ];
Tc_3  = [ -1.223776e+02 ; -2.095567e+02 ; 6.518657e+02 ];
omc_error_3 = [ 2.296173e-03 ; 3.533060e-03 ; 5.561164e-03 ];
Tc_error_3  = [ 2.713487e+00 ; 3.723560e+00 ; 2.351849e+00 ];

%-- Image #4:
omc_4 = [ 2.050091e+00 ; 1.997150e+00 ; -3.165673e-01 ];
Tc_4  = [ -1.178702e+02 ; -2.317174e+02 ; 5.979372e+02 ];
omc_error_4 = [ 2.229778e-03 ; 3.543904e-03 ; 5.924751e-03 ];
Tc_error_4  = [ 2.516437e+00 ; 3.365641e+00 ; 2.121425e+00 ];

%-- Image #5:
omc_5 = [ 2.065080e+00 ; 2.013532e+00 ; -2.847116e-01 ];
Tc_5  = [ -1.131274e+02 ; -1.561238e+02 ; 5.718193e+02 ];
omc_error_5 = [ 2.227369e-03 ; 3.172376e-03 ; 5.761141e-03 ];
Tc_error_5  = [ 2.357878e+00 ; 3.359315e+00 ; 1.955830e+00 ];

%-- Image #6:
omc_6 = [ 2.112237e+00 ; 1.390929e+00 ; -7.828956e-02 ];
Tc_6  = [ -1.388191e+02 ; -8.492964e+01 ; 6.792660e+02 ];
omc_error_6 = [ 4.119251e-03 ; 2.755189e-03 ; 5.000803e-03 ];
Tc_error_6  = [ 2.781326e+00 ; 4.177379e+00 ; 2.223202e+00 ];

%-- Image #7:
omc_7 = [ 2.140752e+00 ; 1.050511e+00 ; 1.466618e-01 ];
Tc_7  = [ -1.736580e+02 ; -6.398354e+01 ; 6.429487e+02 ];
omc_error_7 = [ 4.994381e-03 ; 2.564142e-03 ; 4.682467e-03 ];
Tc_error_7  = [ 2.646759e+00 ; 3.999642e+00 ; 2.094184e+00 ];

%-- Image #8:
omc_8 = [ 2.137336e+00 ; 1.114766e+00 ; 1.050565e-01 ];
Tc_8  = [ -8.815910e+01 ; -4.167023e+01 ; 6.050018e+02 ];
omc_error_8 = [ 4.883021e-03 ; 2.393553e-03 ; 4.864660e-03 ];
Tc_error_8  = [ 2.472118e+00 ; 3.772479e+00 ; 1.937673e+00 ];

%-- Image #9:
omc_9 = [ 1.690216e+00 ; 2.324676e+00 ; -7.779764e-01 ];
Tc_9  = [ -9.436244e+01 ; -1.891757e+02 ; 8.454011e+02 ];
omc_error_9 = [ 1.628070e-03 ; 4.021194e-03 ; 6.574615e-03 ];
Tc_error_9  = [ 3.466055e+00 ; 5.017973e+00 ; 2.892101e+00 ];

%-- Image #10:
omc_10 = [ 1.468808e+00 ; 2.555215e+00 ; -9.817536e-01 ];
Tc_10  = [ -4.409646e+01 ; -2.034015e+02 ; 8.591948e+02 ];
omc_error_10 = [ 1.557252e-03 ; 4.401188e-03 ; 7.306879e-03 ];
Tc_error_10  = [ 3.526480e+00 ; 5.113586e+00 ; 2.797171e+00 ];

%-- Image #11:
omc_11 = [ 1.370071e+00 ; 2.470213e+00 ; -1.034127e+00 ];
Tc_11  = [ -3.477934e+01 ; -2.029600e+02 ; 9.340383e+02 ];
omc_error_11 = [ 1.703588e-03 ; 4.575725e-03 ; 6.928574e-03 ];
Tc_error_11  = [ 3.828174e+00 ; 5.556342e+00 ; 3.125004e+00 ];

%-- Image #12:
omc_12 = [ 1.652383e+00 ; 2.169593e+00 ; -7.662157e-01 ];
Tc_12  = [ -1.505040e+02 ; -1.718256e+02 ; 8.705145e+02 ];
omc_error_12 = [ 1.845674e-03 ; 4.140012e-03 ; 6.015602e-03 ];
Tc_error_12  = [ 3.560139e+00 ; 5.187767e+00 ; 3.129356e+00 ];

%-- Image #13:
omc_13 = [ 2.013000e+00 ; 1.453366e+00 ; -2.223372e-01 ];
Tc_13  = [ -8.755542e+00 ; -1.583324e+02 ; 8.586700e+02 ];
omc_error_13 = [ 4.062441e-03 ; 3.059104e-03 ; 5.159992e-03 ];
Tc_error_13  = [ 3.514474e+00 ; 5.163841e+00 ; 2.886665e+00 ];

%-- Image #14:
omc_14 = [ 2.030592e+00 ; 1.390492e+00 ; -1.806673e-01 ];
Tc_14  = [ 1.340755e+01 ; -9.788254e+01 ; 7.603816e+02 ];
omc_error_14 = [ 4.221435e-03 ; 2.798881e-03 ; 5.124948e-03 ];
Tc_error_14  = [ 3.097198e+00 ; 4.644053e+00 ; 2.593489e+00 ];

%-- Image #15:
omc_15 = [ 1.853997e+00 ; 1.117411e+00 ; -5.680809e-01 ];
Tc_15  = [ -1.562287e+01 ; 1.423338e+01 ; 8.664304e+02 ];
omc_error_15 = [ 4.657556e-03 ; 3.448827e-03 ; 4.674996e-03 ];
Tc_error_15  = [ 3.518430e+00 ; 5.526510e+00 ; 2.942487e+00 ];

%-- Image #16:
omc_16 = [ 1.288165e+00 ; 1.922299e+00 ; -1.036097e+00 ];
Tc_16  = [ -6.987047e+01 ; -5.572467e+01 ; 9.813928e+02 ];
omc_error_16 = [ 3.345478e-03 ; 4.810940e-03 ; 5.551215e-03 ];
Tc_error_16  = [ 3.982413e+00 ; 6.117120e+00 ; 3.233787e+00 ];

%-- Image #17:
omc_17 = [ 1.016208e+00 ; 2.152012e+00 ; -1.150343e+00 ];
Tc_17  = [ -1.635374e+01 ; -6.397285e+01 ; 9.863187e+02 ];
omc_error_17 = [ 3.034986e-03 ; 5.026003e-03 ; 5.842098e-03 ];
Tc_error_17  = [ 4.001015e+00 ; 6.131730e+00 ; 3.225257e+00 ];

%-- Image #18:
omc_18 = [ 1.509028e+00 ; 1.692615e+00 ; -9.206471e-01 ];
Tc_18  = [ -5.424097e+01 ; -6.327684e+00 ; 7.362728e+02 ];
omc_error_18 = [ 3.612225e-03 ; 4.334331e-03 ; 5.144586e-03 ];
Tc_error_18  = [ 2.983501e+00 ; 4.657328e+00 ; 2.362110e+00 ];

%-- Image #19:
omc_19 = [ 1.984417e+00 ; 8.752529e-02 ; 1.601769e-02 ];
Tc_19  = [ -6.996981e+01 ; 2.865917e+01 ; 4.548626e+02 ];
omc_error_19 = [ 5.508232e-03 ; 2.721861e-03 ; 3.979025e-03 ];
Tc_error_19  = [ 1.853711e+00 ; 2.948954e+00 ; 1.647712e+00 ];

%-- Image #20:
omc_20 = [ 1.992981e+00 ; 1.737983e-01 ; 2.270472e-01 ];
Tc_20  = [ -6.981466e+01 ; -3.733053e+01 ; 4.217839e+02 ];
omc_error_20 = [ 5.612189e-03 ; 3.020948e-03 ; 3.979440e-03 ];
Tc_error_20  = [ 1.718161e+00 ; 2.594100e+00 ; 1.589632e+00 ];

%-- Image #21:
omc_21 = [ 1.967883e+00 ; 6.286332e-02 ; 2.231306e-01 ];
Tc_21  = [ -1.374155e+02 ; 1.325873e+02 ; 4.773501e+02 ];
omc_error_21 = [ 5.358578e-03 ; 2.604391e-03 ; 3.893085e-03 ];
Tc_error_21  = [ 2.018178e+00 ; 3.348399e+00 ; 1.889726e+00 ];

%-- Image #22:
omc_22 = [ 1.939145e+00 ; 1.614892e-01 ; 3.861389e-01 ];
Tc_22  = [ -1.190530e+02 ; 8.012676e+01 ; 4.555505e+02 ];
omc_error_22 = [ 5.517152e-03 ; 2.925346e-03 ; 3.901118e-03 ];
Tc_error_22  = [ 1.903318e+00 ; 3.076144e+00 ; 1.753048e+00 ];

%-- Image #23:
omc_23 = [ 1.522175e+00 ; 2.046105e+00 ; -9.768817e-01 ];
Tc_23  = [ -1.307406e+02 ; -1.732792e+02 ; 8.072065e+02 ];
omc_error_23 = [ 2.422080e-03 ; 4.638664e-03 ; 6.046251e-03 ];
Tc_error_23  = [ 3.310921e+00 ; 4.784035e+00 ; 2.846743e+00 ];

%-- Image #24:
omc_24 = [ 1.684034e+00 ; 1.687190e+00 ; -7.646449e-01 ];
Tc_24  = [ -9.801423e+01 ; -9.641899e+01 ; 7.271680e+02 ];
omc_error_24 = [ 3.429782e-03 ; 4.069773e-03 ; 5.263192e-03 ];
Tc_error_24  = [ 2.957741e+00 ; 4.421591e+00 ; 2.485667e+00 ];

%-- Image #25:
omc_25 = [ 2.115681e+00 ; 7.312601e-01 ; -2.371467e-01 ];
Tc_25  = [ -1.320841e+02 ; -4.785151e+01 ; 6.454793e+02 ];
omc_error_25 = [ 4.771627e-03 ; 2.790738e-03 ; 4.595935e-03 ];
Tc_error_25  = [ 2.628382e+00 ; 4.001465e+00 ; 2.294085e+00 ];

