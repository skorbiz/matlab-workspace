% Intrinsic and Extrinsic Camera Parameters
%
% This script file can be directly excecuted under Matlab to recover the camera intrinsic and extrinsic parameters.
% IMPORTANT: This file contains neither the structure of the calibration objects nor the image coordinates of the calibration points.
%            All those complementary variables are saved in the complete matlab data file Calib_Results.mat.
% For more information regarding the calibration model visit http://www.vision.caltech.edu/bouguetj/calib_doc/


%-- Focal length:
fc = [ 1397.708892253199792 ; 1289.751144358146348 ];

%-- Principal point:
cc = [ 467.886537607122477 ; 199.759487970284482 ];

%-- Skew coefficient:
alpha_c = 0.000000000000000;

%-- Distortion coefficients:
kc = [ -0.463010917060898 ; 0.190483472116130 ; 0.053177791113432 ; 0.006255473244304 ; 0.000000000000000 ];

%-- Focal length uncertainty:
fc_error = [ 49.477437741326952 ; 84.751884309421541 ];

%-- Principal point uncertainty:
cc_error = [ 28.082796436891510 ; 161.688843542243802 ];

%-- Skew coefficient uncertainty:
alpha_c_error = 0.000000000000000;

%-- Distortion coefficients uncertainty:
kc_error = [ 0.104172669036319 ; 0.190103613596390 ; 0.050796956222575 ; 0.005355643024477 ; 0.000000000000000 ];

%-- Image size:
nx = 1024;
ny = 768;


%-- Various other variables (may be ignored if you do not use the Matlab Calibration Toolbox):
%-- Those variables are used to control which intrinsic parameters should be optimized

n_ima = 4;						% Number of calibration images
est_fc = [ 1 ; 1 ];					% Estimation indicator of the two focal variables
est_aspect_ratio = 1;				% Estimation indicator of the aspect ratio fc(2)/fc(1)
center_optim = 1;					% Estimation indicator of the principal point
est_alpha = 0;						% Estimation indicator of the skew coefficient
est_dist = [ 1 ; 1 ; 1 ; 1 ; 0 ];	% Estimation indicator of the distortion coefficients


%-- Extrinsic parameters:
%-- The rotation (omc_kk) and the translation (Tc_kk) vectors for every calibration image and their uncertainties

%-- Image #1:
omc_1 = [ 1.612388e+00 ; 1.637573e+00 ; -9.328938e-01 ];
Tc_1  = [ -1.823517e+02 ; 2.621429e+02 ; 1.158292e+03 ];
omc_error_1 = [ 5.383572e-02 ; 5.934929e-02 ; 7.432466e-02 ];
Tc_error_1  = [ 2.310921e+01 ; 1.621578e+02 ; 5.137561e+01 ];

%-- Image #2:
omc_2 = [ 1.695513e+00 ; 1.660658e+00 ; -9.019329e-01 ];
Tc_2  = [ -1.623654e+02 ; 1.536715e+02 ; 1.088367e+03 ];
omc_error_2 = [ 4.751467e-02 ; 5.681360e-02 ; 7.527002e-02 ];
Tc_error_2  = [ 2.167834e+01 ; 1.464753e+02 ; 4.934743e+01 ];

%-- Image #3:
omc_3 = [ 1.762459e+00 ; 1.838526e+00 ; -7.867508e-01 ];
Tc_3  = [ -1.629349e+02 ; 1.950151e+02 ; 1.120423e+03 ];
omc_error_3 = [ 3.288593e-02 ; 4.442837e-02 ; 7.485634e-02 ];
Tc_error_3  = [ 2.234516e+01 ; 1.529854e+02 ; 4.999031e+01 ];

%-- Image #4:
omc_4 = [ 1.343222e+00 ; 2.239050e+00 ; -1.131323e+00 ];
Tc_4  = [ -7.932753e+01 ; 1.325341e+02 ; 1.252563e+03 ];
omc_error_4 = [ 1.881800e-02 ; 6.091463e-02 ; 9.805843e-02 ];
Tc_error_4  = [ 2.498658e+01 ; 1.650697e+02 ; 5.417052e+01 ];

