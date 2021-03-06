%Computing Exterior projection after calibration

%LEFT SIDE
KA_left = [fc_left(1) 0 cc_left(1) 0; 0 fc_left(2) cc_left(2) 0; 0 0 1 0; 0 0 0 0]
H_left = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]

KA_left*H_left

%RIGHT SIDE
R = rodrigues(om);
KA_right = [fc_right(1) 0 cc_right(1) 0; 0 fc_right(2) cc_right(2) 0; 0 0 1 0; 0 0 0 0]
H_right = [R(1,1) R(1,2) R(1,3) T(1); R(2,1) R(2,2) R(2,3) T(2); R(3,1) R(3,2) R(3,3) T(3); 0 0 0 1]

KA_right*H_right
