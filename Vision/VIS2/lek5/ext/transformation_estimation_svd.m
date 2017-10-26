function T = transformation_estimation_svd(P, Q)
%TRANSFORMATION_ESTIMATION_SVD   Estimate 4-by-4 homogeneous 
% transformation matrix.
% 
% T = TRANSFORMATION_ESTIMATION_SVD(P,Q) estimates a transformation matrix 
% from corresponding 3D point sets. At least 3 point correspondences must 
% be available.
%
% Input:
%   P - First n-by-3 point set
%   Q - Second n-by-3 point set
%
% Output:
%   T - Estimated transformation between the point sets
%
% This is an adaption of the algorithm found in PCL's
% TransformationEstimationSVD:
% http://docs.pointclouds.org/1.5.1/classpcl_1_1registration_1_1_transformation_estimation_s_v_d.html
%
% Author: Anders Glent Buch


% Initialize result
T = eye(4);

% Get number of points
n = size(P,1);

% Check number of points
if n < 3
    error('Too few points (< 3) for transformation estimation!')
end

% Check dimensions
if size(P,2) ~= 3 || size(Q,2) ~= 3 ||  size(Q,1) ~= n
    error('Inconsistent sizes for transformation estimation!');
end

% Now both P and Q should have size n-by-3, n >= 3

% Demean all points
meanP = mean(P);
meanQ = mean(Q);
PP = P-repmat(meanP,n,1);
QQ = Q-repmat(meanQ,n,1);

% Correlation matrix, 3-by-3
H = PP' * QQ;

% Compute SVD
[U S V] = svd(H);

% Compute R
if det(U)*det(V) < 0
    V(:,3) = -1;
end
R = V * U';

% Compute t
Rc = R*meanP';
t = meanQ' - Rc;

% Store
T(1:3,1:3) = R;
T(1:3,4) = t;
end