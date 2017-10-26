function Q = transform_points(P, T)
%TRANSFORM_POINTS   Transform a set of 3D points.
%
% Q = TRANSFORM_POINTS(P, T) applies a 4-by-4 homogeneous transformation to
% a set of 3D points.
%
% Input:
%   P - The n-by-3 input point set
%   T - A homogeneous 4-by-4 transformation matrix
% 
% Output:
%   Q - The n-by-3 output point set of transformed points
%
% Author: Anders Glent Buch


[r c] = size(P);
Q = zeros(r, c);

R = T(1:3,1:3);
t = T(1:3,4);

for i=1:r
    Q(i,:) = (R*P(i,:)'+t)';
end
end