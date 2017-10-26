% Clean up
clear;clc;close all;

% Bring 'ext' folder into scope
path('ext', path);

% Read the object and the scene
object = dlmread('rhino.txt');
scene = dlmread('scene3.txt');

% Show the object and the scene
show_obj_scn(object, scene, 'Input object and scene');

% When done with pose estimation, store the aligned model here
object_align = object;

%%%%%%%%%%%%%%%%%%%% YOUR RANSAC CODE GOES HERE %%%%%%%%%%%%%%%%%%%%
% Recommendation: use the spin image descriptors to find correspondences
% between object/scene points. Three point pairs are enough for estimating
% a transformation in each RANSAC  iteration. Then you can evaluate a fit
% error or count the number of inliers to check the sampled transformation.
%
% Functions that you may want to use:
%   spin_images - Find spin image descriptors for a point set
%   randperm - Find random unique indices using sampling without
%              replacement
%   knnsearch - Find nearest neighbor(s) between sets of vectors, both 3D
%               points and descriptors
%   transformation_estimation_svd - Estimate transformation between two 3D
%                                   point sets
%   transform_points - Apply a transformation to a set of points
%
% NOTE: when you are done, remember to apply the resulting transformation
% to object_align for the alignment visualization below to work.

a = 2;

%spin_imgs_object = spin_images(object);
%spin_imgs_sceen  = spin_images(scene);
a


%    IDX = KNNSEARCH(object_align_2,scene2);
%    temp = zeros(length(IDX),3);
%    for i = 1:length(IDX)
%        temp(i,:) = object_align_2(IDX(i),:);
%    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Show RANSAC alignment
show_obj_scn(object_align, scene, 'RANSAC result');