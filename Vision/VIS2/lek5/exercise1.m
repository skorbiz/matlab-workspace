% Clean up
clear;clc;close all;

% Bring 'ext' folder into scope
path('ext', path);

% Read the object and the scenes
object = dlmread('rhino.txt');
scene1 = dlmread('scene1.txt');
scene2 = dlmread('scene2.txt');
scene3 = dlmread('scene3.txt');

% Show the object and each scene
figure
show_obj_scn(object, scene1, 'Input object and scene 1', subplot(131));
show_obj_scn(object, scene2, 'Input object and scene 2', subplot(132));
show_obj_scn(object, scene3, 'Input object and scene 3', subplot(133));

% When done with pose estimation, store the aligned models here
object_align_1 = object;
object_align_2 = object;
object_align_3 = object;

%%%%%%%%%%%%%%%%%%%% YOUR ICP CODE GOES HERE %%%%%%%%%%%%%%%%%%%%
% Functions that you may want to use:
%   knnsearch - Find nearest neighbor(s) between sets of vectors, both 3D
%               points and descriptors
%   transformation_estimation_svd - Estimate transformation between two 3D
%                                   point sets
%   transform_points - Apply a transformation to a set of points
%
% NOTE: when you are done, remember to apply the resulting transformation
% to object_align* for the alignment visualizations below to work.

itterations = 15;
object_align_1 = object;
for j = 1:itterations
    IDX = KNNSEARCH(object_align_1,scene1);
    temp = zeros(length(IDX),3);
    for i = 1:length(IDX)
        temp(i,:) = object_align_1(IDX(i),:);
    end
    object_align_1 = temp;
    T = transformation_estimation_svd(object_align_1,scene1);
    object_align_1 = transform_points(object_align_1,T);
end

object_align_2 = object;
for j = 1:itterations
    IDX = KNNSEARCH(object_align_2,scene2);
    temp = zeros(length(IDX),3);
    for i = 1:length(IDX)
        temp(i,:) = object_align_2(IDX(i),:);
    end
    object_align_2 = temp;
    T = transformation_estimation_svd(object_align_2,scene2);
    object_align_2 = transform_points(object_align_2,T);
end


object_align_3 = object;
for j = 1:itterations
    IDX = KNNSEARCH(object_align_3,scene3);
    temp = zeros(length(IDX),3);
    for i = 1:length(IDX)
        temp(i,:) = object_align_3(IDX(i),:);
    end
    object_align_3 = temp;
    T = transformation_estimation_svd(object_align_3,scene3);
    object_align_3 = transform_points(object_align_3,T);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Show ICP alignments
figure
show_obj_scn(object_align_1, scene1, 'ICP result scene 1', subplot(131));
show_obj_scn(object_align_2, scene2, 'ICP result scene 2', subplot(132));
show_obj_scn(object_align_3, scene3, 'ICP result scene 3', subplot(133));

