% cams=readPJ(filename) (version 0.1)
% loads the camera calibration parameters from an openCV like calibration
% file
% 
% INPUT:
% 
%     * filename: name (and path) to the file to be read
% 
% OUTPUT:
% 
%     * cams: a array which contains for every camera a struct with the
%       fields:
%       width           (1x1)       image width
%       height          (1x1)       image height
%       intrinsic       (3x3)       instrinsic parameters
%       distortion      (4x1)       distortion parameters, 2 rad., 2 tang.
%       rotation        (3,3)       external parameters, rotation matrix
%       translation     (3x1)       external parameters, translation vector
%       transformation  (4x4)       ext. par. as am homogenous transform
%       P               (3x4)       complete projection matrix
% 
% NOTES:
% 
%     * -This function is a prototype. If you find problems please let me
%        know.
% 
% Author: Dirk Kraft (kraft@mmmi.sdu.dk)
%
% Date: 17.03.2010
%
function [cams] = readPJ(filename)
    calFile = fopen (filename);
    numCameras = fscanf (calFile,'%d',1);

    for j=1:numCameras,
        width  = fscanf (calFile,'%d',1);
        height = fscanf (calFile,'%d',1);
        intrinsic = reshape(fscanf (calFile,'%g',9),3,3)';
        distortion = fscanf (calFile,'%g',4)';
        rotation = reshape(fscanf (calFile,'%g',9),3,3)';
        translation = fscanf (calFile,'%g',3);
        transformation = [rotation translation ; 0 0 0 1];
        P = [intrinsic, [0 0 0]'] * transformation;
        cams(j) = struct ('width',width, ...
                          'height',height, ...
                          'intrinsic',intrinsic, ...
                          'distortion',distortion, ...
                          'rotation',rotation, ...
                          'translation',translation, ...
                          'transformation',transformation, ...
                          'P',P);
    end
    fclose(calFile);