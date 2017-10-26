% This File should save relevant parameters to a text file
%
cam_width = nx_cam;
cam_height = ny_cam;
cam_intrinsic = [fc_cam(1)  0           cc_cam(1);  ...
                 0          fc_cam(2)   cc_cam(2);  ...
                 0          0           1];
cam_distortion = [0 0 0 0];
cam_rotation = [1 0 0; 0 1 0; 0 0 1;];
cam_translation = [0 0 0];             

fid = fopen('CamCalib.txt','w');
fprintf(fid,'%i\n', cam_width);
fprintf(fid,'%i\n',             cam_height);
fprintf(fid,'%f %f %f \n',      cam_intrinsic(1,:));
fprintf(fid,'%f %f %f \n',      cam_intrinsic(2,:));
fprintf(fid,'%f %f %f \n',      cam_intrinsic(3,:));
fprintf(fid,'%f %f %f %f \n',   cam_distortion);
fprintf(fid,'%f %f %f \n',      cam_rotation(1,:));
fprintf(fid,'%f %f %f \n',      cam_rotation(2,:));
fprintf(fid,'%f %f %f \n',      cam_rotation(3,:));
fprintf(fid,'%f %f %f',         cam_translation);
fclose('all');



proj_width = nx_proj;
proj_height = ny_proj;
proj_intrinsic = [fc_proj(1)  0           cc_proj(1);  ...
                 0          fc_proj(2)   cc_proj(2);  ...
                 0          0                    1];
proj_distortion = [0 0 0 0];
proj_rotation = R_proj;
proj_translation = T_proj;

fid = fopen('ProjCalib.txt','w');
fprintf(fid,'%i\n', proj_width);
fprintf(fid,'%i\n',             proj_height);
fprintf(fid,'%f %f %f \n',      proj_intrinsic(1,:));
fprintf(fid,'%f %f %f \n',      proj_intrinsic(2,:));
fprintf(fid,'%f %f %f \n',      proj_intrinsic(3,:));
fprintf(fid,'%f %f %f %f \n',   proj_distortion);
fprintf(fid,'%f %f %f \n',      proj_rotation(1,:));
fprintf(fid,'%f %f %f \n',      proj_rotation(2,:));
fprintf(fid,'%f %f %f \n',      proj_rotation(3,:));
fprintf(fid,'%f %f %f',         proj_translation);
fclose('all');
