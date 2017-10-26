% ******************************************************
%   Course:               VIS02
%   Author:               Dirk Kraft
% ******************************************************

FILENAME_IMG_LEFT = 'left.png';
FILENAME_IMG_RIGHT = 'right.png';
FILENAME_CALIBRATION = 'calibration.txt';

%********************************************************************
% 1) LOADING THE STEREO CALIBRATION MATRIX
%********************************************************************

s = readPJ(FILENAME_CALIBRATION);

% intrinsic
KAl =   s(1).intrinsic;

KAr =   s(2).intrinsic;

% extrinsic
Hl =   s(1).transformation;
Hr =   s(2).transformation;
%********************************************************************
% 2) COMPUTE P
%********************************************************************

Pl = [ KAl [0 0 0]'] * Hl;
Pr = [ KAr [0 0 0]'] * Hr;


%********************************************************************
% 3) COMPUTE THE OPTICAL CENTER FOR LEFT AND RIGHT ( C=-pow(-R,-1)*t )
%********************************************************************

PXl = Pl(1:3,1:3);
PXr = Pr(1:3,1:3);

pxl = Pl(1:3,4);
pxr = Pr(1:3,4);

Cl =  [(-inv(PXl) * pxl); 1];
Cr =  [(-inv(PXr) * pxr); 1];

%********************************************************************
% 4) COMPUTE EPIPOLES (e = P * C)
%********************************************************************

el = Pl * Cr;
er = Pr * Cl;

%********************************************************************
% 5) COMPUTE FUNDAMENTAL MATRIX F12
%********************************************************************

erx = [ 0 -er(3) er(2);
    er(3) 0 -er(1);
    -er(2) er(1) 0];

Flr = erx * Pr * pinv(Pl);

Flr

%********************************************************************
% 6) LOADING & DISPLAYING IMAGES
%********************************************************************

IMG_LEFT =imread(FILENAME_IMG_LEFT);
IMG_RIGHT = imread(FILENAME_IMG_RIGHT);

%********************************************************************
% 7) DETECT CLICKED IMAGE COORDINATE IN LEFT IMAGE
%********************************************************************
figure;

imshow(IMG_LEFT);
hold on
[click_left_x,click_left_y] = getpts(); % pixel in hom coordinates

for i = 1:size(click_left_x)
    plot([click_left_x(i)-3 click_left_x(i)+3],[click_left_y(i) click_left_y(i)],'Color','r','LineWidth', 1)
    plot([click_left_x(i) click_left_x(i)],[click_left_y(i)-3 click_left_y(i)+3],'Color','r','LineWidth', 1)
    
    plot([click_left_x(i)-3 click_left_x(i)+3],[click_left_y(i)-3 click_left_y(i)+3],'Color','b','LineWidth', 1)
    plot([click_left_x(i)-3 click_left_x(i)+3],[click_left_y(i)+3 click_left_y(i)-3],'Color','b','LineWidth', 1)
end;



hold off

figure
for i = 1:size(click_left_x)
    %********************************************************************
    % 6) DISPLAY EPIPOLAR LINES IN RIGHT IMAGE
    %********************************************************************
    hold off
    imshow(IMG_RIGHT);
    hold on
    
    click_left(1) = click_left_x(i);
    click_left(2) = click_left_y(i);
    click_left(3) = 1;
    
    click_left(1) = 500;
    click_left(2) = 500;

    
    Minf = [inv(PXl)*click_left'; 0]; % computing point at infinity
    mr = Pr*Minf; % reprojection of point at infinity in right image
    
    %********************************************************************
    % 5) COMPUTE EPIPOLAR LINES (l = e_2 x m_1, inf)
    %********************************************************************
    
    if (er(3) ~= 0.0)
        er_norm = er/er(3);
    else
        er_norm = er/10e-50;
    end
    
    % er_norm = er
    
    mr_norm = mr/mr(3);
    %********************************************************************
    % 6) DRAW EPIPOLAR LINES
    %********************************************************************
    
    image_size = size(IMG_RIGHT);
    height = image_size(1);
    width = image_size(2);
    
    %plot([er_norm(1) mr_norm(1)],[er_norm(2) mr_norm(2)],'--','Color','g','LineWidth', 3);
    
    plot([-14344 718],[158 504],'--','Color','g','LineWidth', 3);
    
    line = Flr * click_left';
    pointL = [0 (-1.0 * line(3)/line(2))];
    pointR = [width (-1.0 * line(1)/line(2) * width - line(3)/line(2))];
    plot([pointL(1) pointR(1)],[pointL(2) pointR(2)],'--','Color','r','LineWidth', 1);
    
    line2 = cross(er,mr);
    point2L = [0 (-1.0 * line2(3)/line2(2))];
    point2R = [width (-1.0 * line2(1)/line2(2) * width - line2(3)/line2(2))];
    plot([point2L(1) point2R(1)],[point2L(2) point2R(2)],'-.','Color','b','LineWidth', 1);
    
    
    % 7) DETECT CLICKED IMAGE COORDINATE IN RIGHT IMAGE
    hold off
    hold on
    [click_right_x,click_right_y] = getpts(); % pixel in hom coordinates
    hold off
    % 8) CREATE 3D LINES FROM CLICKED 2D POINTS
    click_right(1) = click_right_x(size(click_right_x,1));
    click_right(2) = click_right_y(size(click_right_x,1));
    click_right(3) = 1;
    
    M_inf = inv(PXr)*click_right'; % computing point at infinity
    
    Minf = Minf(1:3);
    M_inf = M_inf(1:3);
    
    mu1 = cross(Cl(1:3),Minf) / norm(Minf);
    v1 = Minf / norm(Minf);
    
    mu2 = cross(Cr(1:3),M_inf) / norm(M_inf);
    v2 = M_inf / norm(M_inf);
    
    % 9) COMPUTE INTERSECTION OF THOSE LINES
    
    M1 = ((v1 * cross(v2,mu2)' - (v1 * v2') * v1 * cross(v2,mu1)') / (norm(cross(v1,v2))^2)) * v1 + cross(v1,mu1);
    M2 = ((v2 * cross(v1,mu1)' - (v2 * v1') * v2 * cross(v1,mu2)') / (norm(cross(v2,v1))^2)) * v2 + cross(v2,mu2);
    
    avgM = M1 + (M2-M1)/2;
    
    
    % 10) PRINT RESULT
    
    display('******************************')
    i
    click_left
    click_right
    M1
    M2
    avgM
    display('******************************')
end;

hold off


