function spin_imgs = spin_images(model, varargin)
%SPIN_IMAGES   Calculate 3D spin image descriptors.
%
% spin_imgs = SPIN_IMAGES(model) calculates spin image descriptors for all
% the points in model with default parameters.
%
% Input:
%   model - The n-by-3 input point set
%   radius - The spatial radius for controlling the spin image locality
%   imgw - The 2D histogram width (descriptor dimension is imgW^2)
%   minn - Minimum number of neighbors required for calculating a
%          valid descriptor
% 
% Output:
%   spin_imgs - The n-by-imgW^2 output descriptor set
%
% Johnson, A.E.; Hebert, M.; , "Using spin images for efficient object
% recognition in cluttered 3D scenes," Pattern Analysis and Machine
% Intelligence, IEEE Transactions on , vol.21, no.5, pp.433-449, May 1999

% Defaults
radius = 0.025;
imgw = 10;
minn = 10;
if nargin > 1
    radius = varargin{1};
    if nargin > 2
        imgw = varargin{2};
        if nargin > 3
            minn = varargin{3};
        end
    end
end

imgh = imgw/2;
N = size(model,1);

spin_imgs_tmp = zeros(2*imgh,imgw,N);
neighborRadius = radius * sqrt(2); % neighbor search radius has to be larger because the spin image is cylindrical

[idx dist] = rangesearch(model, model, neighborRadius - 1e-7);

% loop over points
parfor i=1:N
   %[ i N]
   pt = model(i,:);
   %pt = model(:,i)';
   spinImg = spin_imgs_tmp(:,:,i);
   
   % Neighbors of current point
    neighbors = model(idx{i},:);

   if size(neighbors,1) >= minn
      % first we compute the normal vector
      % compute PCA
      coeff = princomp(neighbors);
      
      % third component is surface normal
      normal = coeff(:,3);
      if dot(normal, pt) < 0
         normal = -normal;
      end

      % now we compute the spin image
      nn = length(neighbors);
      diffs = (neighbors - repmat(pt,nn,1))./radius;
      lens = sqrt(dot(diffs',diffs'))';

      % y-coord is dot prod between normal and vector to neighbor
      yvals = dot(repmat(normal',nn,1),diffs,2);
      
      % x-coord is distance of the neighbor to the normal line
      xvals = sqrt(lens.^2 - yvals.^2);

      % only add points if they are actually in the spin image
      xyvals = [xvals yvals];
      xyvals = xyvals(abs(xvals) < 1 & abs(yvals) < 1,:);
      nPts = length(xyvals);
      xvals = xyvals(:,1); yvals = xyvals(:,2);
      xInds = round(xvals.*(imgw-1))+1;
      yInds = imgh + round(yvals.*(imgh-1))+1;

      if(nnz(xInds < imgw & yInds < 2*imgh) == 0)
         pt = pt
         disp('error computing spin image!');
      else
         inds = sub2ind(size(spinImg),xInds,yInds);
         for j=1:length(inds)
            spinImg(inds(j)) = spinImg(inds(j)) + 1;
         end
      end
      
      % if only 1 bin is occupied, then this corresponds to the query point
      % --> neighborhood is too sparse, i.e. no spin image
      if(nnz(spinImg >= 1) == 1)
         spinImg = zeros(2*imgh,imgw);
      else
         % normalize spin image
         spinImg = spinImg./nPts;
      end
      spin_imgs_tmp(:,:,i) = spinImg;
   else
      disp('not enough neighbors!');
   end
end

bins = imgw*imgw;
spin_imgs = zeros(N, bins);
for i=1:N
    spin_imgs(i,:) = reshape(spin_imgs_tmp(:,:,i), 1, bins);
end

end