function [ alpha ] = vlahos_greenscreen( img, a1,a2 )
% this function returns the alpha value for green screen matting
img=double(img)/255;
alpha=1-a1*(img(:,:,2)-a2*img(:,:,1));
alpha=imadjust(alpha,stretchlim(alpha));
%imshow(alpha,[0 1]);
end

