close all;
clear all;
% bgfn=400;
for fgfn=1:11

    fgname=strcat('D:\Dropbox\CVFX\Project\Code\Prob3\FaceMorph\yog_ind_green\Name_',num2str(fgfn),'.jpg');

im=imread(fgname);
im1=imresize(im,[750 940]);
% imshow(im1)
im2=padarray(im1,[328 1],'replicate','pre');
im2=padarray(im2,[1 490],'replicate','pre');
im2=padarray(im2,[1 489],'replicate','post');
% imshow(im2);
bgname=strcat('D:\Dropbox\CVFX\Project\Data\Frames_ffmpeg\MVI_0113\0001.jpg'); 
bg1=imread(bgname);
im=im2;

%%
alpha=vlahos_greenscreen(im,10,1);

im=double(im)/255;
bg1=double(bg1)/255;
res(:,:,1)=alpha.*(im(:,:,1))+(1-alpha).*(bg1(:,:,1));
res(:,:,2)=alpha.*(im(:,:,2))+(1-alpha).*(bg1(:,:,2));
res(:,:,3)=alpha.*(im(:,:,3))+(1-alpha).*(bg1(:,:,3));

imshow(res);

finalname=strcat('D:\Dropbox\CVFX\Project\Code\Prob3\FaceMorph\resized yog ind\comp',num2str(fgfn),'.jpg'); 
imwrite(res,finalname);

end
%% making yogish's shirt opaque in resized images
for fgfn=1:11

    fgname=strcat('D:\Dropbox\CVFX\Project\Code\Prob3\FaceMorph\resized composite\resized yog gs\comp',num2str(fgfn),'.jpg');
im=imread(fgname);
    greenimgname=strcat('D:\Dropbox\CVFX\Project\Code\Prob3\FaceMorph\gs_yog_green\Name_',num2str(fgfn),'.jpg');

greenimg=imread(greenimgname);
im1=imresize(greenimg,[750 940]);
% imshow(im1)
im2=padarray(im1,[328 1],'replicate','pre');
im2=padarray(im2,[1 490],'replicate','pre');
im2=padarray(im2,[1 489],'replicate','post');
greenimg=im2;

imtool(im);
im(794:1050,800:1098,:)=greenimg(794:1050,800:1098,:);
im(1:527,1296:1920,:)=bg1(1:527,1296:1920,:);

imtool(im);

finalname=strcat('D:\Dropbox\CVFX\Project\Code\Prob3\FaceMorph\resized composite\resized yog gs\comp',num2str(fgfn),'.jpg'); 
imwrite(im,finalname);
end

