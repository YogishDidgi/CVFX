close all;
clear all;
bgfn=400;
for fgfn=60:285
    if fgfn<100
    fgname=strcat('D:\Dropbox\CVFX\Project\Data\Frames_ffmpeg\MVI_0077\00',num2str(fgfn),'.jpg');
    else
      fgname=strcat('D:\Dropbox\CVFX\Project\Data\Frames_ffmpeg\MVI_0077\0',num2str(fgfn),'.jpg'); 
    end
im=imread(fgname);
im=imtranslate(im,[-250,100],'fillvalue',[44,145,51]);
imshow(im)

%im=imresize(im,0.1);
bgname=strcat('D:\Dropbox\CVFX\Project\Data\Frames_ffmpeg\hand_out_screen\0',num2str(bgfn),'.jpg'); 
bg1=imread(bgname);
%bg1=imresize(bg1,[size(im,1) size(im,2)]);
%%
alpha=vlahos_greenscreen(im,5,1);
im=double(im)/255;
bg1=double(bg1)/255;
res(:,:,1)=alpha.*(im(:,:,1))+(1-alpha).*(bg1(:,:,1));
res(:,:,2)=alpha.*(im(:,:,2))+(1-alpha).*(bg1(:,:,2));
res(:,:,3)=alpha.*(im(:,:,3))+(1-alpha).*(bg1(:,:,3));
finalname=strcat('D:\Dropbox\CVFX\Project\Results\trial results\creepy_hand2\0',num2str(bgfn),'.jpg'); 
imwrite(res,finalname);
bgfn=bgfn+1;
end

%%
% removing yogish from behind
imwoyogish=imread('D:\Dropbox\CVFX\Project\Data\Frames_ffmpeg\tshirt_g_hand_out_screen\0001.jpg');
imtool(imwoyogish)
for bgfn=400:625

fname=strcat('D:\Dropbox\CVFX\Project\Results\trial results\creepy_hand_ref\0',num2str(bgfn),'.jpg');
final_ref=imread(fname);
imtool(final_ref)
final_ref(1:500,1:400,:)=imwoyogish(1:500,1:400,:);
imtool(final_ref)

finalname=strcat('D:\Dropbox\CVFX\Project\Results\trial results\creepy_hand_ref2\0',num2str(bgfn),'.jpg'); 
imwrite(final_ref,finalname);
end

