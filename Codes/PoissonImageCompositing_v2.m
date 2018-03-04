function PoissonImageCompositing_v2
startFrameFolder1 = 350;
startFrameFolder2 = 258;
Folder1 = 'Input1\';
Folder2 = 'Input2\';
ResultFolder = 'F:\Workspace\Matlab\CVFX_Fall2015\Project\Output\';

Stage1 = true;
Stage2 = false;
Stage3 = false;
Stage4 = false;

if(Stage1)
    startIter = 1;
    endIter = 198;
    numImages = endIter - startIter + 1;
    for i = startIter:endIter
        fprintf('Processing Stage 1: image %d out of %d\t',i,numImages);
        tic
        iterSource = startFrameFolder1 + i - 1;
        iterTarget = startFrameFolder2 + i - 1;
        imageSourceOrig = imread(strcat(Folder1,num2str(iterSource,'%04i'),'.jpg'),'jpg');
        imageTargetOrig = imread(strcat(Folder2,num2str(iterTarget,'%04i'),'.jpg'),'jpg');
        imageFinalName = [ResultFolder num2str(i,'%04i') '.jpg'];
        [imageSource,imageTarget,imageMask] = getImages(imageSourceOrig,imageTargetOrig,i,iterSource);
%         figure,
        subplot(2,2,1),imshow(imageSource);
        subplot(2,2,2),imshow(imageTarget);
        subplot(2,2,3),imshow(imageMask);

        finalImage = Poisson(imageSource,imageTarget,imageMask);
%         imwrite(finalImage, imageFinalName);
        subplot(2,2,4),imshow(finalImage);
        t1 = toc;
        fprintf('Took %f sec\n',t1);
    end
end

if(Stage2)
    startIter = 199;
    endIter = 264;
    numImages = endIter - startIter + 1;
    for i = startIter:endIter
        fprintf('Processing Stage 2-1(Poisson): image %d out of %d\t',i - startIter + 1,numImages);
        tic
        iterSource = startFrameFolder1 + i - 1;
        iterTarget = startFrameFolder2 + i - 1;
        imageSourceOrig = imread(strcat(Folder1,num2str(iterSource,'%04i'),'.jpg'),'jpg');
        imageTargetOrig = imread(strcat(Folder2,num2str(iterTarget,'%04i'),'.jpg'),'jpg');
        imageFinalName = [ResultFolder num2str(i,'%04i') '.jpg'];
        [imageSource,imageTarget,imageMask] = getImages(imageSourceOrig,imageTargetOrig,i,iterSource);
    %     figure,
    %     subplot(2,2,1),imshow(imageSource);
    %     subplot(2,2,2),imshow(imageTarget);
    %     subplot(2,2,3),imshow(imageMask);
    %     
        finalImage = Poisson(imageSource,imageTarget,imageMask);

        imwrite(finalImage, imageFinalName);
    %     subplot(2,2,4),imshow(finalImage);
        t1 = toc;
        fprintf('Took %f sec\n',t1);
    end
end

if(Stage3)
    startIter = 265;
    endIter = 379;
    numImages = endIter - startIter + 1;
    for i = startIter:endIter
        fprintf('Processing Stage 2-2(Hard composite): image %d out of %d\t',i - startIter + 1,numImages);
        tic
        iterSource = startFrameFolder1 + i - 1;
        iterTarget = startFrameFolder2 + i - 1;
        imageSourceOrig = imread(strcat(Folder1,num2str(iterSource,'%04i'),'.jpg'),'jpg');
        imageTargetOrig = imread(strcat(Folder2,num2str(iterTarget,'%04i'),'.jpg'),'jpg');
        imageFinalName = [ResultFolder num2str(i,'%04i') '.jpg'];
        [imageSource,imageTarget,imageMask] = getImages(imageSourceOrig,imageTargetOrig,i,iterSource);
    %     figure,
    %     subplot(2,2,1),imshow(imageSource);
    %     subplot(2,2,2),imshow(imageTarget);
    %     subplot(2,2,3),imshow(imageMask);
    %     
        imageSourceDouble = double(imageSource);
        imageTargetDouble = double(imageTarget);
        imageMaskDouble = double(imageMask/255);
        imageMask3Double = zeros(size(imageMaskDouble, 1), size(imageMaskDouble, 2), 3);
        for i = 1:3
            imageMask3Double(:,:,i) = imageMaskDouble;
        end
        imageMask3Double = imgaussfilt(imageMask3Double,4);
        imageHardComposite = imageMask3Double.*imageSourceDouble + (1 - imageMask3Double).*imageTargetDouble;
        finalImage = uint8(imageHardComposite);

        imwrite(finalImage, imageFinalName);
    %     subplot(2,2,4),imshow(finalImage);
        t1 = toc;
        fprintf('Took %f sec\n',t1);
    end
end

if(Stage4)
    startIter = 380;
    endIter = 551;
    numImages = endIter - startIter + 1;
    for i = startIter:endIter
        fprintf('Processing Stage 3: image %d out of %d\t',i - startIter + 1,numImages);
        tic
        iterSource = startFrameFolder1 + i - 1;
        iterTarget = startFrameFolder2 + i - 1;
        imageSourceOrig = imread(strcat(Folder1,num2str(iterSource,'%04i'),'.jpg'),'jpg');
        imageTargetOrig = imread(strcat(Folder2,num2str(iterTarget,'%04i'),'.jpg'),'jpg');
        imageFinalName = [ResultFolder num2str(i,'%04i') '.jpg'];
        [imageSource,imageTarget,imageMask] = getImages(imageSourceOrig,imageTargetOrig,i,iterSource);
    %     figure,
    %     subplot(2,2,1),imshow(imageSource);
    %     subplot(2,2,2),imshow(imageTarget);
    %     subplot(2,2,3),imshow(imageMask);
    %     
        finalImage = Poisson(imageSource,imageTarget,imageMask);
        imwrite(finalImage, imageFinalName);
    %     subplot(2,2,4),imshow(finalImage);
        t1 = toc;
        fprintf('Took %f sec\n',t1);
    end
end
end

function finalImage = Poisson(imageSource,imageTarget,imageMask)
    imageSourceDouble = double(imageSource);
    imageTargetDouble = double(imageTarget);
    imageMaskDouble = double(imageMask/255);
    imageMaskDoubleCopy = imageMaskDouble;
    imageMaskDoubleCopy(find(imageMaskDoubleCopy)) = 1:sum(imageMaskDoubleCopy(:));
    laplacianMask = delsq(double(imageMaskDoubleCopy));
    imagePoissonComposite = imageTargetDouble;
    [m,n] = size(imageMaskDouble);
    
%     gradXFilter = [0 -1 1];
%     gradYFilter = [0 -1 1]';
%     gradImageXSource = imfilter(imageSourceDouble, gradXFilter);
%     gradImageYSource = imfilter(imageSourceDouble, gradYFilter);
%     gradImageXTarget = imfilter(imageTargetDouble, gradXFilter);
%     gradImageYTarget = imfilter(imageTargetDouble, gradYFilter);
%     normCmpMask = double(((gradImageXTarget.*gradImageXTarget + gradImageYTarget.*gradImageYTarget) > ...
%                     (gradImageXSource.*gradImageXSource + gradImageYSource.*gradImageYSource)));
%     imageMask3Double = zeros(size(imageMaskDouble, 1), size(imageMaskDouble, 2), 3);
%     for i = 1:3
%         imageMask3Double(:,:,i) = imageMaskDouble;
%     end
%     guidanceX = imageMask3Double.*normCmpMask;
%     guidanceY = imageMask3Double.*normCmpMask;
%     gradImageXFinal = guidanceX.*gradImageXTarget + (1 - guidanceX).*gradImageXSource;
%     gradImageYFinal = guidanceY.*gradImageYTarget + (1 - guidanceY).*gradImageYSource;
%     laplaceImageXFinal = imfilter(gradImageXFinal, gradXFilter);
%     laplaceImageYFinal = imfilter(gradImageYFinal, gradYFilter);
    
    %Poisson
    for channelID = 1:3
        RHS = double(zeros(size(laplacianMask,1), 1));
        iter = 1;
        for j = 1:n
            for i = 1:m
%                 if(((i - 1) > 0) && ((i + 1) <= size(imageMaskDouble, 1)) && ((j - 1) > 0) && ((j + 1) <= size(imageMaskDouble, 2)))
                    if(imageMaskDouble(i, j)==1)
                        RHS(iter) = (4*imageSourceDouble(i, j, channelID) - imageSourceDouble(i - 1, j, channelID) ...
                                    - imageSourceDouble(i + 1, j, channelID) - imageSourceDouble(i, j - 1, channelID) ...
                                    - imageSourceDouble(i, j + 1, channelID)) ...
                                    + (1.0 - imageMaskDouble(i - 1, j))*imageTargetDouble(i - 1, j, channelID) ...
                                    + (1.0 - imageMaskDouble(i + 1, j))*imageTargetDouble(i + 1, j, channelID) ...
                                    + (1.0 - imageMaskDouble(i, j - 1))*imageTargetDouble(i, j - 1, channelID) ...
                                    + (1.0 - imageMaskDouble(i, j + 1))*imageTargetDouble(i, j + 1, channelID);
%                         RHS(iter) = -(laplaceImageXFinal(i, j, channelID) + laplaceImageYFinal(i, j, channelID)) ...
%                                     + (1.0 - imageMaskDouble(i - 1, j))*imageTargetDouble(i - 1, j, channelID) ...
%                                     + (1.0 - imageMaskDouble(i + 1, j))*imageTargetDouble(i + 1, j, channelID) ...
%                                     + (1.0 - imageMaskDouble(i, j - 1))*imageTargetDouble(i, j - 1, channelID) ...
%                                     + (1.0 - imageMaskDouble(i, j + 1))*imageTargetDouble(i, j + 1, channelID);
                        iter = iter + 1;
                    end
%                 end
            end
        end
        LHS = laplacianMask\RHS;
        iter = 1;
        for j = 1:n
            for i = 1:m
                if(imageMaskDouble(i, j)==1)
                    imagePoissonComposite(i, j, channelID) = LHS(iter);
                    iter = iter + 1;
                end
            end
        end    
    end
    finalImage = uint8(imagePoissonComposite);
end

function [imageSource,imageTarget,imageMask] = getImages(imageSourceOrig,imageTargetOrig,currentIndex,iterSource)
imageSource = padarray(imageSourceOrig,[1 1],'replicate','both');
imageTarget = padarray(imageTargetOrig,[1 1],'replicate','both');
imageMask = zeros(size(imageTarget,1),size(imageTarget,2));

if(currentIndex <= 198)
    imageMask(460:(end-1),180:570) = 255;
elseif(currentIndex>=380 && currentIndex<=551)
    imageMask(400:(end-1),130:600) = 255;
    imageTemp = imageSource;
    imageSource = imageTarget;
    imageTarget = imageTemp;
elseif(currentIndex>=199 && currentIndex<=264)
%     imageMask3 = imread(['Masks/' num2str(iterSource,'%04i') '.jpg']);
    imageMask3 = imread(['Masks/0603.jpg']);
    imageMask = imageMask3(:,:,1);
    imageMask = padarray(imageMask,[1 1],'replicate','both');
%     imageMask = imresize(imageMask, [size(imageTarget, 1) size(imageTarget, 2)]);
else
    imageMask3 = imread(['Masks/' num2str(iterSource,'%04i') '.jpg']);
    imageMask = imageMask3(:,:,1);
    imageMask = padarray(imageMask,[1 1],'replicate','both');
%     imageMask = imresize(imageMask, [size(imageTarget, 1) size(imageTarget, 2)]);
% imageMask(460:(end-1),180:570) = 255;
end
% imageTarget = imresize(imageTarget, 0.1);
% imageSource = imresize(imageSource, [size(imageTarget, 1) size(imageTarget, 2)]);
% imageMask = imresize(imageMask, [size(imageTarget, 1) size(imageTarget, 2)]);
end
