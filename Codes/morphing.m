% morphing
clear all;
close all
% clear all
I1 = imread('yogish_specs.jpg');
I2 = imread('indrani.jpg');

I1=imresize(I1,[size(im,1) size(im,2)]);
I2=imresize(I2,[size(im,1) size(im,2)]);
imshow(I2);

if 0
    [pts_I1,pts_I2] = cpselect(I1, I2, 'Wait', true);
    save pts_yd2green;
    save pts_ibgreen;
else
    load pts_gsgreen;
    load pts_ydgreen;
end
I1_copy = I1;
I2_copy = I2;
pts_I1_copy = pts_I1;
pts_I2_copy = pts_I2;

% I1 = imresize(I1, 0.1);
% I2 = imresize(I2, 0.1);
% pts_I1 = pts_I1/10;
% pts_I2 = pts_I2/10;
numPoints = size(pts_I1, 1);

%%
tic
[X,Y] = meshgrid(1:size(I1,2),1:size(I1,1));
X=round(X(:));
Y=round(Y(:));

numMorphImages = 10;
morphID = 0:(1/numMorphImages):1;
morphSet = zeros(size(I1,1),size(I1,2),size(I1,3),numMorphImages);
morphSet(:,:,:,1) = I1;
morphSet(:,:,:,numMorphImages) = I2;
for iter = 1:length(morphID)
    fprintf('##### Running iter: %d #####\n',iter);
    pts_I1_intermediate = (1-morphID(iter)).*pts_I1 + morphID(iter).*pts_I2;
    pts_I2_intermediate = pts_I1_intermediate;
    fprintf('\tCalculating Forward Transform');
    %Forward transform
    b = [pts_I1_intermediate; zeros(3, 2)];
    A = zeros(numPoints + 3, numPoints + 3);
    phi = zeros(numPoints, numPoints);
    for i = 1:numPoints
        for j = (i+1):numPoints
            phi(i,j) = norm(pts_I1(j,:) - pts_I1(i,:));
            phi(i,j) = (phi(i,j)^2)*log(phi(i,j));
            phi(j,i) = phi(i,j);
        end
    end
    A(1:numPoints,1:numPoints) = phi;
    A(numPoints+1:numPoints+3, 1:numPoints) = [pts_I1'; ones(1, numPoints)];
    A(1:numPoints, numPoints+1:numPoints+3) = [pts_I1 ones(numPoints, 1)];
    x_fwd = A\b;
    fprintf('\tEnd\n');
    fprintf('\tCalculating Backward Transform');
    %Backward transform
    b = [pts_I2_intermediate; zeros(3, 2)];
    A = zeros(numPoints + 3, numPoints + 3);
    phi = zeros(numPoints, numPoints);
    
    for i = 1:numPoints
        for j = (i+1):numPoints
            phi(i,j) = norm(pts_I2(j,:) - pts_I2(i,:));
            phi(i,j) = (phi(i,j)^2)*log(phi(i,j));
            phi(j,i) = phi(i,j);
        end
    end
    A(1:numPoints,1:numPoints) = phi;
    A(numPoints+1:numPoints+3, 1:numPoints) = [pts_I2'; ones(1, numPoints)];
    A(1:numPoints, numPoints+1:numPoints+3) = [pts_I2 ones(numPoints, 1)];
    x_bwd = A\b;
    fprintf('\tEnd\n');
    fprintf('\tGenerating warped image locations');
    %Generate warped image
    newPoints_I1 = [X Y];
    newPoints_I2 = [X Y];
    numPointsNew = size(newPoints_I1,1);
    phiNew_I1 = zeros(numPointsNew, numPoints);
    phiNew_I2 = zeros(numPointsNew, numPoints);
    for j = 1:numPoints
        mat1_I1 = repmat(pts_I1(j,:),numPointsNew,1);
        mat2_I1 = newPoints_I1;
        mat3_I1 = mat1_I1 - mat2_I1;
        mat3_I1 = sqrt(sum(mat3_I1.^2,2));
        findNonZero_I1 = find(mat3_I1 ~= 0);
        phiNew_I1(findNonZero_I1,j) = (mat3_I1(findNonZero_I1).^2).*log(mat3_I1(findNonZero_I1));
        
        mat1_I2 = repmat(pts_I2(j,:),numPointsNew,1);
        mat2_I2 = newPoints_I2;
        mat3_I2 = mat1_I2 - mat2_I2;
        mat3_I2 = sqrt(sum(mat3_I2.^2,2));
        findNonZero_I2 = find(mat3_I2 ~= 0);
        phiNew_I2(findNonZero_I2,j) = (mat3_I2(findNonZero_I2).^2).*log(mat3_I2(findNonZero_I2));
    end
    
    ANew = zeros(numPointsNew + 3, numPoints + 3);
    ANew(1:numPointsNew,1:numPoints) = phiNew_I1;
    ANew(numPointsNew+1:numPointsNew+3,1:numPoints) = [pts_I1'; ones(1, numPoints)];
    ANew(1:numPointsNew,numPoints+1:numPoints+3) = [newPoints_I1 ones(numPointsNew, 1)];
    bNew_fwd = ANew*x_fwd;
    bNew_fwd = round(bNew_fwd(1:numPointsNew,:));
    
    ANew = zeros(numPointsNew + 3, numPoints + 3);
    ANew(1:numPointsNew,1:numPoints) = phiNew_I2;
    ANew(numPointsNew+1:numPointsNew+3,1:numPoints) = [pts_I2'; ones(1, numPoints)];
    ANew(1:numPointsNew,numPoints+1:numPoints+3) = [newPoints_I2 ones(numPointsNew, 1)];
    bNew_bwd = ANew*x_bwd;
    bNew_bwd = round(bNew_bwd(1:numPointsNew,:));
    
    xWarp_fwd = bNew_fwd(:,1);
    yWarp_fwd = bNew_fwd(:,2);
    xWarp_fwd = max(min(xWarp_fwd,size(I1,2)),1);
    yWarp_fwd = max(min(yWarp_fwd,size(I1,1)),1);
    
    xWarp_bwd = bNew_bwd(:,1);
    yWarp_bwd = bNew_bwd(:,2);
    xWarp_bwd = max(min(xWarp_bwd,size(I1,2)),1);
    yWarp_bwd = max(min(yWarp_bwd,size(I1,1)),1);
    
    holeImage_I1 = ones(size(I1,1),size(I1,2));
    holeImage_I1(sub2ind(size(I1(:,:,1)),yWarp_fwd,xWarp_fwd)) = 0;
    holeImage_I2 = ones(size(I2,1),size(I2,2));
    holeImage_I2(sub2ind(size(I2(:,:,1)),yWarp_bwd,xWarp_bwd)) = 0;
    warp_I1 = zeros(size(I1,1),size(I1,2),3);
    warp_I2 = zeros(size(I2,1),size(I2,2),3);
    fprintf('\tEnd\n');
    fprintf('\tGenerating warped images with holes');
    for i=1:3 %numChannels
        image_I1 = I1(:,:,i);
        image_I2 = I2(:,:,i);
        warpImage_I1 = zeros(size(I1,1),size(I1,2),1);
        warpImage_I2 = zeros(size(I2,1),size(I2,2),1);
        
        warpImage_I1(uint32(sub2ind(size(warpImage_I1),yWarp_fwd,xWarp_fwd))) = image_I1(uint32(sub2ind(size(image_I1),newPoints_I1(:,2),newPoints_I1(:,1))));
        warpImage_I2(uint32(sub2ind(size(warpImage_I2),yWarp_bwd,xWarp_bwd))) = image_I2(uint32(sub2ind(size(image_I2),newPoints_I2(:,2),newPoints_I2(:,1))));
        
        warp_I1(:,:,i) = warpImage_I1;
        warp_I2(:,:,i) = warpImage_I2;
    end
    warp_I1 = uint8(warp_I1);
    warp_I2 = uint8(warp_I2);
    fprintf('\tEnd\n');
    fprintf('\tHole-filling of warped images');
    %Hole filling
    for i=1:3
        warpImage_I1 = warp_I1(:,:,i);
        warpImage_I2 = warp_I2(:,:,i);
        value_I1 = warpImage_I1(uint32(sub2ind(size(warpImage_I1),yWarp_fwd,xWarp_fwd)));
        value_I2 = warpImage_I2(uint32(sub2ind(size(warpImage_I2),yWarp_bwd,xWarp_bwd)));
        F_I1 = scatteredInterpolant(xWarp_fwd,yWarp_fwd,double(value_I1));
        F_I2 = scatteredInterpolant(xWarp_bwd,yWarp_bwd,double(value_I2));
        
        [queryY,queryX] = ind2sub(size(warpImage_I1),find(holeImage_I1));
        warpImage_I1(find(holeImage_I1)) = F_I1(queryX,queryY);
        warp_I1(:,:,i) = warpImage_I1;
        [queryY,queryX] = ind2sub(size(warpImage_I2),find(holeImage_I2));
        warpImage_I2(find(holeImage_I2)) = F_I2(queryX,queryY);
        warp_I2(:,:,i) = warpImage_I2;
    end
    warp_I1 = uint8(warp_I1);
    warp_I2 = uint8(warp_I2);
    fprintf('\tEnd\n');
    fprintf('\tCross-dissolving of warped images');
    %Cross-dissolve
    morphSet(:,:,:,iter) = (1-morphID(iter))*(double(warp_I1)) + morphID(iter)*(double(warp_I2));
    fprintf('\tEnd\n');
    figure(1),
    subplot(2,3,1),imshow(I1);
    subplot(2,3,4),imshow(I2);
    subplot(2,3,2),imshow(warp_I1);
    subplot(2,3,5),imshow(warp_I2);
    subplot(2,3,3),imshow(uint8(morphSet(:,:,:,iter)));
    %     keyboard;
    %%
    %Adding name to images
    name1 = 'Yogish';
    name2 = 'Indrani';
    image1 = morphSet(:,:,:,iter);
    position = [50,50];
    image1_name1 = insertText(image1,position,name1,'FontSize',18,'BoxOpacity',0);
    image1_name2 = insertText(image1,position,name2,'FontSize',18,'BoxOpacity',0);

    h = fspecial('gaussian',[2*length(morphID)-1,1],1.5);
    h = (1/h(length(morphID)))*h;%normalise middle value to 1;other pixels are accordingly scaled
    
    h_name1 = h(length(morphID):length(h));
    h_name2 = h(1:length(morphID));
    %     weight = h_name1(iter) + h_name2(iter);
    if(iter < (length(morphID)/2))
        newImage = h_name1(iter)*image1_name1 + (1-h_name1(iter))*image1;
    else
        newImage = h_name2(iter)*image1_name2 + (1-h_name2(iter))*image1;
    end
    %     newImage = (h_name1(iter)*image1_name1 + h_name2(iter)*image1_name2 + (1-h_name1(iter))*image1 + (1-h_name2(iter))*image1)/2;%/weight;
    fileName = ['D:\Dropbox\CVFX\Project\Code\Prob3\FaceMorph\yog_ind_green\Name_' num2str(iter) '.jpg'];
    imwrite(uint8(newImage),fileName);
    figure(2),
    imshow(uint8(newImage));
end
toc
