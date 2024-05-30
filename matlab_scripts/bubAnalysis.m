function [V, alpha, alphaArea, isSmall, allBlobDiam,LFR]=bubAnalysis(BW,I0, px_per_mm,depth_tot,depth_pos,depth_neg,plotNum, plotDim, smoothing, sigma, plot3D)

% INPUTS:
%
% BW = segmentation mask (1 frame)
%
% I0 = original frame (1 frame)
%
% px_per_mm = conversion factor px/mm
%
% depth_tot = depth map of the channel (sum of positive and negative
% displacement in absolute value)
%
% plotNum = 1 if you want to plot bubble number on top of segmentation mask
%
% plotDim = 1 if you want to plot bubble bubble centroids with different
% symbols according to volume reconstruction method
%
% plot3D = 1 if you want to plot the 3D reconstruction
textFontSize=12;
[labeledImage, numberOfBlobs] = bwlabel(BW, 4);
%coloredLabels = label2rgb (labeledImage, 'hsv', 'k', 'shuffle');
%imshow(coloredLabels);


%Calculate Blobs geometrical properties

props = regionprops(labeledImage, BW, 'EquivDiameter','Centroid','MajorAxisLength','MinorAxisLength','Orientation','Area','PixelIdxList');
numberOfBlobs = numel(props);
allBlobDiam = [props.EquivDiameter];
allBlobCentroids = vertcat(props.Centroid);
allBlobMajAxis = [props.MajorAxisLength];
allBlobMinAxis = [props.MinorAxisLength];
allBlobOrientation = [props.Orientation];
centroidsX = allBlobCentroids(:, 1);			% Extract out the centroid x values into their own vector.
centroidsY = allBlobCentroids(:, 2);			% Extract out the centroid y values into their own vector.


if plotNum==1
    figure(1)
    imshow(BW)
    for k = 1 : numberOfBlobs           % Loop through all blobs.
        % Place the blob label number at the centroid of the blob.
        text(centroidsX(k), centroidsY(k), num2str(k), 'FontSize', textFontSize, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','color','r');
    end
end

if plotDim==1
    figure(2)
    imshow(BW);
    isSmall=0;
    hold on;
    for k = 1 : numberOfBlobs           % Loop through all  blobs.
        % Identify if blob #k is small
        smallBlob = allBlobDiam(k)/px_per_mm < mean(depth_tot(labeledImage==k));%|| allBlobDiam(k)/px_per_mm<2.5; % Criterion for small blob classification

        if smallBlob
            plot(centroidsX(k), centroidsY(k), 'r+', 'MarkerSize', 15, 'LineWidth', 2);
            isSmall(k)=1;
        else
            plot(centroidsX(k), centroidsY(k), 'bx', 'MarkerSize', 15, 'LineWidth', 2);
            isSmall(k)=0;
        end
    end
else
    isSmall=0;
    for k = 1 : numberOfBlobs           % Loop through all  blobs.
        % Identify if blob #k is small
        smallBlob = allBlobDiam(k)/px_per_mm < mean(depth_tot(labeledImage==k)); %|| allBlobDiam(k)/px_per_mm<2.5; % Criterion for small blob classification

        if smallBlob
            isSmall(k)=1;
        else
            isSmall(k)=0;
        end
    end
end

% Volume calculation and 3D plot
cont=0;

if plot3D==1
    figure(1)
    hold on
    V=zeros(numberOfBlobs,1);
    [x, y, z] = size(BW);
    extruded_pos = zeros(x, y, z);
    extruded_neg = zeros(x, y, z);
    for k = 1 : numberOfBlobs
        if isSmall(k)
            cont=cont+1;
            [X Y Z]=ellipsoid(centroidsX(k),centroidsY(k),0,allBlobMajAxis(k)/2,allBlobMinAxis(k)/2,allBlobMinAxis(k)/2);
            V(k)=(4/3*pi*allBlobMajAxis(k)/2*allBlobMinAxis(k)/2*allBlobMinAxis(k)/2)/(px_per_mm)^3;
            alpha=allBlobOrientation(k);
            rotz=[cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
            % Create transformation matrix
            R = [cosd(alpha), -sind(alpha), 0;
                sind(alpha), cosd(alpha), 0;
                0, 0, 1];
            X = X - centroidsX(k);
            Y = Y - centroidsY(k);
            Z = Z;
            % Apply rotation to ellipsoid coordinates
            XYZ = [X(:), Y(:), Z(:)];
            XYZ_rot = (XYZ*R);

            % Reshape rotated coordinates back into ellipsoid shape
            x_rot = reshape(XYZ_rot(:,1), size(X))+centroidsX(k);
            y_rot = reshape(XYZ_rot(:,2), size(Y))+centroidsY(k);
            z_rot = reshape(XYZ_rot(:,3), size(Z));

            s=surf(x_rot, y_rot, z_rot);

            keeperIndexes = k;
            keeperBlobsImage = ismember(labeledImage, keeperIndexes);
            % extrude each white pixel in the mask
            % obtain indices of white pixels
        else
            keeperIndexes = k;
            keeperBlobsImage = ismember(labeledImage, keeperIndexes);
            % extrude each white pixel in the mask
            % obtain indices of white pixels
            [rows, cols] = find(keeperBlobsImage == 1);

            if smoothing == 1
                %SMOOTHING
                % calculate extruded depth
                distMap = bwdist(~keeperBlobsImage);

                delta=mean(depth_tot(distMap==1));
                lc=2.7; %capillary length

                if delta>lc
                    weightFunc=1-exp(-lc.*((distMap-1)./px_per_mm));
                else
                    weightFunc=1-exp(-delta.*((distMap-1)./px_per_mm));
                end

                % Modify depth map by multiplying with weighting function
                extruded_pos(keeperBlobsImage) = depth_pos(keeperBlobsImage)-(abs(depth_pos(keeperBlobsImage).* (1-weightFunc(keeperBlobsImage))));
                extruded_neg(keeperBlobsImage) = depth_neg(keeperBlobsImage)+ (abs(depth_neg(keeperBlobsImage).* (1-weightFunc(keeperBlobsImage))));
            else
                %ORIGINAL
                extruded_pos(keeperBlobsImage) = depth_pos(keeperBlobsImage);
                extruded_neg(keeperBlobsImage) = depth_neg(keeperBlobsImage);
            end
            extruded_depth=abs(extruded_pos-extruded_neg);
            V(k)=(abs(sum(extruded_pos(find(labeledImage==k))) - sum(extruded_neg(find(labeledImage==k)))))/px_per_mm^2;

        end
    end

    % Plot big blobs
    [m, n] = size(BW);
    [X, Y] = meshgrid(1:n, 1:m);
    %extruded_pos(extruded_pos==0)=NaN;
    %extruded_neg(extruded_neg==0)=NaN;

    N = 1;    % number of pixels to grow the borders
    se = ones(2*N + 1, 2*N + 1);
    if size(find(isSmall==0))>0
        alphamask = imdilate(extruded_depth>0,se);
%hold on
        h=surf(X,Y,double(extruded_pos)*px_per_mm);
        h1=surf(X,Y,double(extruded_neg)*px_per_mm);
        h.AlphaData = alphamask*255;
        h.FaceAlpha = 'interp';
        h1.AlphaData = alphamask*255;
        h1.FaceAlpha = 'flat';
    end
    axis equal
    xlim([0,1200])
    ylim([0,460])
    shading interp
    %colormap jet
    colorbar
    lightangle(-60,70)
    %clim([2.5 25])
    %overlab grey image
    set(gca, 'YDir','reverse')
    %hold on
    imshow(cat(3, I0,I0,I0))
    %hold off
   

else
    V=zeros(numberOfBlobs,1);
    [x, y, z] = size(BW);
    extruded_pos = zeros(x, y, z);
    extruded_neg = zeros(x, y, z);
    for k = 1 : numberOfBlobs
        if isSmall(k)
            cont=cont+1;
            [X Y Z]=ellipsoid(centroidsX(k),centroidsY(k),0,allBlobMajAxis(k)/2,allBlobMinAxis(k)/2,allBlobMinAxis(k)/2);
            V(k)=(4/3*pi*allBlobMajAxis(k)/2*allBlobMinAxis(k)/2*allBlobMinAxis(k)/2)/(px_per_mm)^3;
            alpha=allBlobOrientation(k);
            rotz=[cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
            % Create transformation matrix
            R = [cosd(alpha), -sind(alpha), 0;
                sind(alpha), cosd(alpha), 0;
                0, 0, 1];
            X = X - centroidsX(k);
            Y = Y - centroidsY(k);
            Z = Z;
            % Apply rotation to ellipsoid coordinates
            XYZ = [X(:), Y(:), Z(:)];
            XYZ_rot = (XYZ*R);

            % Reshape rotated coordinates back into ellipsoid shape
            x_rot = reshape(XYZ_rot(:,1), size(X))+centroidsX(k);
            y_rot = reshape(XYZ_rot(:,2), size(Y))+centroidsY(k);
            z_rot = reshape(XYZ_rot(:,3), size(Z));

        else
            keeperIndexes = k;
            keeperBlobsImage = ismember(labeledImage, keeperIndexes);
            % extrude each white pixel in the mask
            % obtain indices of white pixels
            [rows, cols] = find(keeperBlobsImage == 1);

            if smoothing == 1
                %SMOOTHING
                % calculate extruded depth
                distMap = bwdist(~keeperBlobsImage);
               delta=mean(depth_tot(distMap==1));
                lc=2.7; %capillary length

                if delta>lc
                    weightFunc=1-exp(-lc.*((distMap-1)./px_per_mm));
                else
                    weightFunc=1-exp(-delta.*((distMap-1)./px_per_mm));
                end

                % Modify depth map by multiplying with weighting function
                extruded_pos(keeperBlobsImage) = depth_pos(keeperBlobsImage)-(abs(depth_pos(keeperBlobsImage).* (1-weightFunc(keeperBlobsImage))));
                extruded_neg(keeperBlobsImage) = depth_neg(keeperBlobsImage)+ (abs(depth_neg(keeperBlobsImage).* (1-weightFunc(keeperBlobsImage))));
            else
                %ORIGINAL
                extruded_pos(keeperBlobsImage) = depth_pos(keeperBlobsImage);
                extruded_neg(keeperBlobsImage) = depth_neg(keeperBlobsImage);
            end
            %             end
            extruded_depth=abs(extruded_pos-extruded_neg);
            V(k)=abs(sum(extruded_pos(find(labeledImage==k))) - sum(extruded_neg(find(labeledImage==k))))/px_per_mm^2;

        end
    end
end

[labeledImageL, numberOfBlobsL] = bwlabel(imcomplement(BW), 4);
propsL = regionprops(labeledImageL, imcomplement(BW), 'Area');

%Calculate Alpha
Vtot=sum(sum(depth_tot))/px_per_mm^2;
Vgas=sum(V);
alpha=Vgas/(Vtot);
alphaArea=sum([props.Area])/(sum([props.Area])+sum([propsL.Area]));


%LFR
% Filter blobs based on major axis length
minMajorAxisLength = 8.9*px_per_mm*2.5;
validBlobs = [props.MajorAxisLength] > minMajorAxisLength;

% Get the pixel indices for valid blobs
validPixelIdx = vertcat(props(validBlobs).PixelIdxList);
LFarea=sum(vertcat(props(validBlobs).Area));
LFR=LFarea/(size(BW,1)*size(BW,2));
% Create a binary image for valid blobs
filteredMask = false(size(BW));
filteredMask(validPixelIdx) = true;

% Display the filtered blobs (uncomment if needed)
% figure();
% imshow(filteredMask);
% title('Filtered Blobs');
allBlobDiam=allBlobDiam./px_per_mm;
