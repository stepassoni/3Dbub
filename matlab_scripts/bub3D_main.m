%% IMPORT AND CLEANUP
%tic
clear all
clc
close all

pattern='coarse';   % Flow pattern analyzed, used only for defining the working directory
smoothing=1;    % 1: use of smoothing function; 0: no smoothing

% Import the original frames and the binary mask, set directories
% according to your environment
inputMask = dir(['../../visualization_data/',pattern,'/segmentation/*.tif']);
inputOrig = dir(['../../visualization_data/',pattern,'/original_cropped/*.tif']);


fileNames = {inputMask.name};
testCases=regexprep(fileNames,'seg_','');
testCases=regexprep(testCases,'_multi.tif','')';
load('px_per_mm.mat');  % px-to-mm conversion factor

for file=1:length(testCases)
    inputMask_path=[inputMask(file).folder,'/',inputMask(file).name];
    inputOrig_path=[inputOrig(file).folder,'/',inputOrig(file).name];

    info = imfinfo(inputMask_path);
    numberOfPages = length(info);

    for k=1:numberOfPages
        BW(:,:,k)=imbinarize(imread(inputMask_path,k));
        filled(:,:,k) = imfill(BW(:,:,k), 'holes');
        BW(:,:,k) = bwareaopen(BW(:,:,k), 4); %remove noise
        holes(:,:,k) = filled(:,:,k) & ~BW(:,:,k);
        bigholes(:,:,k) = bwareaopen(holes(:,:,k), 200);
        smallholes(:,:,k) = holes(:,:,k) & ~bigholes(:,:,k);
        BW(:,:,k) = BW(:,:,k) | smallholes(:,:,k);
        I0(:,:,k)=imread(inputOrig_path,k);
    end


    %% DEPTH MAP CALCULATION

    [depth_pos depth_neg depth_tot]=depthmap(BW(:,:,1),2.5,8.9,-63,px_per_mm,4,0); % Calculate channel depth map

    clearvars info

    %% BLOB ANALYSIS
    clc
    close all

    clear V diam
    target=BW;
    numberOfPages=size(target,3);
    frameLimit=numberOfPages;

    if smoothing
        for i=1:frameLimit
            [V{i}, alpha(i), alphaArea(i), isSmall{i}, diam{i}, LFR(i)]=bubAnalysis(target(:,:,i),I0(:,:,i),px_per_mm,depth_tot,depth_pos, depth_neg,0,0,1,3,0);
        end
    else
        for i=1:frameLimit
            [V{i}, alpha(i), alphaArea(i), isSmall{i}, diam{i}, LFR(i)]=bubAnalysis(target(:,:,i),I0(:,:,i),px_per_mm,depth_tot,depth_pos, depth_neg,0,0,0,3,0);
        end

    end

    alphaAVG=mean(alpha);
    alphaAreaAVG=mean(alphaArea);
    LFR_AVG=mean(LFR);

    % VOLUME STATISTICS
    clear Vol
    close all
    timeAvgVol=zeros(size(target,3),1);
    if size(target,3)>1
        for i=1:size(target,3)
            if i==1
                Vol=V{i};
                allDiam=diam{i}';
                avgVol(i)=mean(V{i});
                timeAvgVol(i)=mean(V{i});
            else
                Vol=[Vol;V{i}];
                allDiam=[allDiam;diam{i}'];
                avgVol(i)=mean(V{i});
                timeAvgVol(i)=mean(V{i});
                timeAvgVol(i)=mean([timeAvgVol(i) timeAvgVol(i-1)]);
            end
        end
    else
        %         minVol=min(cell2mat(V));
        %         maxVol=max(cell2mat(V));
        %         avgVol=mean(cell2mat(V));
        %         maxDev=((maxVol-avgVol)/avgVol);
        %         minDev=(abs(minVol-avgVol)/avgVol);
        %         Vol=cell2mat(V);
    end
    pd=fitdist(allDiam,'lognormal');
    Dmean=exp(pd.mu+(pd.sigma^2)/2);
    Dvariance=exp(2*pd.mu+pd.sigma^2)*(exp(pd.sigma^2)-1);
    % Data tabel
    Results(file).testCase=testCases(file);
    Results(file).alpha=alpha;
    Results(file).alphaAVG=alphaAVG;
    Results(file).alphaArea=alphaArea;
    Results(file).alphaAreaAVG=alphaAreaAVG;
    Results(file).allDiam=allDiam;
    Results(file).D_mean=Dmean;
    Results(file).D_variance=Dvariance;
    Results(file).avgVol=avgVol;
    Results(file).LFR=LFR_AVG;


end

if smoothing
    save(['results_smoothing_',pattern,'.mat'],'Results');
else
    save(['results_',pattern,'.mat'],'Results');
end
%toc
