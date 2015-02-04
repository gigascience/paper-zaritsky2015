function [] = speedKymograph(dirname,distanceFromEdge)
% speedKymograph creates a kymograph visualizing spatiotemporal dynamics of
% cells speed.

% speedKymograph(dirname) - traverses over the .tif/.lsm files (each
% holding raw data of one time-lapse experiment) in the *main* directory
% "dirname" and generates a spatiotemporal kymograph for each corresponding
% time-lapse experiment. The kymograph is calculated for "phase 1" in the
% healing process, until first touch of cells from opposing borders of the
% wound. The kymograph is saved in both the corresponding experiment's
% directory and in a directory named speedKymographs at the main directory.
% Both the kymograph data (.mat) and a visualization (.bmp) files are
% saved.

% Inputs: 
%   dirname - where the data files are. It is expected to have a .tif/.lsm
%               file for each time-lapse experiment + the corresponding directories the
%               hold the actual processed data needed to produce the kymograph.
%   distanceFromEdge - what is the maximal spacial distance from the wound
%                       to consdier for the kymograph (default 300 um)
%
%   Assaf Zaritsky, November 2014
%   Email: assafzar@gmail.com

if nargin < 3
    distanceFromEdge = 300; % um
end

allKymographsDirname = [dirname filesep 'speedKymographs'];
if ~exist(allKymographsDirname,'dir')
    mkdir(allKymographsDirname);
end

fnames = dir(dirname);
nfnames = length(fnames);

for i = 3 : nfnames
    fname = fnames(i).name;
    [pathstr, name, ext] = fileparts(fname);    
    
    if ~strcmp(ext, '.tif') && ~strcmp(ext, '.lsm')
        continue;
    end
    
    if exist([allKymographsDirname filesep name '_speedKymograph.mat'],'file')
        continue;
    end
    
    curDirname = [dirname filesep name filesep];
    
    load([curDirname 'experimentParams.mat']);% pixelSize, timePerFrame, maxTime, timePhase1
    
    imgsDirname = [curDirname 'images' filesep];
    mfDirname = [curDirname 'MF' filesep];
    roiDirname = [curDirname 'ROI' filesep];
    gcDirname = [curDirname 'GC' filesep];
    
    patchSize = ceil(15.0/pixelSize); % 15 um in pixels
    strips =  patchSize : patchSize : distanceFromEdge; % um
    nstrips = length(strips);

    nTime = timePhase1;   
    
    speedKymograph = nan(nstrips,nTime);

    for t = 1 : nTime
        roiFname = [roiDirname sprintf('%03d',t) '_roi.mat']; % ROI is a mask that labels each pixel as cellular (1) or background (0)
        mfFname = [mfDirname sprintf('%03d',t) '_mf.mat']; % dxs, dys define together the velocity at each pixel. The units are in um/hour. 
        
        load(roiFname);
        sizeTH = 3600/(pixelSize.^2); % in pixels based on 60um x 60um
        ROI = cleanROI(ROI,sizeTH); % Discrad small excesses, both in cellular regions and in background
        load(mfFname);
        
        speed = sqrt(dxs.^2 + dys.^2); % um / hour
        DIST = bwdist(~ROI);
    
        for d = 1 : nstrips
            inDist = ...
                (DIST > strips(d)-patchSize) & ...
                (DIST < strips(d)) & ...
                ~isnan(speed);
            speedInStrip = speed(inDist);
            speedKymograph(d,t) = mean(speedInStrip);            
        end
    end

    save([curDirname 'speedKymograph.mat'],'speedKymograph');
    save([allKymographsDirname filesep name '_speedKymograph.mat'],'speedKymograph');
    
    caxisVals = [0 30];
    outFname = [curDirname 'speedKymograph.bmp'];
    
    plotKymograph(speedKymograph,nTime,timePerFrame,patchSize,distanceFromEdge,caxisVals,outFname);
    plotKymograph(speedKymograph,nTime,timePerFrame,patchSize,distanceFromEdge,caxisVals,[allKymographsDirname filesep name '_speedKymograph.bmp']);

    close all;
end
end


function [] = plotKymograph(kymograph,nTime,timePerFrame,patchSize,distanceFromEdge,caxisVals,outFname)

maxTime = nTime * timePerFrame;
maxTimeMod = (maxTime - mod(maxTime,100));
maxDistMod = (distanceFromEdge - mod(distanceFromEdge,100));

xTick = 1:(200/timePerFrame):((maxTimeMod/timePerFrame)+1);
xTickLabel = 0:200:maxTimeMod;
yTick = 1:(100/patchSize):((distanceFromEdge/patchSize)+1);
yTickLabel = 0:100:maxDistMod;

h = figure;
imagescnan(kymograph);
hold on;
caxis(caxisVals); colorbar;
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[1,maxTime/timePerFrame]);
set(haxes,'XTick',xTick);
set(haxes,'XTickLabel',xTickLabel);
set(haxes,'YTick',yTick);
set(haxes,'YTickLabel',yTickLabel);
set(haxes,'FontSize',32);
xlabel('Time (minutes)','FontSize',32); ylabel('Distance from edge (\mum)','FontSize',32);
hold off;
eval(sprintf('print -dbmp16m  %s', outFname));
end

function [ROI_OUT] = cleanROI(ROI,sizeTH)
tmp1 = excludeExcesses(ROI,sizeTH);
tmp2 = excludeExcesses(~ROI,sizeTH);
ROI_OUT = tmp1 & ~tmp2; 
end

function [ROI_OUT] = excludeExcesses(ROI,sizeTH)
[L, num] = bwlabel(ROI,4);
ROI_OUT = false(size(ROI));
for i = 1 : num
    cur = L == i;
    if sum(cur(:)) > sizeTH
        ROI_OUT(cur) = true;
    end
end
end

% function [ROI_OUT] = fillHoles(ROI)
% tmp1 = imfill(ROI,'holes');
% tmp2 = imfill(~ROI,'holes');
% ROI_OUT = tmp1 & ~tmp2; 
% end