% correct turbulence using region-based fusion (pixel-based is also provided) and roi
% Created by Nantheera Anantrasirichai,  The University of Bristol
% 20/04/2012 n.anantrasirichai@bristol.ac.uk
% Please cite
% Anantrasirichai, N.; Achim, A.; Kingsbury, N.G.; Bull, D.R., "Atmospheric
% Turbulence Mitigation Using Complex Wavelet-Based Fusion," Image Processing,
% IEEE Transactions on , vol.22, no.6, pp.2398-2408, June 2013

clear all

addpath('codes');

% -------------------------------------------------------------------------
% input distorted sequence
% -------------------------------------------------------------------------
dirnameroot = 'images\';


dirname = dirnameroot;% [dirnameroot,'distortion',num2str(dist),'\'];
extfile = 'tif';%'png';
startFrame = 1;
totalFrame = 50; %length(dir([dirnameroot,'*.',extfile]));

% -------------------------------------------------------------------------
% control parameters
% -------------------------------------------------------------------------
resizeRatio = 0.5;
fusionMethod = 'pixel';     % 'region', 'pixel'
                            % 'pixel' is faster, 'region' is better for noisy
refFrameType = 'average';   % 'maxGradient', 'maxHP'
doFrameSelection = false;   % true if i) speed up, ii) need stabilisation
selectROI = false;          % true if ROI is significantly shifted between succesive frames
sharpMethod = 'gradient';   %'gradient' or 'dtcwt'
doPostprocess = false;      % true - sharpening and denoising;
doEnhance = true;           % false;

mseGain = 1;
gradGain1 = 1;
gradGain2 = 1;
infoGain2 = 20;

levels = 4;

numFrameRegis = 50;
maxFrameused = 25;
clipLimit = 0.002;
useBigAreaInfo = false;     % true;


% -------------------------------------------------------------------------
% load input sequence
% -------------------------------------------------------------------------
[input, inputU, inputV] = loadInput(dirname, extfile,startFrame, totalFrame, resizeRatio);
% for k = 1:length(mov)
%     img = mov.cdata;
%     if size(img,3)==3
%         img = double(rgb2ycbcr(img));
%         inputU(:,:,k) = img(:,:,2);
%         inputV(:,:,k) = img(:,:,3);
%     end
%     input(:,:,k) = double(img(:,:,1));
% end
disp('Done read inputs')
% -------------------------------------------------------------------------
% find reference frame
% -------------------------------------------------------------------------
avgFrame = findRefFrame(input, refFrameType);

% -------------------------------------------------------------------------
% select good frames
% -------------------------------------------------------------------------
if doFrameSelection
    [height, width, totalFrame] = size(input);
    rangei = 1:height;
    rangej = 1:width;
    rangek = 1:totalFrame;
    [valcostSelect, indcostSelect] = findGoodFrame(avgFrame, input, rangei, rangej, rangek, totalFrame, sharpMethod, 0, ...
    mseGain, gradGain1, 0);
else
    indcostSelect = 1:totalFrame;
end

% -------------------------------------------------------------------------
% select good frames from ROI
% -------------------------------------------------------------------------
if selectROI
    [x_pos, y_pos] = selectPoint(avgFrame,2,'click roughly ROI');
    rangei = round(y_pos(1):y_pos(2));
    rangej = round(x_pos(1):x_pos(2));
    [valcostSelect, indcostSelectupdate] = findGoodFrame(input(:,:,indcostSelect(1)),input,...
        rangei, rangej, indcostSelect, numFrameRegis, sharpMethod, 0, ...
        mseGain, gradGain2, infoGain2);
else
    indcostSelectupdate = indcostSelect;
end
if doFrameSelection || selectROI
    halfval = (mseGain+gradGain2)/2; % half of max costSelect
    [~,numFrametoUse] = min(abs(valcostSelect-halfval));
    numFrametoUse = min(numFrametoUse,maxFrameused);

    % reduce input to only good selected frames
    input = input(:,:,indcostSelectupdate(1:numFrametoUse));
    if ~isempty(inputU)
        inputU = inputU(:,:,indcostSelectupdate(1:numFrametoUse));
        inputV = inputV(:,:,indcostSelectupdate(1:numFrametoUse));
    end
end
disp('Done frame selection')
% -------------------------------------------------------------------------
% registration
% -------------------------------------------------------------------------
if isempty(inputU)
    [xrest] = Nick_regis(input, input, input, levels);
else
    [xrest, xrestU, xrestV] = Nick_regis(input, inputU, inputV, levels);
end
disp('Done registration')
% -------------------------------------------------------------------------
% fusion
% -------------------------------------------------------------------------
if strcmp(fusionMethod, 'pixel')
    [zrest, zrestsmooth] = Nick_pixel_fuse(xrest, levels);
else
    [zrest, zrestsmooth, zrestReg] = fuseRegionROI(xrest, levels);
end
disp('Done fusion')
% -------------------------------------------------------------------------
% contrast enhancement
% -------------------------------------------------------------------------
if doPostprocess
    zrest = postDenoiseSharpen(zrest,levels);
end

% -------------------------------------------------------------------------
% contrast enhancement
% -------------------------------------------------------------------------
if doEnhance
    zrest = adapthisteq(uint8(zrest), 'clipLimit',clipLimit);
end
% -------------------------------------------------------------------------
% final result
% -------------------------------------------------------------------------
figure; imshow([input(:,:,1) zrest]); title('Distorted image and corrected image')
if ~isempty(inputU)
    zrest(:,:,2) = mean(inputU,3);
    zrest(:,:,3) = mean(inputV,3);
    zrest = ycbcr2rgb(uint8(zrest));
end
if strcmp(fusionMethod, 'pixel')
    imwrite(uint8(zrest), 'results\fusionResultPixel.png','png');
else
    imwrite(uint8(zrest), 'results\fusionResultRegion.png','png');
end
