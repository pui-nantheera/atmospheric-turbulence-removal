function [zrestReg, zrestsmooth, zrest, regMap] = fuseRegionROI(input, levels, ...
    segmentFusedType,zgain,lowpassfuse, numframeEachReg, activity, weight, alllevel, phaseActivity)

% get dimension
[height,width,totalFrame] = size(input);

% -------------------------------------------------------------------------
% process parameters
% -------------------------------------------------------------------------
if nargin < 2
    levels = 4;
end
if nargin < 3
    segmentFusedType = 'Joint'; %'Unimode' or 'Joint'
end
if nargin < 4
    % Set subband gains with some pre-emphasis for reconstructing output image, zrest.
    zgain = ones(6,1)*[3 1.8 1.3 1];
else
    zgain = ones(6,1)*zgain;
end
if nargin < 5
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Method to fuse low pass
    lowpassfuse = 'average'; % 'average' or 'variance' or 'median'
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Number of frames to fuse for high pass
    numframeEachReg = 1;%3; % if >1, using weighted average when fusing image
    numframeEachReg = min(numframeEachReg,size(input,3));
    % Method to find decision for high pass - which frame for each region is the best to use for fusion
    activity = 'average'; % 'average' or 'variance' or 'entropy'
    % Method to fuse high pass
    weight = 'maximum'; % 'maximum' or 'average' or 'linear'
    % All level (alllevel = 1) or individual (alllevel = 0) or subband (alllevel = -1)
    alllevel = 1;
    % phase activity
    phaseActivity = 'average'; % 'average' or 'median'
end

% -------------------------------------------------------------------------
% Region segmentation
% -------------------------------------------------------------------------
if totalFrame>5
    if strcmp(segmentFusedType,'Joint')
        [overlay,intolay,regMap,intmap,gsurf1] = Segment_Rob_forManyFrames(input(:,:,min(10,size(input,3))), segmentFusedType, levels-1);
    else
        regMap = ones(height,width);
        for f = 1:totalFrame
            [overlay,intolay,regMap1,intmap,gsurf1] = Segment_Rob(input(:,:,f), segmentFusedType, levels-1);
            regMap = regMap.*prod(regMap1,3);
        end
    end
else
    [overlay,intolay,regMap,intmap,gsurf1] = Segment_Rob(input, segmentFusedType, levels-1);
end

tic
% -------------------------------------------------------------------------
% segment for each wavelet level
% -------------------------------------------------------------------------
if (size(regMap,3)>1)
    % combine all maps
    combineMap = prod(regMap,3);
    % relabel
    regMap = bwlabel(combineMap);
end
regMap = reIndexMap(regMap);
% decomposition of segmented maps
regDecompMap = LabDecompProt(regMap, levels, [height width]);
% reorder
if alllevel~=-1
    for k = 1:levels
        regDecompMap{k} = repmat(regDecompMap{k},[1 1 6]);
    end
end
totalReg = max(regMap(:));

% -------------------------------------------------------------------------
% initialising
% -------------------------------------------------------------------------
Zrg = cell(levels,1); % for region
% store absolute values of each region which combine all levels
magReg = zeros(totalFrame,totalReg);
% low pass accumulators
if strcmp(lowpassfuse,'median')
    Zl = zeros(size(regDecompMap{levels-1},1),size(regDecompMap{levels-1},2),totalFrame);
else
    Zl = zeros(size(regDecompMap{levels-1},1),size(regDecompMap{levels-1},2));
end
% phase accumulators
Zh = cell(levels,1);
% -------------------------------------------------------------------------
% Using DT-CWT for choosing lucking region for high pass bands
% -------------------------------------------------------------------------
% loop for each frame to find avarage magnitude of high frequency of each
% region
tempActReg = zeros(levels,totalReg);
if alllevel==1
    tempActReg = zeros(levels,totalReg);
    if numframeEachReg > 1
        tempActReg = zeros(levels,totalReg,numframeEachReg);
        frameActReg = ones(levels,totalReg,numframeEachReg);
    end
elseif alllevel==-1
    tempActReg = zeros(levels,totalReg,6);
    if numframeEachReg > 1
        tempActReg = zeros(levels,totalReg,6,numframeEachReg);
        frameActReg = ones(levels,totalReg,6,numframeEachReg);
    end
end
clear Xl Xh 
for f = 1:totalFrame
    f
    % transform
    [Xl{f},Xh{f}] = dtwavexfm2(double(input(:,:,f)),levels,'near_sym_b','qshift_d');
    % for high pass
    if f==1
        Zrg = Xh{f};
        Zl  = Xl{f};
    else
        if strcmp(lowpassfuse,'median')
            Zl(:,:,f) = Xl{f};
        else
            % sum low pass for averaging
            Zl = Zl + Xl{f};
        end
    end
    
    % each level
    for k = 1:levels
        % Sum the phase vectors.
        if strcmp(phaseActivity,'median')
            Zh{k} = cat(4,Zh{k},Xh{f}{k});
        else
            if f==1
                Zh{k} = Xh{f}{k};
            else
                Zh{k} = Zh{k} + Xh{f}{k};
            end
        end
        % for each region
        for reg = 1:totalReg
            if alllevel == 1
                % current data
                curdata = Xh{f}{k}.*(regDecompMap{k}==reg);
                curdata(regDecompMap{k}~=reg) = [];
                if strcmp(activity,'average')
                    % find average absolute magnitude in the region
                    magReg(f,reg) = magReg(f,reg) + sum(abs(curdata(:)));
                elseif strcmp(activity,'variance')
                    % find varience
                    magReg(f,reg) = magReg(f,reg) + sum(var(curdata(:)));
                elseif strcmp(activity,'entropy')
                    % find Shannon entropy
                    curdata = curdata + (curdata==0).*(10^-10);
                    magReg(f,reg) = magReg(f,reg) + sum((curdata(:).^2).*(log(curdata(:).^2)));
                end
            else
                if alllevel==-1
                    % for each subband
                    for sub = 1:6
                        curdata = Xh{f}{k}(:,:,sub).*(regDecompMap{k}==reg);
                        curdata(regDecompMap{k}~=reg) = [];
                        curSumAbs = sum(abs(curdata(:)));
                        if numframeEachReg > 1
                            [val, ind] = min(tempActReg(k,reg,sub,:));
                            if curSumAbs > val
                                tempActReg(k,reg,sub,ind) = curSumAbs;
                                frameActReg(k,reg,sub,ind) = f;
                            end
                        else
                            if curSumAbs > tempActReg(k,reg,sub)
                                tempActReg(k,reg,sub) = curSumAbs;
                                ind = find(regDecompMap{k}==reg) + (sub-1)*size(regDecompMap{k},1)*size(regDecompMap{k},2);
                                Zrg{k}(ind) = Xh{f}{k}(ind);
                            end
                        end
                    end
                else
                    % current data
                    curdata = Xh{f}{k}.*(regDecompMap{k}==reg);
                    % choose the maximam
                    curdata(regDecompMap{k}~=reg) = [];
                    curSumAbs = sum(abs(curdata(:)));
                    if numframeEachReg > 1
                        [val, ind] = min(tempActReg(k,reg,:));
                        if curSumAbs > val
                            tempActReg(k,reg,ind) = curSumAbs;
                            frameActReg(k,reg,ind) = f;
                        end
                    else
                        if curSumAbs > tempActReg(k,reg)
                            tempActReg(k,reg) = curSumAbs;
                            Zrg{k}(regDecompMap{k}==reg) = Xh{f}{k}(regDecompMap{k}==reg);
                        end
                    end
                end
            end
        end
    end
end

if numframeEachReg>1
    % find weight coeff
    if strcmp(weight,'average')
        wcoef = 1/numframeEachReg;
        wcoef = repmat(wcoef,[1 numframeEachReg]);
    elseif strcmp(weight,'linear')
        wcoef = 1:(-1/numframeEachReg):0;
        wcoef = wcoef/sum(wcoef);
    end
else
    % choose the region according to magReg
    [~,indseg] = max(magReg,[],1);
end

%%
% fuse method
for k = 1:levels
    if (alllevel==1)||(numframeEachReg > 1)
        Zrg{k} = Zrg{k}*0;
        % each region
        for reg = 1:totalReg
            if numframeEachReg>1
                if (alllevel==-1)
                    for sub = 1:6
                        datacur = tempActReg(k,reg,sub,:);
                        [~,indN] = sort(datacur,'descend');
                        ind = find(regDecompMap{k}==reg) + (sub-1)*size(regDecompMap{k},1)*size(regDecompMap{k},2);
                        for ord = 1:numframeEachReg
                            if strcmp(weight,'average') || strcmp(weight,'linear')
                                Zrg{k}(ind) = Zrg{k}(ind) + wcoef(ord)*Xh{indN(ord)}{k}(ind);
                            else
                                Zrg{k}(ind) = max(Zrg{k}(ind),Xh{indN(ord)}{k}(ind));
                            end
                        end
                    end
                else
                    % find first N max frames
                    if (alllevel==1)
                        datacur = magReg(:,reg);
                    else
                        datacur = tempActReg(k,reg,:);
                    end
                    [~,indN] = sort(datacur,'descend');
                    for ord = 1:numframeEachReg
                        if strcmp(weight,'average') || strcmp(weight,'linear')
                            Zrg{k}(regDecompMap{k}==reg) = Zrg{k}(regDecompMap{k}==reg) + wcoef(ord)*Xh{indN(ord)}{k}(regDecompMap{k}==reg);
                        else
                            Zrg{k}(regDecompMap{k}==reg) = max(Zrg{k}(regDecompMap{k}==reg),Xh{indN(ord)}{k}(regDecompMap{k}==reg));
                        end
                    end
                end
            else
                % use only the one max frame
                Zrg{k}(regDecompMap{k}==reg) = Xh{indseg(reg)}{k}(regDecompMap{k}==reg);
            end
        end
    end
    % Generate unit vectors with correct phases.
    if strcmp(phaseActivity,'median')
        Zh{k} = median(Zh{k},4);
    end
    Zh{k} = Zh{k} ./ abs(Zh{k});
    Zrg1{k} = Zh{k}.*abs(Zrg{k}); % Scale by the max magnitudes.
end

% -------------------------------------------------------------------------
% Average low pass
% -------------------------------------------------------------------------
if strcmp(lowpassfuse,'variance')
    % fusing low pass using variance
    regDecompMapLL = regDecompMap{levels-1}(:,:,1);
    SegEle = unique(regDecompMapLL);
    SegLen = length(SegEle);
    varianceL = zeros(totalFrame,SegLen);
    % find varience
    for i=1:SegLen
        seg = find(regDecompMapLL == SegEle(i));
        for f = 1:totalFrame
            seg1 = Xl{f}(seg);
            varianceL(f,i) = var(seg1);
        end
    end
    sumVarL = sum(varianceL,1);
    % weight Zl by variance
    Zl = zeros(size(regDecompMap{levels-1},1),size(regDecompMap{levels-1},2));
    for i=1:SegLen
        seg = find(regDecompMapLL == SegEle(i));
        for f = 1:totalFrame
            if sumVarL(i)==0
                Zl(seg) = Zl(seg) + Xl{f}(seg)* (1/totalFrame);
            else
                Zl(seg) = Zl(seg) + Xl{f}(seg)* (varianceL(f,i)/sumVarL(i));
            end
        end
    end;
elseif strcmp(lowpassfuse,'median')
    Zl = median(Zl,3);
else
    Zl = Zl/totalFrame;
end
% -------------------------------------------------------------------------
% Perform inverse DTCWT on fused wavelet coefs.
% -------------------------------------------------------------------------
zrest = dtwaveifm2(Zl,Zrg1,'near_sym_b','qshift_d',zgain);
zrestsmooth = dtwaveifm2(Zl,Zrg1,'near_sym_b','qshift_d');
toc
zrestReg = zrest;
tempM = regMap==1;
tempM = bwmorph(tempM,'dilate',2);
h = fspecial('average', [5 5]);
mask = imfilter(double(tempM),h);

prest = Nick_pixel_fuse(input,levels,[2 1.4 1 1]);
zrestReg = prest.*mask + (1-mask).*zrestReg;