function [valcostSelect, indcostSelect] = findGoodFrame(avgFrame, input, rangei, rangej, rangek, numFrameRegis, sharpMethod, useBigAreaInfo, ...
     mseGain, gradGain,infoGain)


totalFrame = size(input,3);
numFrameRegis = min(numFrameRegis,totalFrame);

% compute mse and grad
% ref frame
refROI = avgFrame(rangei,rangej,:);
grad = zeros(1,totalFrame);
mse  = zeros(1,totalFrame)+realmax;
infoArea = zeros(1,totalFrame);
for f = rangek(1:numFrameRegis)
    curROI = input(rangei,rangej,f);
    % using gradient to find sharpness
    if strcmp(sharpMethod,'gradient')
        [fx,fy] = gradient(curROI);
        grad(f) = sum(abs(fx(:))+abs(fy(:)));
    else
        [~,Zh] = dtwavexfm2(curROI,4,'near_sym_b','qshift_d');
        grad(f) = sum(abs(Zh{1}(:))) + sum(abs(Zh{2}(:))) + sum(abs(Zh{3}(:))) + sum(abs(Zh{4}(:)));
    end
    
    if useBigAreaInfo
        infoArea(f) = sum(curROI(:)<otsuth);
    end
    mse(f)  = sum((refROI(:)-curROI(:)).^2)/length(rangei)/length(rangej);
end

valgrad = sort(grad, 'descend');
valmse = sort(mse);
gainGrad = 1/mean(valgrad(1:numFrameRegis));
gainMSE  = 1/mean(valmse(1:numFrameRegis));
% calculate cost for selecting frames
costSelect = mseGain*(1./(1+gainMSE*(mse))) + gradGain*(1 - 1./(1+gainGrad*(grad)));
if useBigAreaInfo
    valinfoArea = sort(infoArea, 'descend');
    infoGain = 1/mean(valinfoArea(1:numFrameRegis));
    costSelect = costSelect + infoGain*(1 - 1./(1+infoGain*infoArea));
end
[valcostSelect, indcostSelect] = sort(costSelect, 'descend');
