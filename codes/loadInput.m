function [input, inputU, inputV] = loadInput(folderName,extfile, totalFrame, resizeRatio)

files = dir([folderName,'*.',extfile]); 
if length(totalFrame)==1
    totalFrame = min(totalFrame,length(files));
    rangek = 1:totalFrame;
else
    rangek = unique(min(totalFrame,length(files)));
    totalFrame = length(rangek);
end

for k = 1:totalFrame
    img = imread([folderName,files(rangek(k)).name]);
    
    if resizeRatio~=1
        img = imresize(img,resizeRatio);
    end
    
    if k==1
        % get dimension
        [height, width] = size(img(:,:,1));
        input = zeros(height, width, totalFrame);
        if size(img,3)==3
            inputU = zeros(height, width, totalFrame);
            inputV = zeros(height, width, totalFrame);
        else
            inputU = [];
            inputV = [];
        end
    end
    
    if size(img,3)==3
        img = double(rgb2ycbcr(img));
        inputU(:,:,k) = img(:,:,2);
        inputV(:,:,k) = img(:,:,3);
    end
    input(:,:,k) = double(img(:,:,1));
end