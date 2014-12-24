% Post process denoising and sharpening

function outputImage = postDenoiseSharpen(inputImage,levels, option)

inputImage = double(inputImage);

% if # of levels of wavelet transform isn't defined.
% ---------------------------------------------
if nargin < 2
    % dimension
    [height,width] = size(inputImage);
    if min(height,width)<=128
        levels = 3;
    else
        levels = 4;
    end
end
if nargin < 3
    option = 1;
end

% DT-CWT
[Xl,Xh] = dtwavexfm2(inputImage,levels,'near_sym_b','qshift_d');

% new highpass
Zh = Xh;

% Noise variance estimation using robust median estimator..
tmp = real(Xh{1}(:,:,:));
Nsig = median(abs(tmp(:)))/0.6745;

% Set the windowsize and the corresponding filter
windowsize  = 7;
windowfilt = ones(1,windowsize)/windowsize;

% sharpening gain
zgain = [1 1 1 1];%[2.3 1.5 1.1 1];
% h = fspecial('average', [5 5]);
maskAreaTh = 0.25; % how many pixels to boost for sharpening

% for each level
for level = 1:levels
    if level > length(zgain)
        gain = 1;
    else
        gain = zgain(level);
    end
    h = fspecial('average', 2*[levels-level+1 levels-level+1]);
    for band = 1:6
        currentBand = Xh{level}(:,:,band);
        
        % -----------------------------------------------------------------
        % denoising
        % -----------------------------------------------------------------
        if level < levels
            % The corresponding noisy parent coefficients
            parentBand = Xh{level+1}(:,:,band);
            % Extend noisy parent matrix to make the matrix size the same as the coefficient matrix.
            parentBand  = expand(parentBand);
        else
            parentBand = Xh{level}(:,:,band);
        end
        
        % Signal variance estimation
        Wsig = conv2(windowfilt,windowfilt,(real(currentBand)).^2,'same');
        Ssig = sqrt(max(Wsig-Nsig.^2,eps));
        
        % Threshold value estimation
        T = sqrt(3)*Nsig^2./Ssig;
        
        % Bivariate Shrinkage
        if all(size(currentBand)==size(parentBand))
            denoiseBand = bishrink(currentBand,parentBand,T);
        else
            denoiseBand = currentBand;
        end
        
        % -----------------------------------------------------------------
        % sharpening
        % -----------------------------------------------------------------
        if option
            % threshold for cleaning the mask map
            meanDetails = mean(abs(currentBand(:)));
            stdDetails  = std(abs(currentBand(:)));
            maskArea = realmax;
            ratio = 0.5;
            while maskArea > maskAreaTh
                thresClean  = meanDetails + ratio*stdDetails;
                
                % create mask
                mask = abs(currentBand) > thresClean;
                mask = bwmorph(mask,'clean');
                mask = bwmorph(mask,'spur');
                mask = bwmorph(mask,'dilate');
                maskArea = sum(mask(:))/size(mask,1)/size(mask,2);
                ratio = ratio + 0.2;
            end
            mask = imfilter(double(mask),h);
            
            % finding gain - erf function is used
            if gain > 1
                gainMat = max(abs(currentBand(:)))*erf(gain*abs(currentBand)/max(abs(currentBand(:)))).*currentBand./abs(currentBand);
            else
                gainMat = currentBand;
            end
            
            % -----------------------------------------------------------------
            % combining
            % -----------------------------------------------------------------
            currentBand = mask.*gainMat + (1-mask).*denoiseBand;
            Zh{level}(:,:,band) = currentBand;
        else
            Zh{level}(:,:,band) = denoiseBand;
        end
    end
end

% inverse DT-CWT
outputImage = dtwaveifm2(Xl,Zh,'near_sym_b','qshift_d');