function shift = affineshift2(avec,sizesh,mode);

% function sh1 = affineshift2(avec,sizesh,mode);
% Create shift vector from the affine vectors in avec(:,:,:).
% mode(4) determines how slices 7 to 12 of avec are used (if at all). 
% 
% sizesh(1:2) is the size required for the shift matrix.
% This version removes the restriction that sizesh(1:2) must be
% at least as big as size(avec,1:2).
% Note that sizesh must be even and an exact multiple or submutliple
% of size(avec) or size(avec)+1.
% It is assumed that xi goes from -1 to +1 across the full width of the
% image (ie to the outer edges of the pixels).  Hence if there are n pixels
% per row, the centres of these will be at xi = (2*k - n - 1)/n where k = 1:n.
% The avec sampling points are at the centre of each 2x2 group of pixels
% and these will then be at xi = (2*k - n)/n where k = 1:(n-1).
%
% Nick Kingsbury, Cambridge University, July 2004.

% Increase size of avec by 1 if dimension is odd.
sav = size(avec);
if rem(sav(1),2) == 1
   t = 1:sav(1);
   avec = 0.5 * (avec([1 t],:,:) + avec([t sav(1)],:,:));
end
if rem(sav(2),2) == 1
   t = 1:sav(2);
   avec = 0.5 * (avec(:,[1 t],:) + avec(:,[t sav(2)],:));
end

av = avec;
sav = size(av);
upfactor = sizesh(1) / sav(1);

if upfactor > 1 
    % Upsample the affine vectors avec.
    av = [];
    for k = 1:sav(3)
        av(:,:,k) = upsample(squeeze(avec(:,:,k)),round(upfactor));
    end
elseif upfactor < 1
    % Downsample the affine vectors avec by averaging each 2x2 group.
    while sav(1) > sizesh(1),
        t1 = 1:2:sav(1);
        t2 = 1:2:sav(2);
        av = 0.25*(av(t1,t2,:) + av(t1+1,t2,:) + av(t1,t2+1,:) + av(t1+1,t2+1,:));
        sav = size(av);
    end
end

% Create x and y indices, going from -1 to +1 across the image.
sc = size(av,1); sr = size(av,2);
yi = ((2*[1:sc]'-sc-1) / sc) * ones(1,sr);
xi = ones(sc,1) * ((2*[1:sr]-sr-1) / sr);
j = sqrt(-1);

% Calculate shift vectors at each location.
shift = av(:,:,2) + av(:,:,4).*yi + av(:,:,6).*xi + ...
   j*(av(:,:,1) + av(:,:,3).*yi + av(:,:,5).*xi);

% Add extra higher order terms depending on mode(4).
if mode(4) == 1,      % xy mode
   shift = shift + (av(:,:,8) + j*av(:,:,7)).*xi.*yi;
elseif mode(4) == 2,  % Linear zoom mode.
   shift = shift + (av(:,:,8) + j*av(:,:,7)).*xi.*yi + ...
      mode(5)*(av(:,:,7).*xi.*xi + j*av(:,:,8).*yi.*yi);
elseif mode(4) == 3,  % Linear zoom mode with separate x^2 and y^2 terms.
   shift = shift + (av(:,:,8) + j*av(:,:,7)).*xi.*yi + ...
      av(:,:,9).*xi.*xi + j*av(:,:,10).*yi.*yi;
end

return;

