function [overlay,intolay,map,intmap, gsurf] = runseg(Im, t1,t2,hminfactor)
%% Segments Image (Im) as a colour image if it is colour RGB image or in
%% gray scale only if it is a gray image (one plane only).  All images
%% should be in the range 0 - 255

Levels = 4;

Im_RGB = Im;
%Im = rgb2ycbcr(Im_RGB);
ImSize = size(Im);
highcoef = cell(Levels, 1); 
%highcoef = cell(Levels, ImSize(3)); %% removed for single im segmentation
% wavelet transform the images
for a = 1%:ImSize(3)  
    [lowcoef(:,:,a),highcoef(:,a)] = dtwavexfm2(double(Im(:,:,a)),Levels,'antonini','qshift_06');
end;  
%segment the images
[overlay,intolay,map,intmap, gsurf]=cmssegmm(Im, highcoef, 'joint', Levels,t1,t2,hminfactor);
for a = 1%:ImSize(3)  
    overlay(:,:,a) = max(double(Im_RGB(:,:,a)), (map==0)*255);
    intolay(:,:,a) = max(double(Im_RGB(:,:,a)), (intmap==0)*255);
end;  

% imwrite(uint8(overlay), 'D:\work\prog\jpg2ksynth_seg_opt\working/16 overlay.png', 'png');
% imwrite(uint8(intolay), 'D:\work\prog\jpg2ksynth_seg_opt\working/17 pre speclust overlay.png', 'png');
% imwrite(uint8(intmap), 'D:\work\prog\jpg2ksynth_seg_opt\working/19 pre sc map.png', 'png');
% imwrite(uint8(map), 'D:\work\prog\jpg2ksynth_seg_opt\working/18 map.png', 'png');

