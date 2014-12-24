function [overlay,intolay,map,intmap,gsurf] =cmssegmm(someImVol, someHighCoefVol, type, levels,t1,t2,hminfactor);
% [overlay,intolay,map,intmap] =cmssegmm(someImVol,type);
% Combined morphological-spectral multimodal segmentation
% Inputs -  someImVol: registered image set
%           someHighCoefVol are the High pass wavelet coefficients for all the images. 
%           type: either 'joint' or 'unimode' segmentation
%
% Outputs - overlay,intoverlay: segmentation overlayed on mean image
%           map, intmap: integer segmentation maps
%
% "int" refers to the intermediate map, before spectral clustering



N=size(someImVol,1);
M=size(someImVol,2);
gsurf = zeros(N,M);

global globalmax;
globalmax = zeros(levels.*6,size(someImVol,3));

tvol = [];

for k=1:size(someImVol,3);
    someIm = someImVol(:,:,k);
    tvold = someHighCoefVol(:,k);
    
    for n=1:levels;
        tvold{n} = abs(tvold{n}(:,:,[1 6 3 4 2 5]));  %rearrange for backward compatibility with older code
        globalmax((n-1)*6+1:n*6,k) = squeeze(max(max(tvold{n},[],2),[],1))./(2.^n);
    end;

    tvold = texmedfiltd(tvold);
    
    if k <=1
        col = 0;
    else 
        col = 1;
    end;
    
    [gsurfpart(:,:,k),tvolpart] =segprotomm(tvold,double(someIm), col);
    tvol=cat(3,tvol,tvolpart);
end;
clear tvolpart tvold someIm someHighCoefVol;

    gsurf = sum(gsurfpart,3);
    gradmed = median(gsurf(:));
    map = watershed(imhmin(gsurf,hminfactor*gradmed));
    clear gradmed
    %shave off nobbly bits in watersheds (trust me, I know what I'm talking about)
    sed = strel('square',3);
    map2 = zeros(size(map));
    for k=1:max(map(:));
        map2(imclose(map==k,sed)) = k;
    end;

    intmap = map2;
    map = mergestatcombmm(map2,tvol,double(someImVol),t1,t2);

    %crappy pixel fusion just to see overlay result
    someIm = uint8(mean(double(someImVol),3));

    edges = (map==0).*255;
    overlay = max(uint8(edges),someIm);
    edges = (intmap==0).*255;
    intolay = max(uint8(edges),someIm);


clear gsurfpart


return;