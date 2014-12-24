function [overlay,intolay,map,intmap,gsurf0] = Segment_Rob(someImVol, stype, levels, selectDirection, marker, nomerge)
% [overlay,intolay,map,intmap] = Segment_Rob(someImVol,stype, levels);
% Combined morphological-spectral multimodal segmentation
% Inputs -  someImVol: registered image set
%           someHighCoefVol are the High pass wavelet coefficients for all the images.
%           stype: either 'Joint' or 'Unimode' segmentation
%
% Outputs - overlay,intoverlay: segmentation overlayed on mean image
%           map, intmap: integer segmentation maps
%
% "int" refers to the intermediate map, before spectral clustering
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                                    %%
%%          by Rob O'Callaghan                                                        %%
%%          edited by John Lewis Sept '03                                             %%
%%          edited by Eduardo Canga Feb '04                                           %%
%%          edited by Artur Loza March '05                                            %%
%%          edited by Pui May '12 - allow directional selection, allow desired marker %%
%%          University Of Bristol, UK                                                 %%
%%          Copyright 2003                                                            %%
%%                                                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% angles = [+15, -15, +75, -75, +45, -45];
if (nargin < 4)|| isempty(selectDirection)
    selectDirection = 1:6;
end
if nargin < 6
    nomerge = 0;
end
numband = length(selectDirection);

if max(someImVol(:)) > 1
    % Assume colormap is gray 255
    scale = 1;
else
    % otherwise pixel values in the range [0,1]
    scale = 255;
    someImVol = 255 * someImVol;
end

[N, M, P] = size(someImVol);

gsurf = zeros(N,M);
global globalmax;
globalmax = zeros(levels*numband,P);
tvol = [];

for k=1:P
    someIm = someImVol(:,:,k);
    [Yl,tvold] = dtwavexfm2(someIm,levels,'antonini','qshift_06');
    for n=1:levels;
        tvold{n} = abs(tvold{n}(:,:,[1 6 3 4 2 5]));  %rearrange for backward compatibility with older code
        temp  = squeeze(max(max(tvold{n},[],2),[],1))./(2.^n);
        globalmax((n-1)*numband+1:n*numband,k) = temp(selectDirection);
    end;
    
    tvold = texmedfiltd(tvold);
    % -----------------------------------------------------------------------
    % selective orientation
    if length(selectDirection)<6
        for n=1:levels;
            tvold{n} = tvold{n}(:,:,selectDirection);
        end
    end
    % -----------------------------------------------------------------------
    [gsurfpart(:,:,k),tvolpart] =segprotomm(tvold,double(someIm));
    tvol=cat(3,tvol,tvolpart);
end;

clear tvolpart tvold someIm someHighCoefVol;

if sum(stype) == 1 || strcmp(stype,'Joint')
    gsurf0 = sum(gsurfpart,3);
elseif sum(stype) == 2 || strcmp(stype,'Unimode')
    gsurf0 = gsurfpart;
end
clear gsurfpart
%someIm = uint8(mean(double(someImVol),3));%crappy pixel fusion just to see overlay result
%%
sed = strel('square',3); %shave off nobbly bits in watersheds (trust me, I know what I'm talking about)
for n = 1 : size(gsurf0,3)
    
    gsurf = gsurf0(:,:,n);
    
    if (nargin<5)||(isempty(marker))
        map_t = watershed(imhmin(gsurf,0.15*median(gsurf(:))));
        
    else
        gradBox = imimposemin(gsurf,marker);
        map_t = watershed(gradBox);
    end
    clear gsurf
    map2 = zeros(size(map_t));
    for k=1:max(map_t(:));
        map2(imclose(map_t==k,sed)) = k;
    end;
    intmap(:,:,n) = map2;
    
    if (nargin<5)||(isempty(marker))
        if nomerge
            map(:,:,n) = map2;
        else
            map(:,:,n) = mergestatcombmm(map2,tvol,double(someImVol));
        end
    else
        map(:,:,n) = map_t;
    end
    
end

for n=1:P
    overlay(:,:,n) = max(255*(map(:,:,min(end,n))==0),someImVol(:,:,n));
    intolay(:,:,n) = max(255*(intmap(:,:,min(end,n))==0),someImVol(:,:,n));
end

for i = 1:size(map,3)
    map(:,:,i) = remapmap(map(:,:,i), someImVol(:,:,i));
    intmap(:,:,i) = remapmap(intmap(:,:,i), someImVol(:,:,i));
end


overlay=overlay/scale;
intolay=intolay/scale;