function [overlay,intolay,map,intmap] = Segment_Rob(someImVol, stype, levels);
% [overlay,intolay,map,intmap] = Segment_Rob(someImVol,stype, levels);
% Combined morphological-spectral multimodal segmentation
% Inputs -  someImVol: registered image set
%           someHighCoefVol are the High pass wavelet coefficients for all the images. 
%           stype: either 'joint' or 'unimode' segmentation
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
%%          University Of Bristol, UK                                                 %%
%%          Copyright 2003                                                            %%
%%                                                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=size(someImVol,1);
M=size(someImVol,2);
P=size(someImVol,3);
gsurf = zeros(N,M);

global globalmax;
globalmax = zeros(levels.*6,P);

tvol = [];

for k=1:P
    someIm = someImVol(:,:,k);
    [Yl,tvold] = dtwavexfm2(someIm,levels,'antonini','qshift_06');
    for n=1:levels;
        tvold{n} = abs(tvold{n}(:,:,[1 6 3 4 2 5]));  %rearrange for backward compatibility with older code
    end;
      
    %%%%%%%%%%%%%%%%%%%%%%%
    for n=1:levels;
        globalmax((n-1)*6+1:n*6,k) = squeeze(max(max(tvold{n},[],2),[],1))./(2.^n);
    end;
    %%%%%%%%%%%%%%%%%%%%%%%
    
    tvold = texmedfiltd(tvold);
    [gsurfpart(:,:,k),tvolpart] =segprotomm(tvold,double(someIm));
    tvol=cat(3,tvol,tvolpart);
end;
clear tvolpart tvold someIm someHighCoefVol;

if stype == 1 | strcmp(stype,'Joint')
    gsurf = sum(gsurfpart,3);
    %gradmed = median(gsurf(:));
    map = watershed(imhmin(gsurf,0.15*median(gsurf(:))));
    clear gsurf gradmed
    %shave off nobbly bits in watersheds (trust me, I know what I'm talking about)
    sed = strel('square',3);
    map2 = zeros(size(map));
    for k=1:max(map(:));
        map2(imclose(map==k,sed)) = k;
    end;

    intmap = map2;
    map = mergestatcombmm(map2,tvol,double(someImVol));

    %crappy pixel fusion just to see overlay result
    someIm = uint8(mean(double(someImVol),3));

    edges = (map==0).*1;
%    overlay = max(uint8(edges),someIm);
    for i = 1:P
        overlay(:,:,i) = max(double(edges),someImVol(:,:,i));
    end;
    edges = (intmap==0).*1;
    for i = 1:P
        intolay(:,:,i) = max(double(edges),someImVol(:,:,i));
    end;
%    intolay = max(uint8(edges),someIm);

    
elseif stype == 2 | strcmp(stype,'Unimode')
    mapvol = zeros(size(gsurfpart));
    for n=1:size(gsurfpart,3);
        thisgsurf = gsurfpart(:,:,n);
        map_t = watershed(imhmin(thisgsurf,0.15*median(thisgsurf(:))));
        %shave off nobbly bits in watersheds (trust me, I know what I'm talking about)
        sed = strel('square',3);
        map2 = zeros(size(map_t));
        for k=1:max(map_t(:));
            map2(imclose(map_t==k,sed)) = k;
        end;

        intmap(:,:,n) = map2;
        map(:,:,n) = mergestatcombmm(map2,tvol,double(someImVol));

        edges = (map(:,:,n)==0).*1;
        overlay(:,:,n) = max(double(edges),someImVol(:,:,n));
        edges = (intmap(:,:,n)==0).*1;
        intolay(:,:,n) = max(double(edges),someImVol(:,:,n));
        
    end;
    clear thisgsurf
end;
clear gsurfpart

for i = 1:size(map,3)
    map(:,:,i) = remapmap(map(:,:,i), someImVol(:,:,i));
    intmap(:,:,i) = remapmap(intmap(:,:,i), someImVol(:,:,i));
end;
return;