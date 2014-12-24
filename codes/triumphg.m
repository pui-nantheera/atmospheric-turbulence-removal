function [overlay,intolay,map,intmap] =triumphc(someIm);
levels=4;

N=size(someIm,1);
M=size(someIm,2);


%[C,S] = dtwavedec2(someIm,levels,'antonini','qshift_06');
[C,S] = dtwavedec2(someIm,levels,'near_sym_a','qshift_a');
tvold = cell(levels,1);
for n=1:levels;
    tvold{n} = abs(cwtband6(C,S,n));%./2^n;
end;
clear C S;

%%%%%%%%%%%%%%%%%%%%%%%
global globalmax;
globalmax = zeros(levels.*6,1);
for n=1:levels;
    globalmax((n-1)*6+1:n*6) = squeeze(max(max(tvold{n},[],2),[],1))./(2.^n);
end;
%%%%%%%%%%%%%%%%%%%%%%%

tvold = texmedfiltd(tvold);

[map,intmap]=segprotog(tvold,double(someIm));
edges = (map==0).*255;
overlay = max(uint8(edges),someIm);
edges = (intmap==0).*255;
intolay = max(uint8(edges),someIm);
%overlay = max(tribandedge,someIm);
return;