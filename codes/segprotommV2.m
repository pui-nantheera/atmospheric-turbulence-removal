function[gradsurf,tvol]=segprotoc(tvolcell,inimage);
doI =any(inimage(:));
doT = ~isempty(tvolcell);

if doT;
levels = length(tvolcell);
tgradim = texgradd(tvolcell);
tgradim = tgradim(1:size(inimage,1),1:size(inimage,2));



% tgnorm=max(tgradim(:));
% if(tgnorm);
%     tgradim = tgradim./tgnorm;
%     %teng = energyim(tgradim);
%     %tgradim = tgradim./teng;
% end;

% erode back texture response so as not to stifle valid gradients at edges
% of texture regions

se = strel('square',3);
for n=1:levels;
    for m=1:6;
        tvolcell{n}(:,:,m)=imerode(tvolcell{n}(:,:,m),se)./2^n;
    end;
end;

%upsample the texture subbands
% szmat = [size(greyimage),levels*6];
% tvol = quickinterp(tvolcell,szmat);


% texen = sum(tvol,3).^2;
% greyen = mean(double(greyimage(:)).^2);
% texmeanen = mean(texen(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global globalvolcell;
%upsample the texture subbands
szmat = [size(tvolcell{1}(:,:,1)).*2,6.*levels];
% tvolcell2 = texgreynorm(tvolcell2,greyimage);  %tvolcell2 was the original before erosion
% tvolcell2 = texrecede(tvolcell2);
% tvol = quickinterp(tvolcell2,szmat);
tvol = quickinterp(tvolcell,szmat);
tvol = tvol(1:size(inimage,1),1:size(inimage,2),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear tvolcell;

texen = gettexen(tvol);
tgradim(tgradim<(0.1*max(tgradim(:)))) = 0.1*max(tgradim(:));
else
    tvol = [];
end;

if doI;
greygradim = gradim(inimage,1.5);
gmed = median(greygradim(:));
greygradim(greygradim<0.5*gmed)=0;  % to kill off noise before we modify and potentially amplify it
if gmed>0;
    greygradim = greygradim./gmed;
end;
end;

if (doT&doI);
%tx2 = (texen/3-3)/2;
tx2 = (texen/2-7);
clear texen;
tx2 = exp(max(tx2,0));
greygradim = greygradim./tx2;  
clear tx2;
end;

if (doT&doI);
    gradsurf = tgradim.*4+greygradim;%./absolutetexfac;
elseif doT;
    gradsurf = tgradim.*4;
else %(doI);
    gradsurf = greygradim;
end;


clear tgradim greygradim;
gradsurf = filter2(ones(3),gradsurf);
