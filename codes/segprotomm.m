function[gradsurf,tvol]=segprotoc(tvolcell,inimage);
levels = length(tvolcell);
tgradim = texgradd(tvolcell);
tgradim = tgradim(1:size(inimage,1),1:size(inimage,2));

numband = size(tvolcell{1},3);

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
    for m=1:numband;
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
szmat = [size(tvolcell{1}(:,:,1)).*2,numband.*levels];
% tvolcell2 = texgreynorm(tvolcell2,greyimage);  %tvolcell2 was the original before erosion
% tvolcell2 = texrecede(tvolcell2);
% tvol = quickinterp(tvolcell2,szmat);
tvol = quickinterp(tvolcell,szmat);
tvol = tvol(1:size(inimage,1),1:size(inimage,2),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear tvolcell;

texen = gettexen(tvol);

greygradim = gradim(inimage,1.5);
gmed = median(greygradim(:));
greygradim(greygradim<0.5*gmed)=0;  % to kill off noise before we modify and potentially amplify it

%texfloor = max(median(texen(:)),3000);   % a little aribitrary with the absolute floor
%texfloor = 3000;
%texen = (max(texen,texfloor)./texfloor).^3; %to avoid divide by zero, and set a ceiling on greygrad enhancement
if gmed>0;
    greygradim = greygradim./gmed;
end;
%greygradim = mean(cat(3,cgradim,greygradim),3);


%tx2 = (texen/3-3)/2;
tx2 = (texen/2-7);
clear texen;
tx2 = exp(max(tx2,0));
greygradim = greygradim./tx2;  
clear tx2;
% ggnorm=median(greygradim(:));
% if(ggnorm);
%     greygradim = greygradim./ggnorm;
% end;

% remove anything below a threshold - this is not noise removal, but low
% amplitude basin removal for watershed, hence we clamp to the threshold not
% zero
tgradim(tgradim<(0.1*max(tgradim(:)))) = 0.1*max(tgradim(:)); 


%lowval = 0.75*median(double(greyimage(:)));
%greygradim = greygradim./filter2(ones(3)./9,max(greyimage,lowval));   %makes grad an indicator of relative gradient / difference
%greygradim = greygradim.*(2.*prod(size(greyimage))./sum(greyimage(:))); %normalise by mean grey level
%norm = median(temp(:));
% Now, clean up, canny style
% Using pixels over a  high threshold as the marker image, use morph.
% reconstruction to include lower, but connected gradient peaks
% marker = greygradim;
% marker(marker<(0.05)) = 0;
% greygradim=imreconstruct(marker,greygradim,4);
% % Now remove values beneath the low threshold
% greygradim(greygradim<0.01)=0.01;
%greygradim = greygradim./median(greygradim(:));
%absolutetexfac = max(min(texen(:))/3000,1);


gradsurf = tgradim.*4+greygradim;%./absolutetexfac;
clear tgradim greygradim;
gradsurf = filter2(ones(3),gradsurf);
