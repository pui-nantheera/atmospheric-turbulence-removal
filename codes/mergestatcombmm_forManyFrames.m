function outmap = mergestatcombmm_forManyFrames(map,numLoop,gvol)
%This version is now based on recursive spectral segmentation
% texdegree = texen./max(texen(:));

binsperdim=64;

origmap = map;
regions=nonzeros(unique(map(:)));
numreg = length(regions);
adjacemat = zeros(numreg);
dists1 = adjacemat;
dists2=dists1;
se1 = strel('line',6,0);
se2 = strel('line',6,90);
se3 = strel('square',3);
se4 = strel('square',5);

%generate eroded region map, which is only used to train the classes for
%region difference calculation
erodemap = zeros(size(map)+2);
map4erode = erodemap;
map4erode(2:end-1,2:end-1) = map;  %need to pad coz otherwise regions on periphery will erode towards edge!
%% for WMV
%meanarray = zeros(numreg,size(featvol,2));
%stdarray = meanarray;
%%
for k=1:numreg;
    template=(map4erode==k);
    prevtemplate=template;
    numerodes=0;
    while((nnz(template)>(6*binsperdim))&&(numerodes<7)); 
        prevtemplate=template;
        template=imerode(template,se4);
        numerodes=numerodes+1;
    end;
    erodemap(prevtemplate)=k;
    %%% for WMV
%     reg = featvol(prevtemplate(:),:);
%     meanarray(k,:) = mean(reg);
%     stdarray(k,:) = std(reg);
    %%%
end;
clear map4erode;
erodemap = erodemap(2:end-1,2:end-1);
%erodemap =map; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*******
%%% for WMV
% mudev = std(meanarray);
% sigmadev = std(stdarray);
%%%
keepinds = find(erodemap);
maxvec = [];
for k = 1:numLoop
    eval(['load tvol',num2str(k)]);
    featvolt = reshape(tvol,size(tvol,1)*size(tvol,2),size(tvol,3));
    lenfeatvolt = size(featvolt,2);
    eval(['save featvolt',num2str(k),' featvolt lenfeatvolt']);
    
    maxvec = [maxvec, max(featvolt,[],1)];
    clear featvolt tvol
end

%featvolt = reshape(tvol,size(tvol,1)*size(tvol,2),size(tvol,3));
featvolg = reshape(gvol,size(gvol,1)*size(gvol,2),size(gvol,3));
%featvol(:,end-1) = featvol(:,end-1) -min(featvol(:,end-1));
%featvol(:,end) = featvol(:,end) - min(featvol(:,end));
clear tvol gvol;



%Quantize texture feature set   (involves scaling into uniform range)
%make the bins equal width, over the range of the image

% tsize = size(featvolt);
% featvolt = featvolt(keepinds,:);

% maxvec = max(featvolt,[],1);
global globalmax;
maxvec = max(globalmax(:).'./2,maxvec);

binwidths = maxvec./binsperdim;
%featvolt = floor(featvolt./repmat(binwidths,size(featvolt,1),1))+1;
for l = 1:numLoop
    eval(['load featvolt',num2str(l)]);
    if l==1
        rangel = 1:lenfeatvolt;
    else
        rangel = rangel(end) + [1:lenfeatvolt];
    end
    for k=1:size(featvolt,1);
        featvolt(k,:) = floor(featvolt(k,:)./binwidths(rangel))+1;
    end;
    featvolt(featvolt>binsperdim)=binsperdim;  %just to force maxima into the bins
    eval(['save featvolt',num2str(l),' featvolt lenfeatvolt']);
end
% tempvol = zeros(tsize);
% tempvol(keepinds,:) = featvolt;
% featvolt = tempvol;
% clear tempvol;

%Quantize intensity feature set   (involves scaling into uniform range)
%make the bins equal width, over the range of the image

tempvol = zeros(size(featvolg));
featvolg = featvolg(keepinds,:);

maxvec = 255.*ones(1,size(featvolg,2));

binwidths = maxvec./binsperdim;
featvolg = floor(featvolg./repmat(binwidths,size(featvolg,1),1))+1;
featvolg(featvolg>binsperdim)=binsperdim;  %just to force maxima into the bins

tempvol(keepinds,:) = featvolg;
featvolg = tempvol;
clear tempvol;


%Find Adjacency
for k=1:numreg;
    %overlap = imdilate((map==regions(k)),strel('line',2,90)); %for no watershed pixels
    %overlap = ovelap | imdilate((map==regions(k)),strel('line',2,0));
    overlap = imdilate((map==regions(k)),se1); %for watershed pixels
    overlap = overlap|imdilate((map==regions(k)),se2);
    map(map==k)=0;
    %adjvec = nonzeros(overlap.*map);
    adjvec = map(overlap);
    adjvecun=unique(adjvec(adjvec>0));
    for i=1:size(adjvecun);
        %adjacemat(adjvecun(i),k)=1;
        adjacemat(adjvecun(i),k)=nnz(adjvec==adjvecun(i));
    end;
end;
map=origmap;
clear origmap;
adjacemat = adjacemat + adjacemat';

% %Find second order adjacency
% diaginds = logical(eye(size(adjacemat)));
% adjacemat(diaginds) =1;
% adjacemat = logical(adjacemat);
% newadjacemat = zeros(size(adjacemat));
% for k=1:numreg;
%     newadjacemat(k,:)=max(adjacemat(adjacemat(k,:),:),[],1);
% end;
% adjacemat = newadjacemat;
% adjacemat = max(adjacemat,adjacemat');
% adjacemat(diaginds)=0;

[xpos,ypos] = find(triu(adjacemat));
%[xpos,ypos] = find(triu(ones(size(adjacemat)))-eye(size(adjacemat))); %full distance matrix
for k=1:length(xpos);
    template1=(erodemap==xpos(k));
    template2=(erodemap==ypos(k));
    %template1f=(origmap==xpos(k));
    %template2f=(origmap==ypos(k));
    reg1t = [];
    reg2t = [];
    for l = 1:numLoop
        eval(['load featvolt',num2str(l)]);
        reg1t = cat(2,reg1t,featvolt(template1(:),:));
        reg2t = cat(2,reg2t,featvolt(template2(:),:));
    end
    reg1g = featvolg(template1(:),:);
    reg2g = featvolg(template2(:),:);
%     texval1 = texdegree(template1);
%     texval1 = sum(texval1)./length(texval1);
%     texval2 = texdegree(template2);
%     texval2 = sum(texval2)./length(texval2);
%     texval = max(texval1,texval2);

%     if((xpos(k)==7)&(ypos(k)==9));
%         rob=1;
%     end;
    %dists1(xpos(k),ypos(k))=distcalcb(reg1,reg2);
    %dists1(xpos(k),ypos(k))=distcalch(reg1,reg2,binwidths);
    %dists1(xpos(k),ypos(k))=log(1+distcalcmh(reg1,reg2,binsperdim));
    %dists1(xpos(k),ypos(k))=distcalcksq(reg1,reg2,binsperdim);
    [texdiff,greydiff]=distcalcmh(reg1t,reg2t,reg1g,reg2g,binsperdim);
    dists1(xpos(k),ypos(k)) = max([greydiff,texdiff]);%.*texval);
    %dists1(xpos(k),ypos(k))=emdtruncdist(reg1,reg2);
    
    %%for wmv
%     deltam = abs(meanarray(xpos(k),:)-meanarray(ypos(k),:));
%     deltas = abs(stdarray(xpos(k),:)-stdarray(ypos(k),:));
%     dists1(xpos(k),ypos(k))= sum(deltam./mudev+deltas./sigmadev);
    %%

%    brokenreg = (template1f|template2f);
 %   superreg = imclose(brokenreg,se3);

  %  edgetemplate = superreg&(~brokenreg);
    
%     if(any(edgetemplate(:)));
%         dists2(xpos(k),ypos(k))=max(gradsurf(edgetemplate));
%     else
%         panic=1e6;
%     end;
end;
clear featvolg featvolg erodemap;
%max1=max(dists1(:));
%max2=max(dists2(:));
%dists2 = (dists2./max2);
%dists2(dists2==0)=max((dists2(:)));  % takes account of non-adjacent comparisons with no "edge"
%dists2(logical(eye(size(dists2))))=0; % zero the diagonals again
%dists1 = (dists1./max1);
dists1 = dists1+dists1.';
%dists2 = dists2+dists2.'; 
dists = dists1;%+dists2;
origdists = dists;


%the business end of the segmentation
aff = dists;
aff(aff==0) = 1;
aff = 1-aff;
% aff(aff==0) = inf;
% aff= aff-min(aff(:));
% sigmaval = mean(aff(aff<inf));
% aff = exp(-(aff./sigmaval).^2);
aff(logical(eye(size(aff))))=1; %try self-affinity =1

%newmap = specsplit4(aff,adjacemat,-ones(length(aff),1),map,0);
newmap = specsplit5(aff.*adjacemat,adjacemat,map,0);
%newmap = specsplit5(aff,adjacemat,0.7,map,0);

list = unique(nonzeros(newmap));
outmap = zeros(size(newmap));
for k=1:length(list);
    outmap(newmap==list(k))=k;
end;

for k=1:max(outmap(:));
    outmap(imclose(outmap==k,se3))=k;
end;
return;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [texturediff,intensediff] = distcalcmh(r1t,r2t,r1g,r2g,binsperdim); %marginal histograms with either EMD or L1 measures


% texture first
r1 = r1t; r2 = r2t;
N=size(r1,1);
M=size(r2,1);
numdims = size(r1,2);

histObj1 = zeros(binsperdim,numdims);
histObj2 = histObj1;

for n=1:N;
    for k=1:numdims;
        histObj1(r1(n,k),k) = histObj1(r1(n,k),k) +1;
    end;
end;
for n=1:M;
    for k=1:numdims;
        histObj2(r2(n,k),k) = histObj2(r2(n,k),k) +1;
    end;
end;
histObj1 = histObj1./N;
histObj2 = histObj2./M;
% 1-D EMD dist now

histObj1 = cumsum(histObj1,1);
histObj2 = cumsum(histObj2,1);
%histObj1(:,1:end-1) = cumsum(histObj1(:,1:end-1),1);
%histObj2(:,1:end-1) = cumsum(histObj2(:,1:end-1),1);
histObj1 = abs(histObj1-histObj2);
dvec = sum(histObj1,1)./(binsperdim-1);
texturediff = min(max(dvec).*2,1);

% now greyscale
r1 = r1g; r2 = r2g;
N=size(r1,1);
M=size(r2,1);
numdims = size(r1,2);

histObj1 = zeros(binsperdim,numdims);
histObj2 = histObj1;

for n=1:N;
    for k=1:numdims;
        histObj1(r1(n,k),k) = histObj1(r1(n,k),k) +1;
    end;
end;
for n=1:M;
    for k=1:numdims;
        histObj2(r2(n,k),k) = histObj2(r2(n,k),k) +1;
    end;
end;
histObj1 = histObj1./N;
histObj2 = histObj2./M;
% 1-D EMD dist now

histObj1 = cumsum(histObj1,1);
histObj2 = cumsum(histObj2,1);
%histObj1(:,1:end-1) = cumsum(histObj1(:,1:end-1),1);
%histObj2(:,1:end-1) = cumsum(histObj2(:,1:end-1),1);
histObj1 = abs(histObj1-histObj2);
dvec = sum(histObj1,1)./(binsperdim-1);
intensediff = min(max(dvec).*2,1);
return;