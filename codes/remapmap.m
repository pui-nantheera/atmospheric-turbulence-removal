function outmap = remapmap(inmap,inim);
%NB if you provide a reference image in RGB format, the current version of
%this function maps the greylevel to correspond to the red channel only!
numreg = max(inmap(:));
meanvec = zeros(numreg,1);
outmap = zeros(size(inmap));

for n=1:numreg;
    meanvec(n) = mean(inim(inmap==n));
end;
[vals,ind] = sort(meanvec);
for n=1:numreg;
    outmap(inmap==ind(n)) = n;
end;