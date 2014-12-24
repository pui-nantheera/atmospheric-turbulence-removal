function Y = LabDecompProt(LabSeg, Levels, ImSize)

% Takes the labled segements and creates a labeling map which covers all
% the levels of the wavelet decomposition of the image
% The shape of the smaller regions is ensured, by ensuring that when a 
% must be chosen pixel from one of two regions, the pixel will always be
% taken from the smaller of the two.

Y = cell(Levels, 1);
levsize = ceil(ImSize./2);
for a = 1:Levels
    Y(a, :) = { zeros(levsize(1), levsize(2), 1) };    
    levsize = ceil(levsize./2);
end;

Y{1} = LabSeg(1:2:end, 1:2:end);
teven = LabSeg(2:2:end, 2:2:end);
sizeY = size(Y{1});
sizeeven = size(teven);
%% deal with level a-1 not being divisible by 2
if (sizeY(1) ~= sizeeven(1))
    teven(sizeY(1), :) = Y{1}(sizeY(1), 1:sizeeven(2));
end;
if (sizeY(2) ~= sizeeven(2))
    teven(:, sizeY(2)) = Y{1}(:, sizeY(2));
end;    

Y{1}(Y{1}==0) = teven(Y{1}==0);



ImSize = ceil(ImSize./2);

%% Remove any more boundary pixels - deal with the corners first
if (Y{1}(1,1)==0)
    Y{1}(1,1) = Y{1}(2,2);
end;
if (Y{1}(ImSize(1),1)==0)
    Y{1}(ImSize(1),1) = Y{1}(ImSize(1)-1,2);
end;
if (Y{1}(1,ImSize(2))==0)
    Y{1}(1,ImSize(2)) = Y{1}(2,ImSize(2)-1);
end;
if (Y{1}(ImSize(1),ImSize(2))==0)
    Y{1}(ImSize(1),ImSize(2)) = Y{1}(ImSize(1)-1,ImSize(2)-1);
end;

%Check for any remaining blocks
while any(any(Y{1}==0))
    Y{1}(find(Y{1}==0)) =  Y{1}(find(Y{1}==0)-1);
end;


%imview(Y{1}, []);

teven(teven==0) = Y{1}(teven == 0);
dif = (Y{1} ~= teven);
%imview(dif);

lab = 1;
while (any(any(LabSeg==lab)))
    seglen(lab) = length(find(LabSeg==lab));
    lab = lab +1;
end;

if seglen(Y{1}(dif==1)) > seglen(teven(dif==1))
    Y{1}(dif==1) = teven(dif==1);
end;
%imview(Y{1}, []);

%figure, imshow(Y{1}, []);
for a = 2:Levels
    clear teven; clear dif;
    Y{a} = Y{a-1}(1:2:end, 1:2:end);
    %imview(Y{a}, []);
    teven = Y{a-1}(2:2:end, 2:2:end);
    sizeY = size(Y{a});
    sizeeven = size(teven);
    %% deal with level a-1 not being divisible by 2
    if (sizeY(1) ~= sizeeven(1))
        teven(sizeY(1), :) = Y{a}(sizeY(1), 1:sizeeven(2));
    end;
    if (sizeY(2) ~= sizeeven(2))
        teven(:, sizeY(2)) = Y{a}(:, sizeY(2));
    end;
    dif = (Y{a} ~= teven);
    mask = zeros(size(teven));
    mask(dif==1) = seglen(Y{a}(dif==1)) > seglen(teven(dif==1));
    Y{a}(mask==1) = teven(mask==1);
    
    %imview(Y{a}, []);
    
%    Y{a}(Y{a}==0) = temp(Y{a}==0);
%    clear temp;

%    figure, imshow(Y{a}, []);
end;






        
  
        
    
        
        
    

