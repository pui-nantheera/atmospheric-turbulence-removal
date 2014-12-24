function ImgMap = reIndexMap(Map)

% reindex the image segments

SegEle = unique(Map);
SegLen = length(SegEle);

ImgMap = zeros(size(Map));

for i=1:1:SegLen
  ImgMap(find(Map == SegEle(i)))= i;
end