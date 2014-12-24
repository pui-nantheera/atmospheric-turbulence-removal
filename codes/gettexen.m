function outim = gettexen(tvol);

levels = size(tvol,3)./6;
% se = strel('square',3);
% % for n=1:levels;
% %     for m=1:6;
% %         tvolcell{n}(:,:,m)=imerode(tvolcell{n}(:,:,m),se);
% %     end;
% % end;
% 
% szmat = [size(tvolcell{1}(:,:,1)).*2,levels*6];
% tvol = quickinterp(tvolcell,szmat);
outim = zeros(size(tvol,1),size(tvol,2));

for m=1:6;
    %outim = outim + prod(tvol(:,:,m:6:end),3).^(1/levels); % on the basis that edges appear at all scales ...
    outim = outim + sum(tvol(:,:,m:6:end),3)./levels;  %this form is a relic of the above line, as I want the values to end up comparable
end;