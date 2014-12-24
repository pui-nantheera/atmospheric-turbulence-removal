function tvolout = texmedfiltd(tvol);
tvolout = tvol;

numlevels = length(tvol);

%angles = [-45/2, 45/2, -3*45/2, 3*45/2, -45, 45];
angles = [+15, -15, +75, -75, +45, -45];
se1 = strel('square',3);
for m=1:6;
    for n=1:numlevels;   
        se_p=strel('line',7+2*n,angles(m));
        se_t=strel('line',7+2*n,angles(m)+90);
        med_hood_p = getnhood(se_p);
        med_hood_t = getnhood(se_t);
        order_p= ceil(nnz(med_hood_p)/2);
        order_t = ceil(nnz(med_hood_t)/2);
        
        tvolout{n}(:,:,m) = ordfilt2(tvolout{n}(:,:,m),order_t,med_hood_t,'symmetric');
        tvolout{n}(:,:,m) = ordfilt2(tvolout{n}(:,:,m),order_p,med_hood_p,'symmetric');
    end;
end;