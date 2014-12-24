function outmap = specsplit5(W,A,map,level);
%now using new formulation for maximised mean association
%and with propagation of external edges for decision making
level =level+1;
numreg=size(W,1);
thresh1 = 0.96;%0.94;%0.85;%
thresh2 = 0.65;%0.65;%0.74;%
%thresh = 1.2;
W(logical(eye(size(W))))=0;

if (numreg>2);  %non-trivial cut computation
    
    d=sum(W,2);
    
    D = diag(d);
    S = diag(sum(A,2));
    
    Nr = (D-W);
    Dr = (S-A);
    [V,evals]=eig(Nr,Dr); %seems faster than using eigs to get just 2 eigenvectors
    [minval,minind] = min(abs(diag(evals)));
    if( (max(V(:,minind))-min(V(:,minind)))<1e-3 );
        evals(minind,minind) = inf;
        [minval,minind] = min(abs(diag(evals)));
    end;
    [Y,inds] = sort(V(:,minind));
    Wsort = W(inds,inds);
    Asort = A(inds,inds);
    clear A W D S Nr Dr V evals;
    [crit,affval,cutpoint] = choosecut(Wsort,Asort);
    carryon = (crit<thresh2)|((crit./affval)<thresh1);  %can get divide by zero errors !!!
    
    finished = ~carryon;
    if(~finished);
        r1 = inds(1:cutpoint);
        r2 = inds(cutpoint+1:end);
        
        aff1 = Wsort(1:cutpoint,1:cutpoint);
        adj1 = Asort(1:cutpoint,1:cutpoint);
        aff2 = Wsort(cutpoint+1:end,cutpoint+1:end);
        adj2 = Asort(cutpoint+1:end,cutpoint+1:end);
        clear Asort Wsort;
       
        map1 = zeros(size(map));
        map2=map1;
        for k=1:length(r1);
            map1(map==r1(k))=k;
        end;
        for k=1:length(r2);
            map2(map==r2(k))=k;
        end;
        clear map;
        omap1=specsplit5(aff1,adj1,map1,level);
        omap2=specsplit5(aff2,adj2,map2,level);
        omap1(omap1>0) = omap1(omap1>0) + 2^level;
        outmap = omap1+omap2;
    else
        outmap = map;
        outmap(outmap>0)=1;
    end;

elseif (numreg==2);  %best cut is trivial if only two regions
    if (W(2,1)>thresh2);   %slightly arbitrary choice of when to merge two single regions!!
        outmap = map;
        outmap(outmap>0)=1;
    else
        outmap = map;
        outmap(outmap==1) = 1;
        outmap(outmap==2) = 1+2^level;
    end;
    
else  %Have been fed a single region by previous level (could avoid such calls altogether really)
    outmap = map;
    outmap(outmap>0)=1;
        
end;

return;

function [breakval,mval,index] = choosecut(affreorder,adjace);
% choose based on mean between 
N=size(affreorder,1)-1;
mmaffvals = zeros(N,1);

affreorder(logical(eye(size(affreorder))))=0;
for k=1:N;
    crossaff=nonzeros(affreorder(1:k,k+1:end));
   
    A1 = sum(sum(affreorder(1:k,1:k)));
    A2 = sum(sum(affreorder(k+1:end,k+1:end)));
    C = sum(sum(affreorder(1:k,k+1:end)));
    N1 = sum(sum(adjace(1:k,1:k)));
    N2 = sum(sum(adjace(k+1:end,k+1:end)));
    NC = sum(sum(adjace(1:k,k+1:end)));
    if A1>0;
        A1= A1./N1;
    end;
    if A2>0;
        A2 = A2./N2;
    end;
    if C>0;
        C = C./NC;
    end;
    mmaff = max(A1,A2);
    Cvals(k) = C;
    mmaffvals(k) = mmaff;

end;
%[breakval,index] = min(cutvals);
[breakval,index] = min(Cvals);
mval = mmaffvals(index);
if mval==0;
    nuts=2;
end;
return;