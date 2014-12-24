function avgFrame = findRefFrame(input, refFrameType)

if ~strcmp(refFrameType,'average') && ~strcmp(refFrameType,'median')
    grad = zeros(1, size(input,3));
    for k = 1:size(input,3)
        if strcmp(refFrameType,'maxGradient')
            [fx,fy] = gradient(input(:,:,k));
            grad(k) = sum(abs(fx(:))+abs(fy(:)));
        else
            [~,Zh] = dtwavexfm2(input(:,:,k),4,'near_sym_b','qshift_d');
            grad(k) = sum(abs(Zh{1}(:))) + sum(abs(Zh{2}(:))) + sum(abs(Zh{3}(:))) + sum(abs(Zh{4}(:)));
        end
    end
    [~, ind] = max(grad);
    avgFrame = input(:,:,ind);
elseif strcmp(refFrameType,'median')
    avgFrame = median(input,3);
else
    avgFrame = sum(input,3)/size(input,3);
end