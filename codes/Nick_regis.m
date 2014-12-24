% turb_correct.m
%
% Register all the frames in the sequence of tiff frames, stored in directory dirname. 
% Jitter is assumed to be caused by atmospheric turbulence.
%
% Nick Kingsbury, Cambridge University, May-Aug 2011.
%
% Apply for colour video and update ref frame - Pui, Bristol University, Oct 2011


function [xrest,xrestU,xrestV] = Nick_regis(input, inputU, inputV, levels, numToAvg, update)


if nargin < 6
    update = 0;
end
if nargin < 5
    numToAvg = size(input,3);
end
if nargin < 4
    levels = 5;
end

[m,n,N] = size(input);

% reference frame
if length(numToAvg) == 1
    xref = sum(input(:,:,1:numToAvg),3)/numToAvg; %sum(input,3)/size(input,3);
    mask = 1;
else
    if size(numToAvg,1) ~= size(input,1)
        blockSize(1) =  size(input,1)/size(numToAvg,1);
        blockSize(2) =  size(input,2)/size(numToAvg,2);
        mask = zeros(size(input));
        for i = 1:blockSize(1):size(input,1)
            for j = 1:blockSize(2):size(input,2)
                if sum(numToAvg(ceil(i/blockSize(1)),ceil(j/blockSize(2)),:))==0
                    numToAvg(ceil(i/blockSize(1)),ceil(j/blockSize(2)),:) = 1;
                end
                mask(i-1 + (1:blockSize(1)), j-1 + (1:blockSize(2)), :) = repmat(numToAvg(ceil(i/blockSize(1)),ceil(j/blockSize(2)),:), [blockSize 1]);
            end
        end
    else
        mask = numToAvg;
        blockSize = [32 32];
    end
    curMask = mask;
    curMask = curMask + circshift(curMask,[blockSize(1)/4 0 0]);
    curMask = curMask + circshift(curMask,[-blockSize(1)/4 0 0]);
    curMask = curMask + circshift(curMask,[0 blockSize(2)/4 0]);
    curMask = curMask + circshift(curMask,[0 -blockSize(2)/4 0]);
    curMask = min(curMask,1);
    mask = imfilter(curMask, fspecial('average',[blockSize/4]));
    if size(numToAvg,1) ~= size(input,1)
        xref = sum(input.*mask,3)./sum(mask,3);
    else
        xref = sum(input.*numToAvg,3)./sum(numToAvg,3);
    end
end

% =========================================================================
% resize by adding padding
[height, width] = size(xref);
bigsize = 2.^(ceil(log2([height width])));
padaddy = 0; padaddx = 0;
addliney = 0; addlinex = 0;
if height~=bigsize(1)
    padaddy = floor((bigsize(1) - height)/2);
    xref = padarray(xref, [padaddy 0], 'symmetric');
    [height, width] = size(xref);
    if height~=bigsize(1)
        xref(end+1,:)  = xref(end,:);
        addliney = 1;
    end
end
if width~=bigsize(2)
    padaddx = floor((bigsize(2) - width)/2);
    xref = padarray(xref, [0 padaddx], 'symmetric');
    [height, width] = size(xref);
    if width~=bigsize(2)
        xref(:,end+1)  = xref(:,end);
        addlinex = 1;
    end
end
% =========================================================================

% CWT of xref
[Zl,Zh,Zsc] = dtwavexfm2(xref,levels,'near_sym_b','qshift_d');

% Define default mode parameters for estshhemm22().
if exist('mode')~=1, mode=[0 0 0 0 0]; end

ff = 1:N;  % Frames to register.
% ff = 1:5; % Reduce the number of frames if required.

xrest = input; % Prepare to save registered frames in xrest.
zref = zeros(m,n); % Prepare to accumulate registered frames in zref.

% Set up constants for registration using estshhemm22().
omega = -[1 3]*pi/2.15; % Centre frequencies of 1D lowpass and hipass filters.
levelsel = [4 4; 4 3; 4 3; 4 2; 4 2]; % DTCWT levels used by each iteration.
if levels < 4
    levelsel = [3 3; 3 3; 3 2; 3 2]; 
end
nit = size(levelsel,1); % No. of iterations.
avlevel = 3; % Resolution of affine blocks (level 3 = 8x8 blocks).
sizeqfilt = [4 2 1 1 1]; % Size of smoothing filter for Q matrices at each iteration.

% Register frames ff to xref.
for f = ff,
    
    x = double(input(:,:,f));
    if ~isempty(inputU)
        xU = double(inputU(:,:,f));
        xV = double(inputV(:,:,f));
    end
    
    % =========================================================================
    % resize by adding padding
    if (padaddy~=0) || (padaddx~=0)
        x = padarray(x, [padaddy padaddx], 'symmetric');
        x(end+addliney,:) = x(end,:);
        x(:,end+addlinex) = x(:,end);
        if ~isempty(inputU)
            xU = padarray(xU, [padaddy padaddx], 'symmetric');
            xV = padarray(xV, [padaddy padaddx], 'symmetric');
            xU(end+addliney,:) = xU(end,:);
            xV(end+addliney,:) = xV(end,:);
            xU(:,end+addlinex) = xU(:,end);
            xV(:,end+addlinex) = xV(:,end);
        end
    end
    % =========================================================================
    
    % figure; draw(x,[],posn); settitle('x image'); pause(0.1)
    ttl = sprintf('%d frame, registered',f);
    % figure; draw(x,[],posn); settitle(ttl); pause(0.1); fref = gcf
    
    % DTCWT of x
    [Xl,Xh] = dtwavexfm2(x,levels,'near_sym_b','qshift_d');
    if ~isempty(inputU)
        [XlU,XhU] = dtwavexfm2(xU,levels,'near_sym_b','qshift_d');
        [XlV,XhV] = dtwavexfm2(xV,levels,'near_sym_b','qshift_d');
    end
   
    % Perform motion estimation.
%     figure(2); % Define figures used within estshhemm22().
    [shift,delta,avec,qvecs] = estshhemm22(Zl,Zh,Xl,Xh,omega,nit,levelsel,0,[],mode,4,avlevel,sizeqfilt);
    
    % Perform registration of frame f by shifting Xl and Xh according to avec.
    Yl = Xl; Yh = Xh;
    for lev = length(Xh):-1:1
        srefh = size(Xh{lev});
        sc = srefh(1);
        sr = srefh(2);
        sh = affineshift2(avec,[sc sr],mode);  
        % Use sh (adjusted for pel units) to shift Xh.
        if lev <= 1, method = 'lin'; else method = 'cub'; end
        Yh{lev} = shift_cwt_bands2(Xh{lev},real(sh)*sr/2 + sqrt(-1)*imag(sh)*sc/2,method);
        if ~isempty(inputU)
            YhU{lev} = shift_cwt_bands2(XhU{lev},real(sh)*sr/2 + sqrt(-1)*imag(sh)*sc/2,method);
            YhV{lev} = shift_cwt_bands2(XhV{lev},real(sh)*sr/2 + sqrt(-1)*imag(sh)*sc/2,method);
        end
        % Shift Xl if at 2nd coarsest resolution.
        if lev == length(Xh) - 1
            Yl = shift_cwt_bands2(Xl,real(sh)*sr/2 + sqrt(-1)*imag(sh)*sc/2,'spl');
            if ~isempty(inputU)
                YlU = shift_cwt_bands2(XlU,real(sh)*sr/2 + sqrt(-1)*imag(sh)*sc/2,'spl');
                YlV = shift_cwt_bands2(XlV,real(sh)*sr/2 + sqrt(-1)*imag(sh)*sc/2,'spl');
            end
        end
    end
    
    % Take inverse DTCWT of registered subbands to obtain final registered frame.
    z = dtwaveifm2(Yl,Yh,'near_sym_b','qshift_d');
    if ~isempty(inputU)
        zU = dtwaveifm2(YlU,YhU,'near_sym_b','qshift_d');
        zV = dtwaveifm2(YlV,YhV,'near_sym_b','qshift_d');
    end
    
    % =========================================================================
    % resize by taking padding off
    if (padaddy~=0) || (padaddx~=0)
        z = z(padaddy + (1:m), padaddx + (1:n));
        if ~isempty(inputU)
            zU = zU(padaddy + (1:m), padaddx + (1:n));
            zV = zV(padaddy + (1:m), padaddx + (1:n));
        end
    end
    % =========================================================================
    
    % Save the registered frame and accumulate the average of these in zref.
    if length(mask)>1
        tempImg = imfilter(z, fspecial('average',blockSize/4));
        blurImg = z;
        blurImg(blockSize(1)/8+2:end-blockSize(1)/8-1,blockSize(2)/8+2:end-blockSize(2)/8-1) = ...
            tempImg(blockSize(1)/8+2:end-blockSize(1)/8-1,blockSize(2)/8+2:end-blockSize(2)/8-1);
        xrest(:,:,f) = z.*mask(:,:,f) + (1-mask(:,:,f)).*blurImg;
    else
        xrest(:,:,f) = z;
    end
    if ~isempty(inputU)
        xrestU(:,:,f) = zU;
        xrestV(:,:,f) = zV;
    end
    zref = zref + z;
    
    if update
        input(:,:,f) = z;
        if length(numToAvg) == 1
            xref = sum(input(:,:,1:numToAvg),3)/numToAvg; % Convert xref to an average frame.
        else
            if size(numToAvg,1) ~= size(input,1)
                xref = sum(input.*mask,3)./sum(mask,3);
            else
                xref = sum(input.*numToAvg,3)./sum(numToAvg,3);
            end
        end
        [Zl,Zh,Zsc] = dtwavexfm2(xref,levels,'near_sym_b','qshift_d');
        
    end
    
end
% zref = zref / length(ff); % Convert zref to an average.
