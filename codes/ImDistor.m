function [I, opt] = ImDistor(I,opt,type,snrstr)

%ImDistor.m, 2004/06/03 Artur Loza
%
%===========================================================================
% 
%     function Iout = ImDistor(I,opt,type,snrstr)
%
% Distortion of the image
%
%  Inputs:  Iin   - Input Image; uint8 or double
%           opt   - options depending on the type of noise
%           type  - 'grid',  opt = [ampl freq width]
%                   'gwn' ,  opt = [mean stdev]
%                   'saltp', opt =  density
%                   'jpg',   opt =  quality [0 100]
%           snrstr- if specified ('SNR' or 'PSNR'), then 
%           opt(1) is the (P)SNR value
%
% Outputs:  Iout  - Output Image
%
%===========================================================================
%
% Example: Iin=imread('pout.tif');
%          Iout = imdistor(Iin,[30 8 4],'grid','PSNR');
%          figure, 
%          subplot(121), imshow(Iin), 
%          subplot(122), imshow(Iout)
%
%===========================================================================    

% opt vectors:
%           'grid':
%           ampl  - constant amplitude of the grid; any value, also
%                   negative
%           freq  - number of lines in the grid, negative values are
%                   ignored
%           width - width of the grid line, negative values are 
%                   ignored
%
%           'gwn':
%           sdev  - standard deviation
%           mn    - mean
%
%           'saltp':
%           dens  - desity

if nargin == 4 % (P)SNR has been specified by the user
  
  snrv = num2str(opt(1));
  
  switch lower(type(1:3))
    case 'gri', [I, opt] = ImGridistor(I,snrv,opt(2),opt(3),snrstr);
    case 'gwn', [I, opt] = ImGWNdistor(I,snrv,opt(2),snrstr);
    case 'sal', [I, opt] = ImSaltPepperdistor(I,snrv,snrstr);
    case 'jpg', [I, opt] = ImJpegdistor(I,opt(1));
  end
  
else
  
  switch lower(type(1:3))
    case 'gri', [I, opt] = ImGridistor(I,opt(1),opt(2),opt(3));
    case 'gwn', [I, opt] = ImGWNdistor(I,opt(1),opt(2));
    case 'sal', [I, opt] = ImSaltPepperdistor(I,opt(1)); 
    case 'jpg', [I, opt] = ImJpegdistor(I,opt(1));
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% helper files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [I, ampl] = ImGridistor(I,ampl,freq,width,snrstr)

%imgridistor.m, 2004/06/03 Artur Loza
%
%===========================================================================
% 
%     function Iout = imgridistor(Iin,ampl,freq,width)
%
% Periodic (horizontal) grid distortion of the image
%
%  Inputs:  Iin   - Input Image; uint8 or double
%           ampl  - constant amplitude of the grid; any value, also
%                   negative [40], if ampl is a string, eg. '24'
%                   amplitude will be calculated to give SNR = 24 dB
%           freq  - number of lines in the grid, negative values are
%                   ignored [N/32]
%           width - width of the grid line, negative values are 
%                   ignored [3]
%
% Outputs:  Iout  - Output Image
%
%===========================================================================
%
% Example: Iin=imread('pout.tif');
%          Iout1 = imgridistor(Iin,50,6,4);
%          Iout2 = imgridistor(Iin,-30);
%          figure, 
%          subplot(131), imshow(Iin), 
%          subplot(132), imshow(Iout1)
%          subplot(133), imshow(Iout2)
%
%===========================================================================    

rc = size(I);

if nargin < 4, width = 3; end
if nargin < 3, freq = round(rc(1)/32); end
if nargin < 2, ampl = 40; end

if (freq < 1) | (width < 1)
  
  ind = [];
  
else
  
  % row indices of the grid
  ind = round(linspace(1,rc(1)-width,freq))';
  ind = repmat(ind,1,width) + repmat([0:width-1],length(ind),1);
  ind = sort(ind(:));
  
end

if isstr(ampl), % SNR was specified as 'ampl'
  
  SNRg = str2num(ampl); 
  snrstr = lower(snrstr);
  Id = double(I);
  Iind = double(I(ind,:));
  nind = rc(2)*length(ind);
  
  if snrstr(1) == 's'
    ampl = norm(Id(:)) / ( 10^(SNRg/20)*sqrt(nind) );
  elseif snrstr(1) == 'p'
    ampl = max(Id(:)) / ( 10^(SNRg/20)*sqrt(nind/prod(rc)) ) ;        
  end
  
  ampl0 = [];
  Imaxclip = 255;
  
  % assuming that images are used  
  %  if strcmp(class(I),'uint8') 
  % correction of clipping and rounding effects
  
  clipim = find((Iind+round(ampl)) >= Imaxclip);
  
  while ~isempty(clipim) % estimation error occurs
    
    rind = nind-length(clipim); 
    if rind < 1, 
      warning('SNR requires wider and/or more dense grid'), 
      break
    end
    
    clipimv = Imaxclip-Iind(clipim);
    
    if snrstr(1) == 's'
      ampl = sqrt( (         norm(Id(:))^2/10^(SNRg/10) - norm(clipimv)^2) / rind );
    elseif snrstr(1) == 'p'
      ampl = sqrt( ( prod(rc)*max(Id(:))^2/10^(SNRg/10) - sum(clipimv.^2)) / rind);
    end
    ampl = round(ampl);
    
    if length(clipim) == length(find((Iind+ampl) >= Imaxclip))
      clipim = [];
    else
      clipim = find((Iind+ampl) >= Imaxclip);           
    end
    
    if any(ampl0 == ampl), 
      break % stop if infinite loop (due to rounding)
    else
      ampl0 = [ampl0 ampl];
    end     
    
  end % while
  
  % assuming that images are used     
  %   end
  % ampl = round(ampl);
  
end

if strcmp(class(I),'uint8')   
  
  I(ind,:) = uint8( double(I(ind,:)) + ampl );
  
else % double ?
  
  I(ind,:) = I(ind,:) + ampl;
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [I, sdev] = ImGWNdistor(I,sdev,mn,snrstr)

%ImGWNdistor.m, 2004/06/03 Artur Loza
%
%===========================================================================
% 
%     function Iout = ImGWNdistor(Iin,sdev,mn)
%
% Additive GWN distortion of the image
%
%  Inputs:  Iin   - Input Image; uint8 or double
%           sdev  - standard deviation [0.1], if sdev is a string, eg. '24'
%                   std will be calculated to give SNR = 24 dB
%           mn    - mean [0]
%
% Outputs:  Iout  - Output Image
%
%===========================================================================
%
% Example: Iin=imread('pout.tif');
%          Iout = ImGWNdistor(Iin);
%          figure, 
%          subplot(121), imshow(Iin), 
%          subplot(122), imshow(Iout)
%
%===========================================================================    

rc = size(I);

if nargin < 3, mn = 0; end
if nargin < 2, sdev = .1; end

if isstr(sdev), % SNR was specified as 'ampl'
  
  SNRg = str2num(sdev);
  Id = double(I(:));
  nind = prod(rc);
  snrstr = lower(snrstr);
  
  if snrstr(1) == 's'
    sdev = norm(Id) / ( sqrt(nind)*10^(SNRg/20) ); 
    if abs(mn) > 0 % the same
      sdev = sqrt(norm(Id)^2 / (nind*10^(SNRg/10)) - mn); 
    end
  elseif snrstr(1) == 'p'
    sdev = max(Id(:)) / 10^(SNRg/20);    
    if abs(mn) > 0 % the same
      sdev = sqrt(max(Id(:))^2 / 10^(SNRg/20) - mn^2 );    
    end
  end
  
end

noi  = sdev*randn(rc)+mn;

  % correction of clipping and rounding effects - tbi

if strcmp(class(I),'uint8')
  
  I = uint8(double(I) + noi);
  
else % double ?
  
  I = I + noi;
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [I, dens] = ImSaltPepperdistor(I,dens,snrstr)

%ImGWNdistor.m, 2004/06/03 Artur Loza
%
%===========================================================================
% 
%     function Iout = ImSaltPepperdistor(Iin,dens)
%
% Salt and Pepperdistor.m distortion of the image
%
%  Inputs:  Iin   - Input Image; uint8 or double
%           dens  - desity [0.1], if dens is a string, eg. '24'
%                   dens will be calculated to give SNR = 24 dB
%
% Outputs:  Iout  - Output Image
%
%===========================================================================
%
% Example: Iin=imread('pout.tif');
%          Iout = ImSaltPepperdistor(Iin);
%          figure, 
%          subplot(121), imshow(Iin), 
%          subplot(122), imshow(Iout)
%
%===========================================================================    

rc = size(I);

if nargin < 3, mn = 0; end
if nargin < 2, dens = .1; end

if isstr(dens), % SNR was specified 
  
  SNRg = str2num(dens);
  Id = double(I(:));
  nind = prod(rc);
  snrstr = lower(snrstr);
  A = max(Id(:))-mean(Id(:));
  
  if snrstr(1) == 's' 
    dens = norm(Id)^2 / ( A^2*nind*10^(SNRg/10) ); 
  elseif snrstr(1) == 'p'
    dens = max(Id(:))^2 / ( A^2*10^(SNRg/10) );     
  end
  
end

I = imnoise(I,'salt & pepper',dens);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [I, qlty] = ImJpegdistor(I,qlty)

%ImJpegdistor.m, 2004/06/03 Artur Loza
%
%===========================================================================
% 
%     function Iout = ImJpegdistor(Iin,qlty)
% JPEG distortion of the image
%
%  Inputs:  Iin   - Input Image; uint8 or double
%           qlty  - quality value [0 100]
%
% Outputs:  Iout  - Output Image
%
%===========================================================================
%
% Example: Iin=imread('pout.tif');
%          Iout = ImJpegdistor(Iin,50);
%          figure, 
%          subplot(121), imshow(Iin), 
%          subplot(122), imshow(Iout)
%
%===========================================================================    

imwrite(uint8(I),'temp','jpeg','Quality',qlty);
I = imread('temp','jpeg');
dos('delete temp');