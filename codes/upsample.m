function y = upsample(x,factor,h1,h2)

% function y = upsample(x,factor,h1,h2)
%
% Up-sample the matrix x by factor(1:2) in each direction, placing each
% factor(1)*factor(2) group of pels in y directly below each pel in x.
% h1 and h2, if given, are the filtering functions for columns and rows.
% factor defaults to 2x2 up sampling, if not given.
% h1 and h2 default to bilinear interpolation. 
%
% Nick Kingsbury, Cambridge University, Dec 2001.

% Default to 2:1 upsampling.
if nargin < 2, factor = [2 2]; end
% Default to equal upsampling in both directions if factor is a scalar.
if length(factor) < 2, factor = factor*[1 1]; end

if all(factor == [1 1])
   y = x;
   return;
end

if nargin < 3,  % Define h1 and h2 as rectangular pulses for bilinear interpolation.
   h1 = ones(factor(1),1) / factor(1);
   h2 = ones(factor(2),1) / factor(2);
   % Adjust to make lengths of h1 and h2 odd to get correct alignment.
   if rem(factor(1),2) == 0, h1 = conv(h1,[1;1]/2); end
   if rem(factor(2),2) == 0, h2 = conv(h2,[1;1]/2); end
end

y = kron(x,ones(factor));
y = colfilter(y,h1).';
y = colfilter(y,h2).';

return


