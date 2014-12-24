function out = gradIm(inIm,sigma)

GaussianDieOff = .001;  

% Design the filters - a gaussian and its derivative
   
pw = 1:30; % possible widths
ssq = sigma*sigma;
width = max(find(exp(-(pw.*pw)/(2*sigma*sigma))>GaussianDieOff));
if isempty(width)
   width = 1;  % the user entered a really small sigma
end
t = (-width:width);
len = 2*width+1;
t3 = [t-.5; t; t+.5];            % We will average values at t-.5, t, t+.5
gau = sum(exp(-(t3.*t3)/(2*ssq))).'/(6*pi*ssq);  % the gaussian 1-d filter
dgau = (-t.* exp(-(t.*t)/(2*ssq))/ ssq).';       % derivative of a gaussian 


% Convolve the filters with the image in each direction
% The canny edge detector first requires convolutions with
% the gaussian, and then with the derivitave of a gauusian.
% I convolve the filters first and then make a call to conv2
% to do the convolution down each column.
   
hei = size(inIm,1); 
wid = size(inIm,2); 

m = length(gau);
m2 = fix(m/2);

ye = wreflect([(1-m2):(hei+m2)], 0.5, hei+0.5); % Use 'reflect' so r < m2 works OK.
xe = wreflect([(1-m2):(wid+m2)], 0.5, wid+0.5); % Use 'reflect' so r < m2 works OK.

ax1 = conv2(inIm(ye,:),gau(:),'valid').'; 
ay1 = conv2(ax1(xe,:),dgau(:),'valid').'; 

ax2 = conv2(inIm(ye,:),dgau(:),'valid').'; 
ay2 = conv2(ax2(xe,:),gau(:),'valid').'; 

mag = sqrt((ay1.*ay1) + (ay2.*ay2));
% magmax = max(mag(:));
% if magmax>0
%   mag = mag / magmax;   % normalize
% end
out = mag;
   