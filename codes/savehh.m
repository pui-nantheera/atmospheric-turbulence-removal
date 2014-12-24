% savehh.m
%
% Save the frames of hillhouse as separate tif images in order to test
% the software for dstl.
% Adjust the frame size to 480 x 512 so we can check that all works ok
% with non-square frames.
%
% Nick Kingsbury, Cambridge University, May 2011.

load xref

xsize = size(xin);
m = xsize(1);
n = xsize(2);
N = xsize(3);
dirname = 'images';
filename = 'frame1000';

t = [33:m]';  % Reduce height of frame.
fn = length(filename) + [-3:0];
for f = 1:N,
    filename(fn) = sprintf('%4d',1000+f); %Adjust last 4 digits of filename.
    x =  uint8(xin(t,:,f));
    imwrite(x,[dirname '\' filename '.tif'],'tiff','compression','none');
end
