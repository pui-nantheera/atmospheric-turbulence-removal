function y=drawrgb(a,posn);
% function y=drawrgb(a,posn);
% Function to display an image in a box matched to the pixel size.
% a  is assumed to be a 3-D true rgb matrix of uint8 values.
% If posn is specified, it determines the x and y position of the
% bottom left corner of the figure.
%
% Nick Kingsbury, Cambridge University, Sept 2001.

sz=get(0,'screensize');

xsize=sz(3);
ysize=sz(4);

sa=size(a);
m=sa(1); n=sa(2);

% If a is just a 2-D array,extend it to 3-D for a greyscale image. 
if length(sa) < 3, a = cat(3,a,a,a); end

% scale=round(min([xsize/(2*n),(ysize)/(15*m/8)]));
scale=round(min([512/n,512/m]));  % Better for 2^n image sizes.

figure(gcf);

if scale==0

clf reset
image(a);
set(gca,'position',[0.01 0.01 .98 .98]);
axis('off');
axis('image');
set(gcf,'position',get(0,'screensize'));

else

clf reset

image(a);
axis('off');
axis('image')

if nargin < 2, posn = []; end

if length(posn) == 2,
   pos=[posn(1)  posn(2)  scale*n+4  scale*m+4];
elseif length(posn) >= 4
   pos=[posn(1:2) posn(3:4)+4];
   scale = min(round(posn(3:4)./[n m]));
else
   figmax = fix((ysize - scale*m - 40)/20);
   ypos = (rem(gcf-1,figmax)+1)*20+14;   
   pos=[xsize-scale*n+12-gcf*20 ypos scale*n+4 scale*m+4];
end

set(gcf,'position',pos);

u=get(gca,'units');
pos=[2 2 scale*n scale*m];
set(gca,'units','pixels');
set(gca,'position',pos);
set(gca,'units',u);

end



