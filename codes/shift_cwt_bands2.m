function Yhsh = shift_cwt_bands2(Yh,v_field,interpmethod)

% function Yhsh = shift_cwt_bands2(Yh,v_field,interpmethod)
%
% Modify the 6 subbands of CWT coefs in Yh(:,:,6) so as to be equivalent
% to shifting the source image by the vector field, v_field. 
% Each subband is processed independently of the others. 
% v_field is specified in units of sample periods at the current 
% wavelet level.
% interpmethod is a string which specifies the rule used by interp2:
%   'nearest','linear','cubic','spline'  are valid methods in
% increasing order of computational complexity (defaults to 'cubic').
% 'linear' is about 4 times as fast as 'cubic'.
%
% If Yh is a 3-D array, then assume that it is the six bandpass bands
% and use complex bandpass interpolation, appropriate for each band.
% If Yh is only a 2-D array, then assume that it is the LoLo band
% in purely real format and use normal lowpass interpolation.
%
% Modified July 2004 to use phase wrapping and unwrapping, and a 
% more efficient interpolator (interp2).
%
% Nick Kingsbury, Cambridge University, July 2004.

if nargin < 3, interpmethod = 'cubic'; end

sy = size(Yh);
sv = size(v_field);
if length(sy) == 3,  % Bandpass band interpolation.
   if sy(3) ~= 6
      error('Yh must be 3-D array 6 subbands deep.');
   end
   if any(sv ~= sy(1:2))
      error('v_field must be same size as each subband in Yh.');
   end
   
   % Set up matrix extension sizes and padding arrays.
   extx = ceil(max(abs(real(v_field(:))))); 
   exty = ceil(max(abs(imag(v_field(:)))));
   exty(isnan(exty)) = 1;
   extx(isnan(extx)) = 1;
   exty(exty>(sv(1)*sv(2))) = 1;
   extx(extx>(sv(1)*sv(2))) = 1;
   
   z1 = zeros(exty,sv(2)+2*extx);
   z2 = zeros(sv(1),extx);

   % Create linear ramp matrices for phase wrapping.
   thx = ones(sv(1),1) * [1:sv(2)] + extx;
   thy = [1:sv(1)]' * ones(1,sv(2)) + exty;
   
   % Create matrices of interpolated point locations in the extended input matrix.
   xs = real(v_field) + thx;
   ys = imag(v_field) + thy;
   
   % Set up expected rotation rates for each subband with freq offset of -1/4 
   % and -3/4 times the sampling rate, to centre the passband of h on the
   % scaling func and wavelet passbands.
   jw0 = -sqrt(-1) * pi/2.15; % Nominally j * pi/2, but reduced a bit due to asymmetry of subband freq responses.
   jwx = jw0 * [1 3 3 3 3 1];
   jwy = jw0 * [3 3 1 -1 -3 -3];
   
   % Loop for each directional subband.
   Yhsh = zeros(sy);
   for d = 1:6
       % Get the subband, unwrap the phases and extend its borders.
       ye = [z1; z2 Yh(:,:,d).*exp(-thx*jwx(d)-thy*jwy(d)) z2; z1];
       % Interpolate ye to the new points, specified in (xs,ys).
       yi = interp2(ye,xs,ys,interpmethod);
       % Rewrap the phases.
       Yhsh(:,:,d) = yi .* exp(xs*jwx(d)+ys*jwy(d));
   end
   
else  % LoLo band interpolation - no need for any phase wrapping here.
   if size(v_field) ~= sy(1:2)
      error('v_field must be same size as each subband in Yh.');
   end
   
   % Set up matrix extension sizes.
   extx = ceil(max(abs(real(v_field(:))))); 
   exty = ceil(max(abs(imag(v_field(:)))));
   if exty>sy(1) exty=sy(1); end
   if extx>sy(2) extx=sy(2); end

   % Create matrices of interpolated point locations in the extended input matrix.
   xs = real(v_field) + ones(sv(1),1) * [1:sv(2)] + extx;
   ys = imag(v_field) + [1:sv(1)]' * ones(1,sv(2)) + exty;
   
   % Get the subband and extend its borders using symmetric extension.
   ye = Yh(max(1,[[exty:-1:1] [1:sv(1)] sv(1)+1-[1:exty]]),max(1,[[extx:-1:1] [1:sv(2)] sv(2)+1-[1:extx]]));
   % Interpolate ye to the new points, specified in (xs,ys).
   Yhsh = interp2(ye,xs,ys,interpmethod);
   
end

return

