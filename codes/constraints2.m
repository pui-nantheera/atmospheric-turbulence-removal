function cvecs = constraints2(refh,prevh,w6,mode,mask)

% function cvecs = constraints2(refh,prevh,w6,mode,mask)
%
% Generate the constraint vectors cvecs from the reference and previous
% subband coefs in refh and prevh.
% One set of 6 constraints is produced for each 2x2 group of subband coefs.
% w6(1:6,1:2) is the expected phase shift per sample for each subband
% which allows proper interpretation of measured phase shifts mod 2*pi.
% Note that the units of shift are assumed to be (half) the image size and
% not the pel spacing (since this changes from level to level).
% 
% Constraints2 uses adjacent 2x2 groups of coefs so that size(cvecs,1:2)
% is just one less than size(refh,1:2).
% size(cvecs,3) = 6 and size(cvecs,4) = 3 (the 1st, 2nd and last elements
% of the c vectors used in estshhemm2.m).
%
% mode is as defined in estshhemm22.m.
% mask is a matrix which multiplies the constraint vector amplitudes,
% typically so that some constraints can be masked out.
% If mask is smaller than (size(refh,[1 2])-1), then it is upsampled as required.
%
% July 2004: changed the normalised image size to be 2 units, as all the
% affine parts assume image offsets go from -1 to +1 across the image.
%
% Nick Kingsbury, Cambridge Univ, July 2004

shi = size(refh);

if any(rem(shi,2))
   error('subbands must be of even size for CONSTRAINTS.M');
end

sc = shi(1)-1;
sr = shi(2)-1;

tc = 1:sc;
tr = 1:sr;

% Increase size of mask as required to match size of 1st two dims of cvecs.
while size(mask,1) < shi(1)-1,
    % First increase no of rows, using [0.5 1 0.5] filter.
    sc = size(mask,1);
    msk = [];
    msk([0:sc]*2+1,:) = 0.5 * (mask([1 1:sc],:) + mask([1:sc sc],:));
    msk([1:sc]*2,:,:) = mask;
    mask =msk;
end
while size(mask,2) < shi(2)-1,
    % Now increase no of columns.
    sr = size(mask,2);
    msk = [];
    msk(:,[0:sr]*2+1) = 0.5 * (mask(:,[1 1:sr]) + mask(:,[1:sr sr]));
    msk(:,[1:sr]*2) = mask;
    mask = msk;
end

% Measure horizontal phase gradients by taking the angle of 
% summed conjugate products across horizontal pairs.
refcp = refh(:,tr+1,:) .* conj(refh(:,tr,:));
prevcp = prevh(:,tr+1,:) .* conj(prevh(:,tr,:));
hcp = refcp(tc,:,:) + refcp(tc+1,:,:) + prevcp(tc,:,:) + prevcp(tc+1,:,:);
% Measure horiz phase shifts to within +/- pi of the expected angles:
hdphi = zeros(size(hcp));
for d = 1:6
    w = w6(d,1); % expected horiz phase shift for band d.
    hdphi(:,:,d) = angle(hcp(:,:,d) * exp(-j*w)) + w;
end

% Measure vertical phase gradients by taking the angle of 
% summed conjugate products across vertical pairs.
refcp = refh(tc+1,:,:) .* conj(refh(tc,:,:));
prevcp = prevh(tc+1,:,:) .* conj(prevh(tc,:,:));
vcp = refcp(:,tr,:) + refcp(:,tr+1,:) + prevcp(:,tr,:) + prevcp(:,tr+1,:);
% Measure vert phase shifts to within +/- pi of the expected angles:
vdphi = zeros(size(vcp));
for d = 1:6
    w = w6(d,2); % expected vert phase shift for band d.
    vdphi(:,:,d) = angle(vcp(:,:,d) * exp(-j*w)) + w;
end

% Measure temporal phase differences between refh and prevh
% by taking the angle of summed conjugate products over 2x2 regions.
tcp = prevh .* conj(refh);
tcp = tcp(tc,tr,:) + tcp(tc+1,tr,:) + tcp(tc,tr+1,:) + tcp(tc+1,tr+1,:);
tdphi = angle(tcp); % Expected phase shifts are zero in this case.

% Calculate scale factors Cmag according to magnitudes of tcp, refh and prevh.
magsq = abs(refh).^3 + abs(prevh).^3;
Cmag = abs(tcp).^2 ./ (magsq(tc,tr,:)+magsq(tc+1,tr,:)+magsq(tc,tr+1,:)+magsq(tc+1,tr+1,:));

% Apply mask to Cmag in each direction.
for d=1:6, Cmag(:,:,d) = Cmag(:,:,d) .* mask; end

% Generate contraint vectors, with units of spatial shift normalised to half the image size
% (since shift goes from -1 to +1 across the image), and then scale all components by Cmag.
cvecs = cat(4,(vdphi*(shi(1)/2)).*Cmag,(hdphi*(shi(2)/2)).*Cmag,tdphi.*Cmag);

% Set bands 2 and 5 to zero.
if mode(3) == 1
	cvecs(:,:,[2 5],:) = 0;
end

% Remove the outer perimeter of constraints to reduce edge effects.
if mode(2) == 1
	cvecs([1 end],:,:,:) = 0; 
    cvecs(:,[1 end],:,:) = 0;
elseif mode(2) == 2  % Remove just bands nearly parallel to edges.
	cvecs([1 end],:,[2 3 4 5],:) = 0; 
    cvecs(:,[1 end],[1 2 5 6],:) = 0;
end

return;



   
      