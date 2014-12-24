function [shift,delta,av,qvecs] = estshhemm22(refl, refh, prevl, prevh, ...
   w, nit, levelsel, search, debug, mode, qscale, avlevel, sizeqfilt, mask)

% function [shift,delta,av,qvecs] = estshhemm22(refl, refh, prevl, prevh, ...
%    w, nit, levelsel, search, debug, mode, qscale, avlevel, sizeqfilt, mask)
%
% Estimate the shift between elements of the reference lowpass and
% highpass subimages, refl and refh, and the previous subimages, 
% prevl and prevh, at given levels of the dual-tree CWT.
% This version used Magnus Hemmendorf's affine vector field model
% from his 2002 Trans Med Imaging paper: 'Phase-based multi-dim volume
% registration'.
%
% refl and prevl are real lowpass subimages oversampled 2:1 each way.
% refh and prevh are cells of sets of 6 complex bandpass subimages in 3-D arrays.
% w(1:2) contains the expected phase rotations per sample 
% for the low and high band 1-D filters (in radians).
% nit is the number of iterations of the algorithm.
% levelsel selects the range of CWT levels used at each iteration:
% such that  levelsel(iter,1:2)  are the coarsest and finest levels
% used at iteration  iter.
% The default nit = 1.
% search  defines which offsets from the zero shift point will be tested
% in the initial search phase of the algorithm.  It is a column vector
% of complex elements representing the shift offsets.
%
% debug is an m x 2 matrix, specifying that element(s) debug(:,1),debug(:,2)
% in the subimages should be printed out for debugging. 
%
% shift is a complex matrix with 1+1j representing the lower right
% corner of a unit square (since vertical distances are measured
% downwards in images).  The shift vector points from a sampling
% point in ref to the equivalent point in prev.
% delta is a matrix giving the height of the error surface minimum
% at each point.
% av(:,:,1:6) is a 3-D array of affine vectors at each point.
% qvecs(:,:,1:28) is a 3-D array of q-vectors at each point. 
%
% mode(1)>0 causes a single affine vector to be used for the whole image
%  for iterations up to mode(1).
% mode(2)=1 removes the constraints for the outer edge of avecs
% mode(2)=2 removes the constraints for the outer edge of avecs on only the
%  subbands which are not normal to the edge of the image.
% mode(3)=1 removes constraints from bands 2 and 5 everywhere.
% mode(4)=0 for pure affine model, 1 for (affine + xy) model, 2 for
%  (affine + linear zoom) model, and 3 for 10-term affine + linear zoom model.
% mode(5)=0.5 for modified zoom (1.0 otherwise).
% qscale defines the scale factor for the quiver plots, and equals 1 by default.
% avlevel is the level of the DTCWT subbands which have the same resolution
%  as the avec matrix.
% sizeqfilt is the size of the Q matrix smoothing filter 
%  (no. of taps = 2*sizeqfilt + 1).
% mask is a gain matrix for scaling the contraint vectors.  It should be
% the size of the coarsest subbands - 1.  It allows parts of the image to be ignored
% where mask = 0, and is usually = 1 elsewhere.
%
% esthemmm2 uses constraints2 for a more dense constraint vector field.
% The July 2004 version allows the avec field to be at finer resolution
%  than the coarsest DT CWT subband resolution.
%
% Nick Kingsbury, Cambridge University, July 2004.

if nargin < 14, mask = 1; end
if nargin < 13, sizeqfilt = 2; end
if nargin < 12, avlevel = length(refh); end
if nargin < 11, qscale = 1; end
if nargin < 10, mode = 0; end
if length(mode) < 5, mode = [mode 0 0 0 0]; end 

if nargin < 9, debug = []; end

if nargin < 8, search = 0; end

% fig = gcf + 1;

% Check matrix sizes.
slo = size(refl);
shi = size(refh{end});
% if any(size(prevl)~=slo | slo~=2*shi(1:2)) | any(size(prevh{end})~=shi) | (shi(3)~=6),
%   error('Input matrices not correct sizes.');
% end

% Mask matrix should be size of coarsest high bands - 1.
if length(mask)==1, mask = mask*ones(shi(1:2)-1); end
 
% Subimages are listed in order of increasing angle of rotation
% from 0 to pi, starting at 15 deg anticlockwise from a horizontal edge,
% and increasing by 30 deg each time.

nd = size(debug,1);
if nd > 0, di = debug(:,1) + (debug(:,2) - 1) * sc;
else       di = [];
end
mat2scr([debug di],'%8d','ESTSHIFT debug at pel:')
% mat2scr(r6(di,:),'%8.2f','r6:')

% p6 = prevh;

% Calc the expected horiz and vert phase shifts for each subband.
w6 = [w(1)  w(2); w(2)  w(2); w(2)  w(1); w(2) -w(1); w(2) -w(2); w(1) -w(2)];

nsch = length(search);

if mode(4) == 0, lenavec = 6;
elseif mode(4) == 3, lenavec = 10;
else lenavec = 8;
end

% Loop for each search offset.
for sch = 1:nsch,
   
   srefh = size(refh{avlevel});
   sc = srefh(1);
   sr = srefh(2);
   
   % Initialise affine vector to search offset.
   av = zeros(sc-1,sr-1,lenavec);
   av(:,:,2) = real(search(sch));
   av(:,:,1) = imag(search(sch));
   delta1 = [];
   
   % Perform nit iterations on shift and del so that del is more accurate.
   for it = 1:nit,
      levels = levelsel(min(it,size(levelsel,1)),:);
      qvecsum = 0;
      % Loop for each level at the current iteration.
      for lev = levels(1):-1:levels(2)
         
         % Calc. size of bandpass subimages in pels.
         srefh = size(refh{lev});
         sc = srefh(1);
         sr = srefh(2);
         ss = sc * sr;
         
         sh = affineshift2(av,[sc sr],mode);  
         
         % pui ----------------------------------
         if (size(sh,1)~=sc)||(size(sh,2)~=sr)
             temp = sh;
             sh = zeros(sc,sr);
             sh(1:min(sc,size(temp,1)),1:min(sr,size(temp,2))) = temp(1:min(sc,size(temp,1)),1:min(sr,size(temp,2))); clear temp
         end
         % end ----------------------------------
         
         % Use sh (adjusted for pel units instead of half-image units) to shift prevh.
         if lev <= 1, method = 'lin'; else method = 'cub'; end
         prevhsh = shift_cwt_bands2(prevh{lev},real(sh)*(sr/2) + imag(sh)*(sqrt(-1)*sc/2),method);
         
         % Calculate motion constraints at each point.
         % cvecs is a 4-D array in (y,x,band,cvec element).
         cvecs = zeros(sc-1,sr-1,6,lenavec+1);
         cvecs(:,:,:,[1 2 lenavec+1]) = constraints2(refh{lev},prevhsh,w6,mode,mask);
         
         % Add extra components which are equivalent to multiplying 
         % the constraint vectors cvecs by K.' matrices:
         % K = [1 0 y 0 x 0 xy 0 0; 0 1 0 y 0 x 0 xy 0; 0 0 0 0 0 0 0 0 1] 
         % for 2-D affine transform with extra xy terms.
         % yi and xi go from -1 to +1 across the image.
         % For the linear zoom model, columns 7 and 8 of K (the xy cols) are 
         % replaced by terms in [xy  x^2  0]' and [y^2  xy  0]' .
         s2 = sc - 1;
         for k=1:s2,
            yi = (2*k-sc) / sc;
            cvecs(k,:,:,[3 4]) = cvecs(k,:,:,[1 2]) * yi; 
            if mode(4) == 3, 
               cvecs(k,:,:,10) = cvecs(k,:,:,3) * yi; 
            elseif mode(4) == 2, 
               cvecs(k,:,:,8) = cvecs(k,:,:,3) * (mode(5)*yi); 
            end
         end
         s2 = sr - 1;
         for k=1:s2,
            xi = (2*k-sr) / sr;
            cvecs(:,k,:,[5 6]) = cvecs(:,k,:,[1 2]) * xi; 
            if mode(4) == 3, % Linear zoom with separate x^2 and y^2 terms.
               cvecs(:,k,:,7) = cvecs(:,k,:,3) * xi;
               cvecs(:,k,:,8) = cvecs(:,k,:,4) * xi;
               cvecs(:,k,:,9) = cvecs(:,k,:,6) * xi;
            elseif mode(4) == 2, % Linear zoom
               cvecs(:,k,:,7) = cvecs(:,k,:,3) * xi + cvecs(:,k,:,6) * (mode(5)*xi);
               cvecs(:,k,:,8) = cvecs(:,k,:,4) * xi + cvecs(:,k,:,8);
            elseif mode(4) == 1,
               cvecs(:,k,:,[7 8]) = cvecs(:,k,:,[3 4]) * xi;
            end
            
         end
         % cvecs now contains the terms  K' * c .
         
         % Form the Q-matrices as vectors of the lower left half
         % of the symmetric Q-matrix.  Q = (K'*c*c'*K)
         % Combine all the constraints (subbands) at each location, ie all 6 directions.
         scv = size(cvecs);
         c4 = scv(4);
         qvecs = zeros([scv(1:2) c4*(c4+1)/2]);
         iQ = zeros(c4,c4);
         k = 1;
         for k1 = 1:c4
            for k2 = k1:c4
               qvecs(:,:,k) = sum(cvecs(:,:,:,k1) .* cvecs(:,:,:,k2),3); % Sum over all 6 directions.
               iQ(k1,k2) = k; iQ(k2,k1) = k; % indices for Q matrix and q vector.
               k = k + 1;
            end
         end
         
         % Change the sample rate of qvecs, back to that of av, remembering
         % that the sizes are always odd and one less than the subband sizes.
         if lev < avlevel,
             % Reduce the sample rate.
             for k = lev:(avlevel-1)
                 % First reduce no of rows, using [0.5 1 0.5] filter.
                 t1 = 2:2:(size(qvecs,1)-1);  
                 qvecs = 0.5*qvecs(t1-1,:,:) + qvecs(t1,:,:) + 0.5*qvecs(t1+1,:,:);
                 % Now reduce no of columns.
                 t2 = 2:2:(size(qvecs,2)-1);  
                 qvecs = 0.5*qvecs(:,t2-1,:,:) + qvecs(:,t2,:) + 0.5*qvecs(:,t2+1,:);
             end
         elseif lev > avlevel,
             for k = (avlevel+1):lev
                 % First increase no of rows, using [0.5 1 0.5] filter.
                 sc = size(qvecs,1);
                 qv = [];
                 qv([0:sc]*2+1,:,:) = 0.5 * (qvecs([1 1:sc],:,:) + qvecs([1:sc sc],:,:));
                 qv([1:sc]*2,:,:) = qvecs;
                 qvecs = qv;
                 % Now increase no of columns.
                 sr = size(qvecs,2);
                 qv = [];
                 qv(:,[0:sr]*2+1,:) = 0.5 * (qvecs(:,[1 1:sr],:) + qvecs(:,[1:sr sr],:));
                 qv(:,[1:sr]*2,:) = qvecs;
                 qvecs = qv;
             end
         end
          
         qvecsum = qvecsum + qvecs; % Accumulate qvecs across levels of CWT.
         
      end  % End of loop for each level of CWT.
      
      % As a simple test, combine all qvecs into a single vector.
      qtot = squeeze(sum(sum(qvecsum)));
      
      % Filter the q vectors with a linear ramp filter.
      % Set up the 1-D filter impulse response.
      shq = sizeqfilt(min(it,end)) + 1; % 
      hq = [1:shq  (shq-1):-1:1].'; 
      hq = hq/sum(hq);
      for k=1:size(qvecsum,3),
         qvecsum(:,:,k) = colfilter(colfilter(qvecsum(:,:,k),hq).',hq).';
         % qvecsum(:,:,k) = conv2(hqc,hqr,squeeze(qvecsum(:,:,k)),'same');
         if mode(1) >= it,
            qvecsum(:,:,k) = qtot(k); % Include this to use single Q matrix.
         end
      end
      
      % Now extract Q matrix, q vector and q0 scalar from each qvec vector, and
      % solve for the affine vector, a, which minimises the constraint error energy.
      iQmat = iQ(1:(end-1),1:(end-1));  % Indices to convert qvec to Q.
      iqv = iQ(1:(end-1),end);  % Indices to convert qvec to q.
      iq0 = iQ(end,end);  % Index to convert qvec to q0.
      sqv = size(qvecsum);
      del = zeros(sqv(1:2));
      avec = zeros([sqv(1:2) lenavec]);
      qf = zeros(sqv(3),1);
      % Loops for each Q matrix locality.
      for row = 1:sqv(2)
         for col = 1:sqv(1)
            % Extract Q, q and q0.
            qf(:) = qvecsum(col,row,:);
            Qmat = qf(iQmat);
      		qv = qf(iqv);
      		q0 = qf(iq0);
            % Solve for a_min: 
      		a = -Qmat \ qv; 
            avec(col,row,:) = a;
            % Calculate the error energy:
            del(col,row) = a.'*Qmat*a + 2*qv.'*a + q0;
         end
      end
      
      % Update av.
      % pui -------------------------------
      if sum(size(av)~=size(avec))
          newsize = max(size(av),size(avec));
          if sum(size(av)<newsize)
              temp = av;
              av = zeros(newsize);
              av(1:size(temp,1),1:size(temp,2),1:size(temp,3)) = temp; clear temp
          end
          if sum(size(avec)<newsize)
              temp = avec;
              avec = zeros(newsize);
              avec(1:size(temp,1),1:size(temp,2),1:size(temp,3)) = temp; clear temp
          end
      end
      % end -------------------------------
      av = av + avec;
      
      % Create shift vector at same resolution as av.
      sav = size(av) + [1 1 0];
      sh1 = affineshift2(avec,sav,mode);
	  
      
      % Update shift.
      sh = affineshift2(av,sav,mode);

%       % Plot shift vectors.
%       figure(fig); set(gcf,'DefaultLineLineWidth',1.5);
%       quiver(-real(sh)*sav(2)/2*qscale,-imag(sh)*sav(1)/2*qscale,0); 
%       set(gcf,'position',[700  30  500  672]);
%       set(gca,'position',[0.01 0.01 .98 .98]);
%       axis equal;
%       axis ij;
%       axis([-1 sav(2)+2 -1 sav(1)+2]);
%       settitle(sprintf('%d iteration, output shift vectors',it));
%       drawnow
      % fig = fig + 1;
      
      % Monitor the mean squared error and the mean affine vector (in pels).
      delta1(it)=mean(del(:));
      avec1=squeeze(mean(mean(av,1),2));
      avec1(1:2:end) = avec1(1:2:end)*size(refh{1},1);
      avec1(2:2:end) = avec1(2:2:end)*size(refh{1},2);
      % levels,delta1(it),avec1
      
      % Debug info:
      mat2scr([real(sh(di))*16 imag(sh(di))*16 del(di)],'%10.4f',...
         're(sh)*16  im(sh)*16  del:')
      
   end  % End of iterations.
   
   % for k=1:6,
   %   t=[28:29 36:37];
   %   k,[r6(t,k)  p6i(t,k)  p61(t,k)]
   % end
   
   % Update shift, curv and delta at points where del < delta and sh is
   % close enough to shi.
   if sch == 1,
      shift = sh;
      delta = del;
   else
      if rem(mode,2) == 1,
         dmin = find(del./(sum(esurf(:,1:2)')'+1e-6) < ...
            (delta(:)./(sum(curv(:,1:2)')'+1e-6)-1e-5) & absl1(sh(:) - shi(:)) < 0.6);
      else
         dmin = find(del < (delta(:) - 1e-2) & absl1(sh(:) - shi(:)) < 0.6);
      end
      %    if ~isempty(dmin),
      %      s=search(sch)
      %      [dmin(:) shift(dmin)*16 sh(dmin)*16]
      %      [dmin(:) delta(dmin) del(dmin)]
      %    end
      shift(dmin) = sh(dmin);
      delta(dmin) = del(dmin);
   end
   
   mat2scr([real(shift(di))*16 imag(shift(di))*16 delta(di)],'%10.4f',...
      'After min, re(shift)*16  im(shift)*16  delta:')
   
end  % End of search offsets

% keyboard
% r6a=r6(28,:).';p6a=p6(4,[4:8:48]).';[r6a p6a r6a.*conj(p6a)/10]

% Uncomment next line to plot mean value of delta vs iteration.
% figure; semilogy(delta1); grid on; pause(0.1)

return




