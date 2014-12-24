function [zrest, zrestsmooth] = Nick_pixel_fuse(xrest, levels, hpgain, phasenorm)


% Restore the registered images by performing a DTCWT
% and averaging the phases and taking max of magnitudes of the DTCWT coefs.
%
% Nick Kingsbury, Cambridge University, May-Aug 2011.

sx = size(xrest);

if nargin < 2
    levels = 4;
end

if nargin < 3
    hpgain = [2 1.4 1 1];
end
zgain = ones(6,1)*hpgain; % [2.5 1.8 1.2 1];%Set subband gains with some pre-emphasis for reconstructing output image, zrest.

if nargin < 4
    phasenorm = 0; % Use phase vector amplitude normalisation if phasenorm = 1.
end

Zh = cell(4,1); % Phase accumulators.
Zm = Zh;        % Magnitude accumulators.
Zm1 = Zh;        % Magnitude accumulators.
disp('Fusing frames together in DTCWT domain:');

% Do first frame separately from the loop.
[Xl,Xh] = dtwavexfm2(double(xrest(:,:,1)),levels,'near_sym_b','qshift_d');
Zl = Xl;
for k = 1:levels,
    xm = abs(Xh{k});
    if phasenorm,
        Zh{k} = Xh{k} ./ max(xm,1e-4);
    else
        Zh{k} = Xh{k};
    end
    Zm{k} = xm;
    Zm1{k} = xm;
%     Zang{k} = Zh{k}./ abs(Zh{k});
end

% Loop for frames 2 to the end.
for f = 2:sx(3),
    fprintf(' %d',f);
    [Xl,Xh] = dtwavexfm2(double(xrest(:,:,f)),levels,'near_sym_b','qshift_d');
    Zl = Zl + Xl; % Sum the lowpass bands
    for k = 1:levels,
        xm = abs(Xh{k});
        if phasenorm,
            Zh{k} = Zh{k} + Xh{k} ./ max(xm,1e-4); % Sum the unit phase vectors.
        else
            Zh{k} = Zh{k} + Xh{k}; % Sum the phase vectors.
        end
         Zm{k} = Zm{k} + xm; % Sum the magnitudes.
         changeMap = xm>Zm1{k};
         Zm1{k} = max(Zm1{k},xm); % Max the magnitudes.
    end
end
fprintf('\n');

% Fuse phases and magnitudes.
Zl = Zl / sx(3);


for k = 1:levels,
    ang{k} = Zh{k} ./ abs(Zh{k}); % Generate unit vectors with correct phases.
    Zh2{k} = ang{k} .* abs(Zm1{k}); % Scale by the max magnitudes.
end

% Perform inverse DTCWT on fused wavelet coefs.
zrest = dtwaveifm2(Zl,Zh2,'near_sym_b','qshift_d',zgain);
zrestsmooth = dtwaveifm2(Zl,Zh2,'near_sym_b','qshift_d');
