clear all, clc, close all

Im = LoadIm('D:\Mfiles\images\OctecIR1.gif');

Im = Im(:,:,1);

Imout = imdistor(Im,[30 16 3],'grid','SNR');

SNR  = 20 * log10( norm(Im(:)) / norm(Im(:)-Imout(:)) )

figure
subplot(121), imshow(Im)
subplot(122), imshow(Imout)

% note: it is assumed that image pixels are [0 ... 255]
% if it is [0 1] then line 
% Imaxclip = 255;
% has to be changed into
% Imaxclip = 1;
% UNLESS there is no clipping in conversion to uint8, ie proper scaling is
% applied