function [unwrappedPhase]  = phase_unwrap_laplacian_2D(dataPhase, RefVox)
% [unwrappedPhase]  = phase_unwrap_laplacian_2D(dataPhase, RefVox)
%================================
% HEADER
%================================
%
% Name:             phase_unwrap_laplacian_2D.m
% Author:           Xu Li, PhD
% Email:            xuli@mri.jhu.edu
%
%================================
% PURPOSE
% unwrap phase according to the laplacian method
% i.e. del2(theta) = cos(theta)*del2(sin(theta)) + sin(theta)*del2(cos(theta))
%
%  Updated: 2012-08-17: Xu Li
%  Updated: 2014-01-03: JvB Added waitbars

warning off all

%% zero padding if necessary
N0 = size(dataPhase);

if nargin < 2
    RefVox = [floor(N0(1)/2), floor(N0(2)/2)];
end

RefFlag = 1;                % Do Ref in general
if RefVox == 0         
    RefFlag = 0;            % do not do Ref
end

if N0(1) < 256
    padsize = [32,32];
else     % high res
    padsize = [64, 64];
end

if RefFlag == 1
    c = - dataPhase(RefVox(1), RefVox(2));
else
    c = 0;
end


dimsOdd = mod(size(dataPhase), 2);
dataPhase = padarray(dataPhase, dimsOdd, 'replicate', 'post');  % padding to even dim first     
dataPhase = padarray(dataPhase, padsize, 'replicate');

N = size(dataPhase);                        % new dimentions after padding, now all even

ksize = [3, 3];               
khsize = (ksize-1)/2;

kernelcore = [0 1 0; 1 -4 1; 0 1 0];

Kernel = zeros(N);
Kernel( 1+N(1)/2 - khsize(1) : 1+N(1)/2 + khsize(1), ...
    1+N(2)/2 - khsize(2) : 1+N(2)/2 + khsize(2)) = kernelcore;

lap_opK = fftn(fftshift(Kernel));                                     % make real delta function, only do fftn

inv_lap_opK = zeros(size(lap_opK));
inv_lap_opK(lap_opK~=0) = 1./lap_opK(lap_opK~=0);                   % inverse of kernel in k-space

PhaseSin = sin(dataPhase);
PhaseCos = cos(dataPhase);    
PhaseLaplacian  = PhaseCos.*ifftn(  lap_opK.* fftn(sin(dataPhase))) - PhaseSin.*ifftn( lap_opK.* fftn(cos(dataPhase)));  %  

clear PhaseSin PhaseCos

unwrappedPhase  =  ifftn( inv_lap_opK.* fftn(PhaseLaplacian));   %  

unwrappedPhase = unwrappedPhase(padsize(1)+(1:N0(1)), padsize(2)+(1:N0(2)));               % crop back

%% shift to a known constant (find the largest phase region and select the middle point for reference)

unwrappedPhase = unwrappedPhase + c;        % phi' as in Schofield & Zhu, Optics Letters, 2003