function [ val ] = Fh_conti(t,sigma)
% [ val ] = Fh( omega, T )
%
% Compute values of the Fourier transform of a Gaussian window
%
% 
% INPUT:
% t          : frequency values
% sigma      : window time spread parameter
%
% OUTPUT:
% val    : Fh(t) values


val = sigma^2*exp(-2*pi*t.^2*sigma.^2);


end