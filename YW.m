function [h] = YW(y)
% [h] = YW(y)
%
% Solve a Yulle-Walker system to find the annihilating filter h.
%
% INPUT:
% y             : observation associated to a weighted sum of complex
%                 explonential. the weight are those of the Dirac, and the 
%                 modes to their locations.
%
% OUTPUT:
% h    : annihilating filter
%

Ly = floor(length(y)/2);

Xyw = toeplitz(y(Ly:end-1), y(Ly:-1:1));
Yyw = -1*(y(Ly+1:end)).';
A = Xyw \ transpose(Yyw);
h = [1 A.'];