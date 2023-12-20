clear all
close all
clc
folder = './';
addpath(folder);
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'FRI_lib']));


N  = 1024;         %% signal length
t  = (0:N-1)/N;

%% Time-frequency representation parameters
%M       = 500;       %% Number of frequential bin
%L       = 20;        %% analysis window size (in bin)

% First component - sinusoidal
A = 1;
freq1 = 340.5;
freq2 = 240.5;
freq3 = 220.5;
X(:,1) = A*exp(2*1i*pi*freq1*t);
% second component - sinusoidal
B = 2;
X(:,2) = B*exp(2*1i*pi*freq2*t);
C = 3;
X(:,3) = C*exp(2*1i*pi*freq3*t);
x0 = sum(X,2);

P = 3;

gamma=0; 
sigma = 0.04;
%Nfft = N/4;
Nfft = N; 
M = Nfft;
%on calcule la FSST 
[STFT,~,~,~,~] = sst2(x0,sigma,Nfft,gamma);
%[tfr] = tfrgab2(x0, M, L); 
spect = (abs(STFT)/N).^2; %% compute the discrete spectrogram, be careful of to the normalization

%% approximation of the spectrogram
%computattion of the window

%we translate val1 by a factor of first_index
val1 = Fh_conti((0:M-1)-freq1,sigma);
%we translate val1 by a factor of second_index
val2 = Fh_conti((0:M-1)-freq2,sigma);
val3 = Fh_conti((0:M-1)-freq3,sigma);

val_approx = A.^2*val1+B.^2*val2+C.^2*val3; %neglecting interference

spect_appr = transpose(val_approx)*ones(1,N);

%we approximation (9) of the paper
M_0 = 20; %order of truncation
%To compute the Fourier coefficients, we compute the Fourier coeffcients
Cp =  sigma/(M*sqrt(2))* exp(-pi*((-M_0:M_0).^2/(2*(sigma*M).^2)))';
%then we recompute approaximation of the spectrogram obtained
Mat1 = exp(2*pi*1i*((0:M-1)-freq1)/M);
Mat2 = exp(2*pi*1i*((0:M-1)-freq2)/M);
Mat3 = exp(2*pi*1i*((0:M-1)-freq3)/M);
A_acc = [];
B_acc = [];
C_acc = [];
for p = -M_0:M_0 
  A_acc = [A_acc ; Mat1.^p];  
  B_acc = [B_acc ; Mat2.^p];
  C_acc = [C_acc ; Mat3.^p];
end

approx2 = A.^2*transpose(A_acc)*Cp + B.^2*transpose(B_acc)*Cp + B.^2*transpose(C_acc)*Cp;
spect_approx2 = abs(approx2)*ones(1,N);

%% computation of the localisation of the localisation of the IFs 

V = exp(2*pi*1i*((0:M-1)'/M*(-M_0:M_0)));
D = diag(Cp); 
W = V*D;
L = zeros(2*M_0+1,N);
for n = 1:N
 L(:,n) = W\spect_appr(:,n);
end 

%annihilating filter technique
FoundDiracsLocations = zeros(P,N);
Amp = zeros(P,N);

for n = 1:N
 y = L((M_0+1)-(P-1):(M_0+1)+P,n); %we pick 2P values from L, considering centered frequencies 
 [h] = YW(y);

 TZroot = roots(h(end:-1:1));
 X =(angle(TZroot)/(2*pi)) > 0;
 Y = angle(TZroot)/(2*pi).*X+(1+angle(TZroot)/(2*pi)).*(1-X);
 FoundDiracsLocations(:,n) = sort(Y*M,'descend');
 
 %amplitude computation
 if (n == 200)
  X = FoundDiracsLocations(:,n)
  pause
 else
  X = FoundDiracsLocations(:,n);
 end
 W = exp(-2*1i*pi*(0:P-1)'*X'/M);
 Amp(:,n) = real(W\L((M_0+1):(M_0+1)+P-1,n)); 
end

% figure
phi1prim = freq1*ones(1,length(t));
phi2prim = freq2*ones(1,length(t)); 
phi3prim = freq3*ones(1,length(t));
plot(t,phi1prim,t,phi2prim,'--',t,phi2prim,'-.',t,FoundDiracsLocations(1,:),'-.',t,FoundDiracsLocations(2,:),'-o',...
t,FoundDiracsLocations(3,:),'-*');
figure;
imagesc((0:N-1)/N,(0:M-1)/M*N,spect_appr);
set(gca,'ydir','normal');
% comparison spectrogram and spect_approx
figure
imagesc((0:N-1)/N,(0:M-1)/M*N,spect)
set(gca,'ydir','normal');
figure 
plot(t,Amp(1,:),'-.',t,Amp(2,:),'-o',t,Amp(3,:),'-*');