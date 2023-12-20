function test_FRI_classique (cas,P)
 %INPUT 
 % cas : 1 we deal with pure harmonics, 2 linear chirp, 3 hyperbolic chirp + linear chirp
 % P   :  number of mode to be detected
%  folder = './';
%  addpath(folder);
%  addpath(strcat([folder 'synchrosqueezedSTFT']));
%  addpath(strcat([folder 'tools']));
%  addpath(strcat([folder 'FRI_lib']));
% 

N  = 1024;         %% signal length
t  = (0:N-1)/N;

if cas == 1
 
 N  = 1024;         %% signal length
 t  = (0:N-1)/N;
 freq1 = 340.5;
 freq2 = 240.5;
 freq3 = 220.5;
 A = 3;
 X(:,1) = A*exp(2*1i*pi*freq1.*t);
 B = 2;
 X(:,2) = B*exp(2*1i*pi*freq2.*t);
 C = 1;
 X(:,3) = C*exp(2*1i*pi*freq3.*t);
 x0 = sum(X,2);
 phi1prim = 340.5*ones(1,N);
 phi2prim = 240.5*ones(1,N);
 phi3prim = 220.5*ones(1,N);
 phiprim = zeros(3,N); 
 phiprim(1,:) = phi1prim;
 phiprim(2,:) = phi2prim;
 phiprim(3,:) = phi3prim;
 gamma=0; 
 sigma = 0.04;
elseif cas == 2
 N  = 1024;         %% signal length
 t  = (0:N-1)/N;
 freq1 = 260.5+100*t;
 freq2 = 240.5+100*t;
 A = 2;
 X(:,1) = A*exp(2*1i*pi*freq1.*t);
 B = 2;
 X(:,2) = B*exp(2*1i*pi*freq2.*t);
 x0 = sum(X,2);
 phi1prim = 260.5+200*t;
 phi2prim = 240.5+200*t;
 phiprim = zeros(2,N); 
 phiprim(1,:) = phi1prim;
 phiprim(2,:) = phi2prim;
 gamma=0; 
 sigma = 0.04;
else
 %we consider the hyperbolic chirp plus the linear chirp
 load('signal_hyp_lin.mat');
 x0 = sig;
 N = length(x0);
 Nfft = N;
 t  = (0:N-1)/N;
 gamma = 0;
 sigma = 0.03;
end 


Nfft = N; 
%One compute the STFT 
[STFT,~,~,~,~] = sst2(x0,sigma,Nfft,gamma);

%then the spectrogram
spect = (abs(STFT)/N).^2; %% compute the discrete spectrogram, be careful of to the normalization

%% approximation of the spectrogram, we consider the approximation (9) of the paper
 M_0 = 20; %order of truncation
 Cp =  sigma/(Nfft*sqrt(2))* exp(-pi*((-M_0:M_0).^2/(2*(sigma*Nfft).^2)))';

 %% computation of the localisation of the localisation of the IFs 

 V = exp(2*pi*1i*((0:Nfft-1)'/Nfft*(-M_0:M_0)));
 D = diag(Cp); 
 W = V*D;
 L = zeros(2*M_0+1,N);
 for n = 1:N
  L(:,n) = W\spect(:,n);
 end 

 %annihilating filter technique
 FoundDiracsLocations = zeros(P,N);
 Amp = zeros(P,N);
 
 jump = 20; %this value has to be adapted to the frequency resolution

 for n = 1:N
  y = L((M_0+1)-(P-1):(M_0+1)+P,n);
  [h] = YW(y);
  
  %computation of the Diracs Location 
  TZroot = roots(h(end:-1:1));
  X =(angle(TZroot)/(2*pi)) > 0;
  Y = angle(TZroot)/(2*pi).*X+(1+angle(TZroot)/(2*pi)).*(1-X);
  FoundDiracsLocations(:,n) = sort(Y*Nfft,'descend');
 
  %computation of the amplitude
  X = FoundDiracsLocations(:,n);
  W = exp(-2*1i*pi*(0:P-1)'*X'/Nfft);
  Amp(:,n) = real(W\L((M_0+1):(M_0+1)+P-1,n));  
 end
 FoundDiracsLocations_amp= FoundDiracsLocations;

 %elimination of the point that correspond to jumps in the IF estimator or
 %to a too small amplitude
 for n=1:N    
  if (n > 1)&&(n < N)  
   for p = 1:P   
    %The Dirac Location cannot jump in neighboring locations
    if (abs(FoundDiracsLocations(p,n-1)- FoundDiracsLocations(p,n)) + ... 
       abs(FoundDiracsLocations(p,n+1)- FoundDiracsLocations(p,n)) > jump)
     FoundDiracsLocations_amp(p,n) = 0; 
    end
   end
  end
  %for each n, we only consider the amplitudes that are larger than some theshold
  val_thresh = 10^(-2)*max(abs(Amp(:,n)));
  index =  (abs(Amp(:,n)) > val_thresh);
  FoundDiracsLocations_amp(:,n) = (FoundDiracsLocations_amp(:,n).*index);
 end

 % we compute the indices where non-zero Dirac locations are found
 Index_nonzeros= cell(P,3);
 for p = 1: P
  X = (1:N).*(FoundDiracsLocations_amp(p,:) > 0);
  X = X(X>0);
  Index_nonzeros{p,1} = X;
  Index_nonzeros{p,2} = FoundDiracsLocations_amp(p,X);
  Index_nonzeros{p,3} = Amp(p,X);
 end    

 %We use polynomial interpolation to fill the holes 
 FoundDiracsLocation_rec = zeros(P,N);
 Amp_rec = zeros(P,N);
 remove_mode = ones(1,P); %arrey used to tell which mode to discard in the end
 for p = 1:P
  X = Index_nonzeros{p,1};
  if (length(X) > N/2) % we only consider the ridges with sufficiently many points to perform the interpolation
   Y = Index_nonzeros{p,2};
   FoundDiracsLocation_rec(p,:) = pchip(X,Y,1:N);
   Y = Index_nonzeros{p,3};
   Amp_rec(p,:)   = pchip(X,Y,1:N);
   remove_mode(p) = 0; 
  end 
 end

 %elimination of the modes that are associated with negative amplitude for a certain duration 
 for p = 1:P
  Y = Amp_rec(p,:);
  Y = Y(Y<0);
  if (length(Y) > N/10) % the amplitude is oscillating 
   remove_mode(p) = 1;
  end    
 end
 
 %output figures 
 if cas <= 2
  figure
  imagesc((0:N-1)/N,(0:Nfft-1)/Nfft*N,spect);
  set(gca,'ydir','normal');
  hold on;
  time_index = 100:N-100;
  if (cas == 1)
   ylim([200 360])
   xlim([t(time_index(1)) t(time_index(end))])
  elseif cas == 2
   ylim([250 450])
   xlim([t(time_index(1)) t(time_index(end))])  
  end
  xlabel('time','FontSize',20);
  ylabel('frequency','FontSize',20);
  ax = gca;
  ax.FontSize = 20;
 end
 
 %computation ot the final ridges
 index = (1:P).*(remove_mode == 0);
 index = index(index > 0);
 FoundDiracsLocations_new = FoundDiracsLocation_rec(index,:);

 if cas <= 2
  time_index = 100:N-100;
  %figure;
  %subplot(2,1,1), 
%   if cas == 1
%    plot(t(time_index),phi1prim(time_index),'-o',t(time_index),phi2prim(time_index),'-c',t(time_index),phi3prim(time_index),'-s','Linewidth',2);
%    ylim([200 360])
%    xlim([t(time_index(1)) t(time_index(end))])
%    xlabel('time','FontSize',20);
%    ylabel('frequency','FontSize',20);
%    title('ground truth IF','FontSize',20);
%    ax = gca;
%    ax.FontSize = 15;
%   elseif cas == 2
%    plot(t(time_index),phi1prim(time_index),'-o',t(time_index),phi2prim(time_index),'-c','Linewidth',2);  
%    xlim([t(time_index(1)) t(time_index(end))])
%    xlabel('time','FontSize',20);
%    ylabel('frequency','FontSize',20);
%    title('ground truth IF','FontSize',20);
%    ax = gca;
%    ax.FontSize = 15;
%   end     
%   subplot(2,1,2),
%   hold on;
%   length(index)
%   for p = 1:length(index)
%    mean_error_FRI(p)      = sqrt(mean((phiprim(p,time_index)-FoundDiracsLocations_new(p,time_index)).^2))
%   end
  
  for k=1:length(index)
   plot(t(time_index),FoundDiracsLocations_new(k,time_index),'Linewidth',2);
   xlim([t(time_index(1)) t(time_index(end))])
   xlabel('time','FontSize',20);
   ylabel('frequency','FontSize',20);
   title('spectrogram + estimated IFs','FontSize',20);
   ax = gca;
   ax.FontSize = 20;
  end
  if (cas == 1)
   ylim([200 360])
   xlim([t(time_index(1)) t(time_index(end))])
  end
  hold off;
 end
 if cas == 3
  time_index = 30:N-70;
  figure
  imagesc(t(time_index),(0:Nfft-1)/Nfft*N,spect(:,time_index)); 
  set(gca,'ydir','normal');
  hold on;
  for k=1:length(index)
   plot(t(time_index),FoundDiracsLocations_new(k,time_index),'Linewidth',2);
   ylim([70 400])
   xlabel('time','FontSize',20);
   ylabel('frequency','FontSize',20);
   title('estimated IFs','FontSize',20);
   ax = gca;
   ax.FontSize = 20;
  end
  hold off;
 end
end