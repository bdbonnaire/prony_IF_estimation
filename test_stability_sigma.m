function test_stability_sigma(P)
 %INPUT 
 % P   :  number of mode to be detected
%  folder = './';
%  addpath(folder);
%  addpath(strcat([folder 'synchrosqueezedSTFT']));
%  addpath(strcat([folder 'tools']));
%  addpath(strcat([folder 'FRI_lib']));
% 

 N  = 1024;         %% signal length
 t  = (0:N-1)/N;

 freq1 = 240.5;
 freq2 = 220.5;
 A = 2;
 X(:,1) = A*exp(2*1i*pi*freq1.*t);
 B = 2;
 X(:,2) = B*exp(2*1i*pi*freq2.*t);
 x0 = sum(X,2);

 phiprim(1,:) = 240.5*ones(1,N);
 phiprim(2,:) = 220.5*ones(1,N);

 Nfft = N; 
 sigma = [0.024 0.028 0.032 0.036 0.040 0.044 0.048];
 gamma = 0;
 
 [STFT,FSST,~,omega,~] = sst2(x0,sigma(6),Nfft,gamma);
 %then the spectrogram
 spect = (abs(STFT)/N).^2;
 figure
 imagesc((0:N-1)/N,(0:Nfft-1)/Nfft*N,spect)
 time_index = 100:N-100;
 ylim([200 260])
 xlim([t(time_index(1)) t(time_index(end))])
 xlabel('time','FontSize',20);
 ylabel('frequency','FontSize',20);
 set(gca,'ydir','normal');
 ax = gca;
 ax.FontSize = 20;
 
 for k = 1:length(sigma)
  %One compute the STFT 
  [STFT,FSST,~,omega,~] = sst2(x0,sigma(k),Nfft,gamma);

  %then the spectrogram
  spect = (abs(STFT)/N).^2; %% compute the discrete spectrogram, be careful of to the normalization

  %We extract the ridges on the STFT and  
  [Cs0,] = exridge_mult(abs(STFT),2, 0, 10);
  [Cs,]  = exridge_mult(abs(FSST),2, 0, 10);
 
  %we compute the reassignment operator omega on the ridge of the FSST
  for q = 1:N
   for p = 1:2
    freq(p,q) = omega(Cs(p,q),q);
   end
  end
  
  %We compute the Dirac pulses using an FRI approach
  %% approximation of the spectrogram, we consider the approximation (9) of the paper
  M_0 = 20; %order of truncation
  Cp =  sigma(k)/(Nfft*sqrt(2))* exp(-pi*((-M_0:M_0).^2/(2*(sigma(k)*Nfft).^2)))';

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
  time_index = 100:N-100;
  for p = 1:P
   Y = Amp_rec(p,time_index);
   Y = Y(Y<0);
   if (length(Y) > length(time_index)/10) % the amplitude is oscillating 
    remove_mode(p) = 1;
   end    
  end
  
  index = (1:P).* (remove_mode == 0);
  index = index(index > 0);
  FoundDiracsLocations_new = FoundDiracsLocation_rec(index,:);
 
  for p = 1:2
   %error of IF estimation with the spectrogram
   mean_error_spec(k,p)     = sqrt(mean((phiprim(p,time_index)-Cs0(p,time_index)-1).^2));
   mean_error_FSST(k,p)     = sqrt(mean((phiprim(p,time_index)-Cs(p,time_index)-1).^2));
   mean_error_reassign(k,p) = sqrt(mean((phiprim(p,time_index)-freq(p,time_index)).^2));
 
   %error using FRI
   mean_error_FRI(k,p)      = sqrt(mean((phiprim(p,time_index)-FoundDiracsLocations_new(p,time_index)).^2));
  end
 end  

 %output figures 
 figure;
%  subplot(2,1,1),
%  plot(sigma,mean_error_spec(:,1),sigma,mean_error_FSST(:,1),'--',sigma,mean_error_reassign(:,1),'-.',...
%  sigma,mean_error_FRI(:,1),'-o','LineWidth',2);
%  xlim([0.024 0.048]);
%  xlabel('\sigma','FontSize',20);
%  ylabel('estimation error','FontSize',20);
%  legend('spectrogram ridge','FSST ridge','$\widehat{\omega_f}$','FRI','Interpreter','latex');
%  set(gca,'ydir','normal');
%  ax = gca;
%  ax.FontSize = 20;
%  ylim([0 20])
%  subplot(2,1,2)
 plot(sigma,mean_error_spec(:,2),sigma,mean_error_FSST(:,2),'-s',sigma,mean_error_reassign(:,2),'-d',...
 sigma,mean_error_FRI(:,2),'-o','LineWidth',2,'MarkerSize',10);
 xlim([0.024 0.048]);
 xlabel('\sigma','FontSize',20);
 ylabel('estimation error','FontSize',20);
 legend('IF-SR','IF-FSSTR','IF-FSST-OG','Prony');
 set(gca,'ydir','normal');
 ax = gca;
 ax.FontSize = 20;
 ylim([0 20])
end 