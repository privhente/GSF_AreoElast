% function [ft_peak,f,P] = fft_x(fig,u,t,range)
    L  = length(u); % Length of signal
    T  = t(2)-t(1);      % Sampling period
    Fs = 1/T;            % Sampling frequency
    NFFT = 2^nextpow2(L);
%     Y  = fft(u,NFFT);
    Y  = fft(u);
    
    P2 = abs(Y/L);
    P1 = P2(1:fix(L/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L *2*pi;
    
%     f = Fs/2*linspace(0,1,NFFT/2+1);
%     SSAS = 2*abs(Y(1:NFFT/2+1));
    
    [val,inzm] = max(P1);
    P1 = P1/val;
    [c,inzf,val] = find(f>=range(1) & f<=range(2));
    
%     figure(fig); hold on; grid on; %axis([0 5 0 1.1*max(P1)]);
%     title('Single-Sided Amplitude Spectrum');
%     xlabel('\omega rad/s'); ylabel('|fft(u)|');
%     rgb_color = RGB_Color;
%     plot(f(inzf),P1(inzf),'-','linewidth',1.0); 

    [val,inz] = findpeaks(P1(inzf));
    f = f(inzf);
    P = P1(inzf);
%     P = P/max(P);
    ft_peak = f(inz);
return