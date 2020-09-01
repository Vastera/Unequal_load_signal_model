function [ Amp,F ] = MyFFT( Sig,Fs )
%MYFFT  calculation of Fourier transform 
%   input: 
%          Sig is input signal;
%          Fs is sampling frequency
%   Output:
% %          amp is the normalized amplitude of FFT's output
%          f is the array of output frequency 
N=length(Sig);
Amp=abs(fft(Sig))/N*2;
Amp=Amp(1:fix(N/2));
F=Fs/N*(1:fix(N/2));
end

