function [ResonantFrequency]=Resonant_frequency(M,k,varargin)
% Copyright@ vastera@163.com
% General introduction: Calculate the resonance frequency (circular frequency) of the system with given parameters M,k,c.
%% ====================== INPUT ========================
% M:          Type: a number
%                           M description: Mass
% k:          Type: a number
%                           k description: stiffness
% ---------------------OPTIONAL:
% optional arg:  c           Type: number
%                            description: damper
%% ====================== OUTPUT =======================
% ResonantFrequency:          Type: a number
%                           ResonantFrequency description: the system's natural frequency (resonance frequency)
%% =====================================================
narginchk(2, 3);
if nargin<=2
    ResonantFrequency=sqrt(k/M);
end
if nargin==3
    c=varargin{1};
    omega_n=sqrt(k/M);
    epsilon=c/(2*M*omega_n);
    ResonantFrequency=omega_n*sqrt(1-epsilon^2);
end
end

