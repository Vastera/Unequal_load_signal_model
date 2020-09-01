%% CopyRight@ vastera@163.com
% Code for the paper "A vibration signal model of planetary gearboxes under
% uneven load sharing among planets" author: Haoqun Ma, Zhipeng Feng.
% 《行星轮不均载的行星齿轮箱振动信号模型》论文信号模型代码 作者：马浩群，冯志鹏
%% %%%%%%%%%%%%%% Paramters 模型参数设置 %%%%%%%%%%%%%%
% input torque 输入负载大小 单位Nm
T_s=3e2;
% position error of i-th planet 行星轮轴的周向位置偏差
% vector length denotes the planet number, each element denotes the corresponding circumferential error 
% 向量长度代表行星轮数，每一个位置的值分别对应一个行星轮周向偏差大小 单位是弧度
% epsilon_i=[0,0,0,0/180*pi,0,0];%6个行星轮的正常情况
epsilon_i=[0,0.5/180*pi,0,0];%4个行星轮中第二个有误差
%% Intialization
fs=10000;% sampling frequency
t=1/fs:1/fs:10;% time seiries
N=61;%ring gear tooth number 
%sun gear tooth number is 
f_c=3;%carrier frequency
R_p=35*0.9e-3/2;% gear module is 0.9e-3
R_r=108*0.9e-3/2;% gear module is 0.9e-3
R_s=36*0.9e-3/2;% gear module is 0.9e-3
F=T_s/300;
%% unequal load sharing case
% perfect case
% epsilon_i=[0,0,0,0];% position error of i-th planet
M=length(epsilon_i);%planet number
a=0.1*2*pi*R_r/M;
%% configuration parameters
i=1:1:M;% planet sequential order
Theta_i=2*pi*(i-1)./M;%nominal position
phi=0.05*pi/180;% phase difference between sun and planet
%% 
e_i=2*(R_r-R_p)*sin(epsilon_i/2).*cos(epsilon_i/2);
k_e=1/(1/3.2e7+1/(1.67e8+2.57e8));
%% Calculate the load sharing coefficients L_i using tranlational analogy
L_i = Load_sharing_coef(epsilon_i,T_s,R_p,R_r,R_s,k_e);
%% when certain planets are unloaded
f_n=200+100*sum(L_i~=0);% f_n increases with the contacting planet number
%% Transfer path effect
eta=exp(-2*R_p/a);%the amplitude factor between planet-ring and planet-sun
f_m=N*f_c;
f_s=f_m/36+f_c;
x=zeros(1,length(t));
l=zeros(M,length(t));
sigma_i=zeros(M,length(t));
for j=1:M
    T_ri=f_m*(t-(Theta_i(j)+epsilon_i(j))/(2*pi*f_c));% actual time sequence of planet-ring
    T_si=f_m*(t-(Theta_i(j)+epsilon_i(j))/(2*pi*f_s-2*pi*f_c))+phi/2/pi;% time sequence of planet-sun
    l(j,:)=Transfer_length(t,M,f_c,R_r,epsilon_i(j),Theta_i(j));
    sigma_i(j,:)=F*exp(-l(j,:)/a);% sigma_i is the time-varying transfer path effect
    x=L_i(j)*sigma_i(j,:).*(III(T_ri)+eta*III(T_si))+x;
end
%% Natural vibration
lambda=exp(-150*t).*sin(2*pi*f_n*t);%normal case
x1=conv(x,lambda);x1=x1(1:length(t));
%% Fourier spectrum with natural frequency
[Amplitude,f]=MyFFT(x1,fs);
[~,f_c0]=III(f_m*t);
f_c0=fs/f_c0/N;% the actual carrier frequency in signal
figure('Name','Fourier spectrum');
plot(f/f_c0,Amplitude,'b');
[Amp_natura1,f]=MyFFT(lambda,fs);
% Natural frequency envelope in spectrum
hold on;
plot(f/f_c0,Amp_natura1*F*20,'g-.');
% transfer path envelope in spectrum
sigma_0=zeros(1,length(t));
sigma_0(t<=Theta_i(2))=F*exp(-2*pi*R_r*f_c0*abs(t(t<=Theta_i(2))-Theta_i(2)/2)/a);
plot(f/f_c0,1.5e2*MyFFT(sigma_0.*III(f_m*t),fs),'r:');
legend('Synthesized signal','Natural vibration','Tranfer path effect','Location','northeast');
% figure property settings
switch M
    case 3
        xlim([110 140]);
    case 4
        xlim([170 200]);
    case 5
        xlim([100 300]);%ylim([0 0.03]);
    case 6
        xlim([100 300]);ylim([0 0.02]);
end
ylabel('Amplitude');xlabel('Carrier order');
FigureSettings;
%% Time domain figure
figure('Name','Time domain');
plot(t*f_m,x1,'b');hold on;
xlim([0,80]);
xlabel('Gear meshing interval');
ylabel('Amplitude');
FigureSettings;

