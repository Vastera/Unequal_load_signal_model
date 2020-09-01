function l=Transfer_length(t,M,f_c,R_r,epsilon_i,Theta_i)
% Copyright@ vastera@163.com
% General introduction:Calculate the length of the time varing transfer path
%% ====================== INPUT ========================
% t:          Type:vector
%                           t description:input time series
% M:          Type: integer
%                           M description: planet numbers
% f_c:        Type:number
%                           f_c description:carrier Frequency
% R_r:        Type: integer
%                           R_r description: ring gear radius
% ---------------------OPTIONAL:
% optional arg:              Type:
%                            description:
%% ====================== OUTPUT =======================
% l:          Type:vector with the same length of t
%                           l description: the time-varying length of transfer path
%% =====================================================
l=Inf*ones(1,length(t)); % otherwise occasion
%% Compare first left boundary with zero
left= -1/(2*M*f_c)+(epsilon_i+Theta_i)/2/pi/f_c;
right= 1/(2*M*f_c)+(epsilon_i+Theta_i)/2/pi/f_c;
if left>0% if left bound is positive
    index=(t>=left & t<right);
else% if left bound is negative
    index=(t<right);
end
l(index) = 2*pi*R_r*f_c*abs(t(index)-(epsilon_i+Theta_i)/2/pi/f_c);
%% Residual n occassions 
n=1;
left= n/f_c-1/(2*M*f_c)+(epsilon_i+Theta_i)/2/pi/f_c;
right= n/f_c+1/(2*M*f_c)+(epsilon_i+Theta_i)/2/pi/f_c;
while right<=max(t)%% while the right bound is lower than the max time
    index=(t>=left & t<right);
    l(index) = 2*pi*R_r*f_c*abs(t(index)-n/f_c-(epsilon_i+Theta_i)/2/pi/f_c);% the first piece function.
    n=n+1;
    left= n/f_c-1/(2*M*f_c)+(epsilon_i+Theta_i)/2/pi/f_c;
    right= n/f_c+1/(2*M*f_c)+(epsilon_i+Theta_i)/2/pi/f_c;
end
end