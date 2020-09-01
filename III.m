function [y,interval]=III(t,varargin)
% Copyright@ vastera@163.com author£ºmahaoqun 
% General introduction:generate a Dirac delta comb function corresponding to input time series at each integer point
%% ====================== INPUT ========================
% t:          Type:vector
%                          t description:input time series (mayber scaled version ) 
% ---------------------OPTIONAL:
% optional arg:  round_type            Type: an integer
%               description:a flag denoting weather to round the impulsive interval,
%                             1 is to find the closest point the nominal position,
%                             2 (default) is to make sure all interval are identical
%% ====================== OUTPUT =======================
% y:          Type: a vector with the same length of t
%                           y description: the output dirac comb series
% interval:  Type: an integer
%                         interval description: actual interval (points not time) of impulses
%% =====================================================
%% check the input arguments
narginchk(1, 2);
round_type=2;
if nargin>=2
    round_type = varargin{1};
end
%% shift towards right by 0.5 (add 0.5)
Reference=t+0.5;
%% floor(t) to get the integer reference
Reference=floor(Reference);
%% minus the integer reference to get the difference to each integer point
Difference=t-Reference;
%% get pairs of points closest to the reference integer point though capturing the upward edges
Difference1=circshift(Difference,-1);% shift the valuses
Index=find((Difference<=0) &(Difference1>0));
%% compare every pair of points to get the every closest points close to the integers
y=zeros(1,length(t));
if round_type==1 % the default round_type
    interval=Index(2)-Index(1);
    for j=1:length(Index)
        if abs(Difference(Index(j)))<=abs(Difference( min([Index(j)+1,length(Difference)]))) % min here in case of exceed Difference length
            y(Index(j))=1;
        else
            y(Index(j)+1)=1;
        end
    end
elseif round_type==2 % Make sure every interval to be identical
    interval=round((Index(end)-Index(1))/(length(Index)-1));
    for j=1:length(Index)
        y(min([Index(1)+interval*(j-1),length(t)]))=1;
    end
end
end

