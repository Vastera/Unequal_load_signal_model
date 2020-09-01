%Calulate the load sharing coefficients of each planet, referring to the Singh's paper A Closed-Form ... Using a Translational Analogy. 
function [L,P_support] = Load_sharing_coef(epsilon,T_s,R_p,R_r,R_s,k_e)
%% CopyRight vastera@163.com mahaoqun ÂíºÆÈº
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% epsilon:        Type: vector with length of planet number
%            desciption:   Postion angular errors of planets, differing from the tangential error e.
% T_s:    Type: double
%             description: input torque to sun gear
% R_p:     Type: double
%            description: radius of the planet pitch circle
% R_r:     Type: double
%            description: radius of the ring gear pitch circle
% R_s:     Type: double
%            description: radius of the sun gear pitch circle
%%%%%%%%%%%%%%%%%%%%%%%%%%% Output:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P_support: Type: a vector of 3*1
%             description: the three support planet at first for input load
% L         Type: vector with the same length of e
%            description: load sharing coefficients of planets
% Translational analogy 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize parameters%%%%%%%%%%%%%%%%%%%%%
M=length(epsilon);
i=1:1:M;
Theta=2*pi*(i-1)/M;
r_b=R_s+R_p;
e=2*(R_r-R_p)*sin(epsilon/2).*cos(epsilon/2);
z=e;
x=-r_b*sin(Theta+epsilon);
y=r_b*cos(Theta+epsilon);
W=2*T_s/R_s;
%% %%%%%%%%%%%%% Determine initial configuration%%%%%%%%%%%%%%%%%
% 1. The first contact point is : 
[~,P_J]=max(z);
% Define a line  tangent to the circle on plane (parallele with the x-y plane), anther point  also on the line 

x_J=x(P_J);y_J=y(P_J);
x_Q=(x_J^3+x_J*y_J^2-sqrt(r_b*x_J^2*y_J^2+r_b*y_J^4))/(x_J^2+y_J^2);% another point Q having the distance of R_b with x_J 
y_Q=-x(P_J)/y(P_J)*(x_Q-x(P_J))+y(P_J);
z_Q=z(P_J);
%% 2. The second contact point 
%Define the plane  by three points :
% plane_eqn definition
% tranverse all i~=J, find i to make 
K_options=[1:P_J-1,P_J+1:M];% tranverse the set without J
P_K = search_plane(K_options, P_J, [x_Q,y_Q,z_Q], x, y, z);
%% 3. The third contact point M
% Define the plane  by the three points 
% tranverse all i~=J && i~=K to make 
M_options=K_options(K_options~=P_K);
P_M=search_plane(M_options,P_J,[x(P_K),y(P_K),z(P_K)],x,y,z);
% Verify whether the load  application point  (origin point) is within the triangle 
% Satisfy the conditions: make each of three pairs of vectors has the same direction.
% if any of three condition not satisfied, then start from the closed edge to repeat the step 3.
Candidates=M_options(M_options~=P_M);% chose from the residual rather than J , K, M points
P_support3=find_3_support(Candidates,P_J, P_K, P_M,  x, y, z);
P_J=P_support3(1);P_K=P_support3(2);P_M=P_support3(3);
% Calculate the forces
% Calculate the forces of three points () by applying the force and moment balance equations
A=[1,1,1;
   y(P_J),y(P_K),y(P_M);
   x(P_J),x(P_K),x(P_M)];
BB=[W;0;0];
F0=A\BB;% find the solution of linear equations
P_support0=[P_J,P_K,P_M];
if min(F0)/min(F0(F0~=min(F0)))<1e-3% if the third spring carries no load
    [~,spare_spring]=min(F0);
    F0(spare_spring)=[];% only two spring carrying load initially
    P_support0(spare_spring)=[];
end
%% %%%%%%%%%%%%%%%%%%%% Load sharing coefficients when support planet >3 %%%
% correponding spring deflection 
F=zeros(1,M);
F(P_support0)=F0;
delta=F/k_e;
% update the spring coordinates 
zd=z-delta;
% If three spring support the plane 
% Define the new plane  by three points ,,
% Calculate the plane coefficients 
xx=[x(P_J),x(P_K),x(P_M)];
yy=[y(P_J),y(P_K),y(P_M)];
if length(P_support0)==3% three spring support load
    zz=[zd(P_J),zd(P_K),zd(P_M)];
else % two spring support load
    zz=[zd(P_support0),mean(zd(P_support0))];
end
    [a,b,c,d]=Plane_eqn(xx,yy,zz);
% Calculate the coordinate in plane for all planets
zd=(-d-a*x-b*y)/c;
% For all H candidates 
H_options=1:M;
H_options=H_options(H_options~=P_J & H_options~=P_K & H_options~=P_M);
P_support=P_support3;
i=1;
while ~isempty(H_options) && i<=length(H_options)
    H=H_options(i);
% if H point protrude the plane 1
    if z(H)>zd(H)
        P_support=[P_support,H];
% Chose  from  to satisfy load application point C is within triangel 
% the other point than  of  is denoted by 
        for G_options=transpose(nchoosek(P_support3,2))%chose 2 from the points J K M
            if IsinTriangle(H,G_options(1),G_options(2),x,y,z)
                G1=G_options(1);G2=G_options(2);B=P_support3(P_support3~=G1 & P_support3~=G2);
                break;
            end
        end
        eta=((x(G2)-x(G1))*(y(B)-y(G1))-(x(B)-x(G1))*(y(G2)-y(G1)))/((x(G2)-x(G1))*(y(H)-y(G1))-((x(H)-x(G1))*(y(G2)-y(G1))));
        A1=(eta-1)*y(G2)+y(B)-eta*y(H);
        A2=(1-eta)*y(G1)-y(B)+eta*y(H);
        A3=y(G2)-y(G1);
        A4=eta*(y(G2)-y(G1));% Here is corrected (A4 in origin paper is wrong)
        Y=[A1,A2,A3,-A4]*[z(G1);z(G2);z(B);z(H)];
        % match the  order to the , for example ,then A's order:
        ind_order=[find(P_support3==G1,1),find(P_support3==G2,1),find(P_support3==B,1)];
        [~,ind_order]=sort(ind_order);
        sub_A=[A1,A2,A3];sub_A=sub_A(ind_order);
        A=[A,[1;y(H);x(H);zeros(size(A,1)-3,1)]];
        A=[A;sub_A,zeros(1,size(A,2)-4),-A4];
        BB=[BB;k_e*Y];
        % update the  for all spring
        F0=A\BB;
        F=zeros(1,M);
        F(P_support)=F0;
        delta=F/k_e;
        zd=z-delta;
        xx=[x(P_J),x(P_K),x(P_M)];
        yy=[y(P_J),y(P_K),y(P_M)];
        zz=[zd(P_J),zd(P_K),zd(P_M)];
        [a,b,c,d]=Plane_eqn(xx,yy,zz);
        zd=(-d-a*x-b*y)/c;
        H_options=H_options(H_options~=H);% remove H from H_options
        i=0;
    end
    i=i+1;
end
if size(A,1)>3% more than 3 planet support
    F=A\BB;
    L0=F/W;
else
    L0=F0/W;
end
%% adjust the output arguments' order
L=zeros(1,M);
L(P_support)=L0;
P_support=sort(P_support);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Support function %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% plane parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a,b,c,d]=Plane_eqn(x, y, z)
a=(y(2)-y(1))*(z(3)-z(1))-(y(3)-y(1))*(z(2)-z(1));
b=(x(3)-x(1))*(z(2)-z(1))-(x(2)-x(1))*(z(3)-z(1));
c=(x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1));
d=x(2)*(y(1)*z(3)-z(1)*y(3))+x(1)*(y(3)*z(2)-z(3)*y(2))+x(3)*(z(1)*y(2)-y(1)*z(2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% search the third point given two points
function P_K = search_plane(K_options, P_J,P_Q, x, y, z)
% NOTICE: P_J is the order number in all candidates points( planet number
% 1:M) BUT the P_Q is the vector formed by the Cartesian coordinates [x_Q, y_Q, z_Q]
% OUPUT: P_K is a order number of all planet candinates
P_K=K_options(1);% the first candidate
x_Q=P_Q(1);y_Q=P_Q(2);z_Q=P_Q(3);%the seconde point coordinates
xx=[x(P_J),x_Q,x(P_K)];
yy=[y(P_J),y_Q,y(P_K)];
zz=[z(P_J),z_Q,z(P_K)];
[a,b,c]=Plane_eqn(xx, yy, zz);
alpha=acos(abs(c)/sqrt(a^2+b^2+c^2));% the baseline alpha value of first candidate
for j=K_options(2:end)
    xx=[x(P_J),x_Q,x(j)];
    yy=[y(P_J),y_Q,y(j)];
    zz=[z(P_J),z_Q,z(j)];
    [a,b,c]=Plane_eqn(xx, yy, zz);
    if acos(abs(c)/sqrt(a^2+b^2+c^2))<alpha% if the i-th candidate is smaller than current
        P_K=j;
    end
end
end
%%%%%%%%%%%%%%%%%% whether the point C is within triangle %%%%%%%%%%%%%
function Inside=IsinTriangle(P_1,P_2,P_3,x,y,z)% judge whether C is within trigangle P_1,P_2,P_3
xx=[x(P_1),x(P_2),x(P_3)];
yy=[y(P_1),y(P_2),y(P_3)];
zz=[z(P_1),z(P_2),z(P_3)];
[~,~,c,d]=Plane_eqn(xx,yy,zz);
C=[0,0,-d/c];% C is on the plane JKM, and x=0,y=0,0*x+0*y+c*z+d=0;
Inside=(IsSameSide(P_1,P_2,P_3,C,x,y,z) && IsSameSide(P_2,P_3,P_1,C,x,y,z) && IsSameSide(P_3,P_1,P_2,C,x,y,z));
end
%%%%%%%%%%%%%%%% used in function IsinTriangle %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SameSide=IsSameSide(P_1,P_2,P_3,C,x,y,z)
% judge whether the two  points  at the same side of line
P_A=[x(P_1),y(P_1),z(P_1)];
P_B=[x(P_2),y(P_2),z(P_2)];
P_C=[x(P_3),y(P_3),z(P_3)];
AB=P_B-P_A;
AC=P_C-P_A;
AP=C-P_A;
v1=cross(AB,AC);
v2=cross(AB,AP);
SameSide=(dot(v1,v2)>0);
end
%%%%%%%%%%%%%%%%%%%%%% find three support springs for the input load W
function P_support=find_3_support( Candidates,P_J, P_K, P_M, x, y, z)
xx=[x(P_J),x(P_K),x(P_M)];
yy=[y(P_J),y(P_K),y(P_M)];
zz=[z(P_J),z(P_K),z(P_M)];
[~,~,c,d]=Plane_eqn(xx,yy,zz);
C=[0,0,-d/c];% C is on the plane JKM, and x=0,y=0,0*x+0*y+c*z+d=0;
if ~IsinTriangle(P_J,P_K,P_M,x,y,z) % if C is not in triangel JKM
    if ~IsSameSide(P_J,P_K,P_M,C,x,y,z)% if C is outside of JK
        Candidates=Candidates(Candidates~=P_M);% remove M from Candidates
        P_M=search_plane(Candidates,P_J,[x(P_K),y(P_K),z(P_K)],x,y,z);% continue searching
    elseif ~IsSameSide(P_K,P_M,P_J,C,x,y,z)% if C is outside of KM
        Candidates=Candidates(Candidates~=P_J);% remove J from Candidates
        P_J=search_plane(Candidates,P_K,[x(P_M),y(P_M),z(P_M)],x,y,z);% continue searching
    else% if C is outside of MJ
        Candidates=Candidates(Candidates~=P_K);% remove J from Candidates
        P_K=search_plane(Candidates,P_M,[x(P_J),y(P_J),z(P_J)],x,y,z);% continue searching
    end
    P_support=find_3_support(Candidates,P_J,P_K,P_M,x,y,z);% Recursion of judging and searching
else
    P_support=[P_J,P_K,P_M];
end
end
