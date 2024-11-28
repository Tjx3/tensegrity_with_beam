clc;
close all;

%% add E I l A
E=1e5;
I=1e4;
l=1;
A=1e-7;

%% beam with v \theta
K_v=[12 6 -12 6;
    6 4 -6 2;
    -12 -6 12 -6;
    6 2 -6 4]

[U,S,V] = svd(K_v);
r=rank(K_v)                      % rank of (A_1g)
U1=U(:,1:r)
U2=U(:,r+1:end);        % U1 is C(A_1g); U2 is N(A_1g') mechanism mode
% S1=S(1:r,1:r);                      % S1 is singular value of A_1g
V1=V(:,1:r);V2=V(:,r+1:end);        % V1 is C(A_1g'); V2 is N(A_1g) self stress mode


%% truss 
K_u=[-1;1]*[-1 1]
[U,S,V] = svd(K_u);
r_u=rank(K_u)                      % rank of (A_1g)
U1=U(:,1:r_u)
U2=U(:,r_u+1:end);        % U1 is C(A_1g); U2 is N(A_1g') mechanism mode
% S1=S(1:r,1:r);                      % S1 is singular value of A_1g
V1=V(:,1:r_u);V2=V(:,r_u+1:end);       % V1 is C(A_1g'); V2 is N(A_1g) self stress mode

%%beam
K=[1 0 0 -1 0 0
     0 12 6 0 -12 6;
     0 6 4 0 -6 2;
    -1 0 0 1 0 0
    0 -12 -6 0 12 -6;
    0 6 2 0 -6 4]
[U,S,V] = svd(K);
r=rank(K)                      % rank of (A_1g)
U1=U(:,1:r)
U2=U(:,r+1:end);        % U1 is C(A_1g); U2 is N(A_1g') mechanism mode
% S1=S(1:r,1:r);                      % S1 is singular value of A_1g
V1=V(:,1:r);V2=V(:,r+1:end);       % V1 is C(A_1g'); V2 is N(A_1g) self stress mode

%% multipule members: structures;