function [C,B_pc,Dz,H,eigenvector]= optrobust_SIM(P) 
%  optrobust_SIM (P): Optimal Robust Stabilization (This is a simplified version)
%
%    optrobust_SIM (P)  is to design a discrete controller for a given discrete system P such that
%    the stability margin b_pc is maximized.
% 
%    Input "P" is the given system (created with either TF, ZPK, SS, or FRD).
%
%    [C,B_OPT]=OPTROBUSTDS(P) returns the stabilizing controller C and the 
%    maximal stability margin B_OPT.
%
% Example : 
% P=tf([0 2 1],[1 4 4],1)
% [C,B_pc]= optrobust_SIM(P)
%
% C =
% 
%  -6.678 z - 8.542
%  ----------------
%  4.576 z + 2.288
% 
% Sample time: 1 seconds
% Discrete-time transfer function.
%
%
% B_pc =
%
%    0.1018
%
syms c;
[tempB,tempA]=tfdata(P);
az=tempA{1,1};
az_i= fliplr(az);
bz=tempB{1,1};
bz_i= fliplr(bz);
m=length(az); n=m-1;%n is the order of az
[Dz,Dzi] = Spectral_Factorization(az,bz);
%%% Matrix construction
L_a=zeros(n,n);L_b=zeros(n,n);L_d=zeros(n,n);
U_a=zeros(n,n);U_b=zeros(n,n);U_d=zeros(n,n);
for i=1:n
    L_a(i:end,i)=az(1:m-i); U_a(i,i:end)=az_i(1:m-i);
    L_b(i:end,i)=bz(1:m-i); U_b(i,i:end)=bz_i(1:m-i);
    L_d(i:end,i)=Dz(1:m-i); U_d(i,i:end)=Dzi(1:m-i);
end
J=fliplr(eye(n));
H=J*inv(L_d)*J*[U_b*J -U_a*J]*inv([L_a L_b; U_a U_b])*[L_d; U_d];
%%%%%%%%%%%%%%%%%%%%%
% Finding the largest eigenvale(hankel norm) and its corresponding
% eigenvector
[eigenvector,eigvalues]=eig(H);
eigvalues=abs(eig(H));
[value,position]=max(eigvalues);
eigenvector=eigenvector(:,position)/eigenvector(1);
%%%%%%%%%%%%%%%%%%%%%%
C=inv([L_a L_b; U_a U_b])*[L_d; U_d]*eigenvector(:,1);
denc=C(1:n);
numc=C(n+1:end);
%%%%%%%% Verify whether numc and denc are coprime!
m=(gcd(poly2sym(numc',c),poly2sym(denc',c)));
if isa(m,'sym')
    m0=sym2poly(m);
    m1=m0(1);
    [numc1,non]=deconv(numc',m0); numc1= numc1*m1;
    [denc1,non]=deconv(denc',m0); denc1= denc1*m1;
end
C=tf(numc1,denc1,1);
%%% Compute for stability margin
B_pc=1/(sqrt(1+(value)^2));
end