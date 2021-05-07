function [H]=Hankel_Stan(P)
%  Hankel_Stan(P) is to compute the Hankel matrix H, of a stable discrete systme P=b(z)/a(z)
%  over its standard basis function [z^{n} z^{n-1} ... z 1]/a(z);
% 
%    Input "P" is the given system (created with either TF, ZPK, SS, or FRD).
%
%    [H]=Hankel_Stan(P) returns the corresponding Hankel matrix H
%
%     An Example:
%     [H]=Hankel_Stan(P)
%     P = tf([0 2 1],[1 4 4],1);
%    H =
%
%   10.4721   -1.0000
%    8.2361   -2.0000
%
syms c;
[tempB,tempA]=tfdata(P);
az=tempA{1,1};
az_i= fliplr(az);
bz=tempB{1,1};
bz_i= fliplr(bz);
[dz,na]=Spectral_Factorization(az,bz);
n=length(az)-1;
m=length(bz)-1;
% Generates the matrix Ad.
Ad=zeros(n);
for k=1:n
        Ad(k,1)=-dz(1,k+1)/dz(1,1);
end
if n>1
   for k=1:n-1
       Ad(k,k+1)=1;
   end
end
% Compute the Hankel matrix H.
% Generates the matrix J.
J=fliplr(eye(n));
a1=polyvalm(az_i,Ad); %
b1=polyvalm(bz_i,Ad); %
a2=polyvalm(az,Ad); %a(Ad)
b2=polyvalm(bz,Ad); %b(Ad)
H=inv([a1;b1]'*[a1;b1])*[a1;b1]'*[b2;-a2]*J;
end
