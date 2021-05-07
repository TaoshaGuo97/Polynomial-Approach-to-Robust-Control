function [Dz,Dzi] = Spectral_Factorization(az,bz)
%  Spectral_Factorization(az,bz): The Spectral Factorization of a stable discrete system bz/az.
%  Spectral_Factorization(az,bz) is to find the polynomial Dz and Dzi such
%  that az(z)az(z^{-1})+bz(z)bz(z^{-1})=Dz(z)Dzi(z^{-1}).
% 
%    Input "az" and "bz" are the coefficients of their corresponding polynomials (created with vetors).
%
%    [Dz,Dzi,erro] = Spectral_Factorization(az,bz) returns the Dz and Dzi.
%   
%    An Example:
%
%     az=[1 4 4 ]; 
%     bz=[0 2 1 ];
%     [Dz,Dzi] = Spectral_Factorization(az,bz);
%     Dz =
%     4.5765    4.0363    0.8740
%     Dzi =
%     0.8740    4.0363    4.5765
syms c
az_i=fliplr(az);
bz_i=fliplr(bz);
n=length(az);
% Spactral Factorization
DzDzi=conv(az,az_i)+conv(bz,bz_i);
r=roots(DzDzi);
%Verify the acuurancy of poles, not bad
%mm=(c+root(1))*(c+root(2))*(c+root(3))*(c+root(4));
%mm=sym2poly(mm)
% mm-DzDzi
% ans =s
%
%   1.0e-13 *
%
%         0   -0.1066   -0.2842   -0.2132   -0.0444
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the stable poles and unstable poles
ppoles = zeros(1,n-1); % Positive poles
npoles = zeros(1,n-1); % Negetive poles
j=1;
for i=1:2*(n-1)
    if abs(r(i))<1 
    ppoles(j)=r(i);
    npoles(j)= 1/r(i);
        j=j+1; 
    end
end
%%%%%%%%%%%%%%% polesplacement to DzDz
Dz=1; 
for i =1:n-1
    Dz=Dz*(c-ppoles(i));
end
%%%%%%%%%% sscal it;
Dz=sym2poly(Dz);
DzDzif=conv(Dz,fliplr(Dz));
sca=DzDzi(1)/DzDzif(1);
Dz=Dz*(sqrt(sca));  % d(z)
Dzi=fliplr(Dz);     % d(z^-1)
%erro= conv(Dz,Dzi)-DzDzi;
end


