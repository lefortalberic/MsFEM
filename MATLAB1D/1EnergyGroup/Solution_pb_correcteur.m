function W = Solution_pb_correcteur(Nbpt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% INPUT - x : coordonnee du point ou on veut evaluer la fonction.
% OUTPUT - val: valeur de la matrice sur ce point.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = 1/(Nbpt -1);
% ----------------------
[~,Psi] = Solution_pb_spectral(Nbpt);

KK = matK_per_Psi(Nbpt-1,Psi); % matrice de rigidite avec A remplacé par Psi**2*A


H = zeros(Nbpt-1,1);
for i=2:(Nbpt-1)
    H(i)= 0.5*(Psi(i)+Psi(i-1))*0.5*(Psi(i)+Psi(i-1))*A((i-2+0.5)*h,1) - 0.5*(Psi(i+1)+Psi(i))*0.5*(Psi(i+1)+Psi(i))*A((i-1+0.5)*h,1);
end
H(1)= -0.5*(Psi(2)+Psi(1))*0.5*(Psi(2)+Psi(1))*A(0.5*h,1) + 0.5*(Psi(Nbpt)+Psi(Nbpt-1))*0.5*(Psi(Nbpt)+Psi(Nbpt-1))*A((Nbpt-1-0.5)*h,1);

AA = KK ;

FF = -H;
warning('off', 'all');
W = AA\FF;
warning('on', 'all');
W = W - mean(W);
W = [W ; W(1)];

%  W_prime = [W(2:end)-W(1:end-1)]./h;
%  plot(W_prime)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
