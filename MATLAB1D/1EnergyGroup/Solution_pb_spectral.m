function [lambda,Psi] = Solution_pb_spectral(Nbpt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% INPUT - x : coordonnee du point ou on veut evaluer la fonction.
% OUTPUT - val: valeur de la matrice sur ce point.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------
% calcul des matrices EF
% ----------------------
MM = matM_per(Nbpt-1);
MM_sigma = matM_sigma_per(Nbpt-1); % matrice de masse avec sigma
KK = matK_per(Nbpt-1); % matrice de rigidite

AA = MM_sigma + KK ;

%r�solution probl�me valeurs propres g�n�ralis�
[Psi,lambda] = eigs(AA,MM,1,'smallestabs');
mean_Psi = mean(Psi);
Psi = Psi./mean_Psi;
Psi = abs([Psi ; Psi(1)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
