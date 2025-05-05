function [lambda,Psi_1,Psi_2] = Solution_pb_spectral(Nbpt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% INPUT - x : coordonnee du point ou on veut evaluer la fonction.
% OUTPUT - val: valeur de la matrice sur ce point.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------
% calcul des matrices EF
% ----------------------
MM = matM_per(Nbpt-1);
Zeros = sparse(Nbpt-1,Nbpt-1); % matrice de masse

MM_sigma11 = matM_sigma_per(Nbpt-1,1); % matrice de masse avec sigma
MM_sigma12 = matM_sigma_per(Nbpt-1,2); % matrice de masse avec sigma
MM_sigma21 = matM_sigma_per(Nbpt-1,3); % matrice de masse avec sigma
MM_sigma22 = matM_sigma_per(Nbpt-1,4); % matrice de masse avec sigma
KK1 = matK_per(Nbpt-1,1); % matrice de rigidite
KK2 = matK_per(Nbpt-1,2); % matrice de rigidite

KK = [KK1 , Zeros ; Zeros ,KK2 ];
MM_sigma = [MM_sigma11 , MM_sigma12 ; MM_sigma21 , MM_sigma22];
MM = [MM , Zeros ; Zeros ,MM ];

AA = MM_sigma + KK ;

%résolution problème valeurs propres généralisé
[Psi,lambda] = eigs(AA,MM,1,'smallestabs');
Psi_1 = Psi(1:Nbpt-1);
Psi_2 = Psi(Nbpt:2*(Nbpt-1));
mean_Psi = sqrt(mean(Psi_1.*Psi_1 + Psi_2.*Psi_2));
mean_Psi_1_carre = sqrt(mean(Psi_1.*Psi_1));
Psi_1 = Psi_1./mean_Psi;
mean_Psi_2_carre = sqrt(mean(Psi_2.*Psi_2));
Psi_2 = Psi_2./mean_Psi;

Psi_1 = abs([Psi_1 ; Psi_1(1)]);
Psi_2 = abs([Psi_2 ; Psi_2(1)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
