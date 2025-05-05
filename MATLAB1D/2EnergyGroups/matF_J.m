function LL = matF_J(Nbpt, vecteur_J,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matM :
% calcul la matrice de masse 
%          %x_gauche le point le plus à gauche à partir duquel on calcule
%          la matrice
%           % h le pas du maillage
% INPUT * 
%
% OUTPUT - Mel matrice de masse 
%
% NOTE (1) calcul de l'intégrale par une quadrature 
%      
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% declarations
% ------------

LL = zeros(Nbpt,1); % matrice de masse

% boucle pour les matrices EF
% ------------------------
for l=2:(Nbpt-1)   %La première et dernière ligne doivent rester nulle
    
    vecteur_J_avant = 0.5*(vecteur_J(l)+vecteur_J(l-1));
    vecteur_J_apres = 0.5*(vecteur_J(l)+vecteur_J(l+1));

    LL(l)=(h/2)*(vecteur_J_avant+ vecteur_J_apres);

end % for l
end
