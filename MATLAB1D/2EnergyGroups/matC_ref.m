function CC = matC_ref(Nbpt, vecteur_J)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matM :
% calcul la matrice de masse 
%          
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

CC = sparse(Nbpt,Nbpt); % matrice de masse

% boucle pour les matrices EF
% ------------------------
for l=2:(Nbpt-2)   %La première et dernière ligne doivent rester nulle
    vecteur_J_avant = 0.5*(vecteur_J(l)+vecteur_J(l-1));
    vecteur_J_apres = 0.5*(vecteur_J(l)+vecteur_J(l+1));

    CC(l,l)=(1/2)*(vecteur_J_avant-vecteur_J_apres);
    CC(l,l+1)=(1/2)*(vecteur_J_apres);
    CC(l+1,l)=(-1/2)*(vecteur_J_apres);
            
end % for l
vecteur_J_avant = 0.5*(vecteur_J(Nbpt-1)+vecteur_J(Nbpt-1-1));
vecteur_J_apres = 0.5*(vecteur_J(Nbpt-1)+vecteur_J(Nbpt-1+1));

CC(Nbpt-1,Nbpt-1)=(1/2)*(vecteur_J_avant-vecteur_J_apres);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
