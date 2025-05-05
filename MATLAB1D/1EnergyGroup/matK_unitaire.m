function KK = matK_unitaire(Nbpt,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT - Mel matrice de masse 
%
% NOTE (1) calcul de l'intégrale par une quadrature 
%      
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% declarations
% ------------

KK = sparse(Nbpt,Nbpt); % matrice de masse

% boucle pour les matrices EF
% ------------------------
for l=1:(Nbpt-1)   %La première et dernière ligne doivent rester nulle

    KK(l,l)=(2/h); %approximation de l'intégrale de a(.)
    KK(l,l+1)=(-1/h);
    KK(l+1,l)=(-1/h);
        
end % for l

KK(1,1)=(1/h);
KK(Nbpt,Nbpt)=(1/h);
end
