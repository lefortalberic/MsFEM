function KK = matK_homo(Nbpt,A_etoile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT - Mel matrice de masse 
%
% NOTE (1) calcul de l'int�grale par une quadrature 
%      
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = 1/(Nbpt-1);
% declarations
% ------------

KK = sparse(Nbpt,Nbpt); % matrice de masse

% boucle pour les matrices EF
% ------------------------
for l=2:(Nbpt-2)   %La premi�re et derni�re ligne doivent rester nulle

    KK(l,l)=(2/h)*A_etoile; %approximation de l'int�grale de a(.)
    KK(l,l+1)=(-1/h)*A_etoile;
    KK(l+1,l)=(-1/h)*A_etoile;
        
end % for l

KK(Nbpt-1,Nbpt-1)=(2/h)*A_etoile;
end
