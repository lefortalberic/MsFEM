function MM = matM_ref_Psi(Nbpt, Psi , Psi_adjoint, h)
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

MM = sparse(Nbpt,Nbpt); % matrice de masse

% boucle pour les matrices EF
% ------------------------
for l=2:(Nbpt-2)   %La première et dernière ligne doivent rester nulle
    psi_avant = 0.5*(Psi(l)+Psi(l-1));
    psi_apres = 0.5*(Psi(l)+Psi(l+1));

    psi_avant_adjoint = 0.5*(Psi_adjoint(l)+Psi_adjoint(l-1));
    psi_apres_adjoint = 0.5*(Psi_adjoint(l)+Psi_adjoint(l+1));

    MM(l,l)=(h/3)*(psi_avant*psi_avant_adjoint+psi_apres*psi_apres_adjoint);
    MM(l,l+1)=h*(1/6)*psi_apres*psi_apres_adjoint;
    MM(l+1,l)=h*(1/6)*psi_apres*psi_apres_adjoint;
            
end % for l
psi_avant = 0.5*(Psi(Nbpt-1)+Psi(Nbpt-1-1));
psi_apres = 0.5*(Psi(Nbpt-1)+Psi(Nbpt-1+1));

psi_avant_adjoint = 0.5*(Psi_adjoint(Nbpt-1)+Psi_adjoint(Nbpt-1-1));
psi_apres_adjoint = 0.5*(Psi_adjoint(Nbpt-1)+Psi_adjoint(Nbpt-1+1));

MM(Nbpt-1,Nbpt-1)=(h/3)*(psi_avant*psi_avant_adjoint+psi_apres*psi_apres_adjoint);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
