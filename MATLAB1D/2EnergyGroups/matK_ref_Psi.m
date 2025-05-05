function KK = matK_ref_Psi(Nbpt,epsilon, Psi , Psi_adjoint, Num_coeff, x_gauche, h)
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

KK = sparse(Nbpt,Nbpt); % matrice de masse

% boucle pour les matrices EF
% ------------------------
for l=2:(Nbpt-2)   %La première et dernière ligne doivent rester nulle
    
    psi_avant = 0.5*(Psi(l)+Psi(l-1));
    psi_apres = 0.5*(Psi(l)+Psi(l+1));

    psi_avant_adjoint = 0.5*(Psi_adjoint(l)+Psi_adjoint(l-1));
    psi_apres_adjoint = 0.5*(Psi_adjoint(l)+Psi_adjoint(l+1));

    KK(l,l)=(1/h)*(psi_avant*psi_avant_adjoint*A(x_gauche+(l-1-0.5)*h,epsilon,Num_coeff) ...
        +psi_apres*psi_apres_adjoint*A(x_gauche+(l-1+0.5)*h,epsilon,Num_coeff)); %approximation de l'intégrale de a(.)
    KK(l,l+1)=(-1/h)*psi_apres*psi_apres_adjoint*A(x_gauche+(l-1+0.5)*h,epsilon,Num_coeff);
    KK(l+1,l)=(-1/h)*psi_apres*psi_apres_adjoint*A(x_gauche+(l-1+0.5)*h,epsilon,Num_coeff);
        
end % for l
psi_avant = 0.5*(Psi(Nbpt-1)+Psi(Nbpt-1-1));
psi_apres = 0.5*(Psi(Nbpt-1)+Psi(Nbpt-1+1));

psi_avant_adjoint = 0.5*(Psi_adjoint(Nbpt-1)+Psi_adjoint(Nbpt-1-1));
psi_apres_adjoint = 0.5*(Psi_adjoint(Nbpt-1)+Psi_adjoint(Nbpt-1+1));

KK(Nbpt-1,Nbpt-1)=(1/h)*(psi_avant*psi_avant_adjoint*A(x_gauche+(Nbpt-2-0.5)*h,epsilon,Num_coeff)+psi_apres*psi_apres_adjoint*A(x_gauche+(Nbpt-2+0.5)*h,epsilon,Num_coeff));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
