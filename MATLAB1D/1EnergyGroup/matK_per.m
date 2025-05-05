function KK = matK_per(Nbpt_per)
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
%h = 1/(Nbpt-1);
h = 1/(Nbpt_per);   %attention ici h est modifié car on précalcul une sous matrice
% declarations
% ------------

KK = sparse(Nbpt_per,Nbpt_per); % matrice de masse

% boucle pour les matrices EF
% ------------------------
for l=2:(Nbpt_per-1)   %La première et dernière ligne doivent rester nulle

    KK(l,l)=(1/h)*(A((l-1-0.5)*h,1)+A((l-1+0.5)*h,1)); %approximation de l'intégrale de a(.)
    KK(l,l+1)=(-1/h)*A((l-1+0.5)*h,1);
    KK(l+1,l)=(-1/h)*A((l-1+0.5)*h,1);
        
end % for l

KK(1,1)=(1/h)*(A(0.5*h,1)+A(1-0.5*h,1));
KK(1,2)=(-1/h)*A(0.5*h,1);
KK(2,1)=(-1/h)*A(0.5*h,1);
KK(1,Nbpt_per)=(-1/h)*A(1-0.5*h,1);
KK(Nbpt_per,1)=(-1/h)*A(1-0.5*h,1);

KK(Nbpt_per,Nbpt_per)=(1/h)*(A((Nbpt_per-1-0.5)*h,1)+A((Nbpt_per-1+0.5)*h,1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
