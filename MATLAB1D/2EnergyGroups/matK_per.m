function KK = matK_per(Nbpt_per,Num_coeff)
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

    KK(l,l)=(1/h)*(A((l-1-0.5)*h,1,Num_coeff)+A((l-1+0.5)*h,1,Num_coeff)); %approximation de l'intégrale de a(.)
    KK(l,l+1)=(-1/h)*A((l-1+0.5)*h,1,Num_coeff);
    KK(l+1,l)=(-1/h)*A((l-1+0.5)*h,1,Num_coeff);
        
end % for l

KK(1,1)=(1/h)*(A(0.5*h,1,Num_coeff)+A(1-0.5*h,1,Num_coeff));
KK(1,2)=(-1/h)*A(0.5*h,1,Num_coeff);
KK(2,1)=(-1/h)*A(0.5*h,1,Num_coeff);
KK(1,Nbpt_per)=(-1/h)*A(1-0.5*h,1,Num_coeff);
KK(Nbpt_per,1)=(-1/h)*A(1-0.5*h,1,Num_coeff);

KK(Nbpt_per,Nbpt_per)=(1/h)*(A((Nbpt_per-1-0.5)*h,1,Num_coeff)+A((Nbpt_per-1+0.5)*h,1,Num_coeff));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
