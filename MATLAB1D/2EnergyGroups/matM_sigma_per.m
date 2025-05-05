function MM_sigma = matM_sigma_per(Nbpt_per,Num_coeff)
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
h = 1/(Nbpt_per);
% declarations
% ------------

MM_sigma = sparse(Nbpt_per,Nbpt_per); % matrice de masse

% boucle pour les matrices EF
% ------------------------
for l=2:(Nbpt_per-1)   %La première et dernière ligne doivent rester nulle

    MM_sigma(l,l)=(h/3)*(Sigma((l-1-0.5)*h,1,Num_coeff)+Sigma((l-1+0.5)*h,1,Num_coeff));
    MM_sigma(l,l+1)=h*(1/6)*Sigma((l-1+0.5)*h,1,Num_coeff);
    MM_sigma(l+1,l)=h*(1/6)*Sigma((l-1+0.5)*h,1,Num_coeff);
            
end % for l

MM_sigma(1,1)=(h/3)*(Sigma(0.5*h,1,Num_coeff)+Sigma(1-0.5*h,1,Num_coeff));
MM_sigma(1,2)=h*(1/6)*Sigma(0.5*h,1,Num_coeff);
MM_sigma(2,1)=h*(1/6)*Sigma(0.5*h,1,Num_coeff);
MM_sigma(1,Nbpt_per)=h*(1/6)*Sigma(1-0.5*h,1,Num_coeff);
MM_sigma(Nbpt_per,1)=h*(1/6)*Sigma(1-0.5*h,1,Num_coeff);

MM_sigma(Nbpt_per,Nbpt_per)=(h/3)*(Sigma((Nbpt_per-1-0.5)*h,1,Num_coeff)+Sigma((Nbpt_per-1+0.5)*h,1,Num_coeff));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
