function MM_sigma = matM_sigma_ref(Nbpt,epsilon)
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
h = 1/(Nbpt-1);
% declarations
% ------------

MM_sigma = sparse(Nbpt,Nbpt); % matrice de masse

% boucle pour les matrices EF
% ------------------------
for l=2:(Nbpt-2)   %La première et dernière ligne doivent rester nulle

    MM_sigma(l,l)=(h/3)*(Sigma((l-1-0.5)*h,epsilon)+Sigma((l-1+0.5)*h,epsilon));
    MM_sigma(l,l+1)=h*(1/6)*Sigma((l-1+0.5)*h,epsilon);
    MM_sigma(l+1,l)=h*(1/6)*Sigma((l-1+0.5)*h,epsilon);
            
end % for l

MM_sigma(Nbpt-1,Nbpt-1)=(h/3)*(Sigma((Nbpt-2-0.5)*h,epsilon)+Sigma((Nbpt-2+0.5)*h,epsilon));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
