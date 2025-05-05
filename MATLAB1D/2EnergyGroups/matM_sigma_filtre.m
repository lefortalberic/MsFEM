function MM_sigma = matM_sigma_filtre(Nbpt,epsilon, x_gauche, h, R, Num_coeff)
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

MM_sigma = sparse(Nbpt,Nbpt); % matrice de masse

% boucle pour les matrices EF
% ------------------------
for l=1:(Nbpt-1)   %La première et dernière ligne doivent rester nulle

    MM_sigma(l,l)=(h/3)*(Sigma(x_gauche+(l-1-0.5)*h,epsilon, Num_coeff)*Phi_R((l-1-0.5)*h,R)+Sigma(x_gauche+(l-1+0.5)*h,epsilon, Num_coeff)*Phi_R((l-1+0.5)*h,R));
    MM_sigma(l,l+1)=h*(1/6)*Sigma(x_gauche+(l-1+0.5)*h,epsilon, Num_coeff)*Phi_R((l-1+0.5)*h,R);
    MM_sigma(l+1,l)=h*(1/6)*Sigma(x_gauche+(l-1+0.5)*h,epsilon, Num_coeff)*Phi_R((l-1+0.5)*h,R);
            
end % for l

MM_sigma(1,1)= (h/3)*(Sigma(x_gauche+(0.5)*h,epsilon, Num_coeff))*Phi_R((0.5)*h,R);
MM_sigma(Nbpt,Nbpt)=(h/3)*(Sigma(x_gauche+(Nbpt-1-0.5)*h,epsilon, Num_coeff))*Phi_R(1-(0.5)*h,R);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
