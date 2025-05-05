function MM = matM_filtre(Nbpt,h,R)
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
for l=1:(Nbpt-1)   %La première et dernière ligne doivent rester nulle

    MM(l,l)=h*(1/3)*(Phi_R((l-1-0.5)*h,R) + Phi_R((l-1+0.5)*h,R));
    MM(l,l+1)=h*(1/6)*Phi_R((l-1+0.5)*h,R);
    MM(l+1,l)=h*(1/6)*Phi_R((l-1+0.5)*h,R);
        
end % for l

MM(1,1)=h*(1/3)*Phi_R((0.5)*h,R);
MM(Nbpt,Nbpt)=h*(1/3)*Phi_R(1-(0.5)*h,R);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
