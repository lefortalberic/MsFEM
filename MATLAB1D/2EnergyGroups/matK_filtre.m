function KK = matK_filtre(Nbpt,epsilon, x_gauche, h, R, Num_coeff)
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

KK = sparse(Nbpt,Nbpt); % matrice de masse

% boucle pour les matrices EF
% ------------------------
for l=1:(Nbpt-1)   %La première et dernière ligne doivent rester nulle

    KK(l,l)=(1/h)*(A(x_gauche+(l-1-0.5)*h,epsilon, Num_coeff)*Phi_R((l-1-0.5)*h,R)+A(x_gauche+(l-1+0.5)*h,epsilon, Num_coeff)*Phi_R((l-1+0.5)*h,R)); %approximation de l'intégrale de a(.)
    KK(l,l+1)=(-1/h)*A(x_gauche+(l-1+0.5)*h,epsilon, Num_coeff)*Phi_R((l-1+0.5)*h,R);
    KK(l+1,l)=(-1/h)*A(x_gauche+(l-1+0.5)*h,epsilon, Num_coeff)*Phi_R((l-1+0.5)*h,R);
        
end % for l
KK(1,1)=(1/h)*(A(x_gauche+0.5*h,epsilon, Num_coeff))*Phi_R((0.5)*h,R);
KK(Nbpt,Nbpt)=(1/h)*(A(x_gauche+(Nbpt-1-0.5)*h,epsilon, Num_coeff))*Phi_R(1-(0.5)*h,R);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
