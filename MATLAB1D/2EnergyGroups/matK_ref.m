function KK = matK_ref(Nbpt,epsilon, Num_coeff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matM :
% calcul la matrice de masse 
%          
% INPUT * 
%
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

    KK(l,l)=(1/h)*(A((l-1-0.5)*h,epsilon, Num_coeff)+A((l-1+0.5)*h,epsilon, Num_coeff)); %approximation de l'int�grale de a(.)
    KK(l,l+1)=(-1/h)*A((l-1+0.5)*h,epsilon, Num_coeff);
    KK(l+1,l)=(-1/h)*A((l-1+0.5)*h,epsilon, Num_coeff);
        
end % for l

KK(Nbpt-1,Nbpt-1)=(1/h)*(A((Nbpt-2-0.5)*h,epsilon, Num_coeff)+A((Nbpt-2+0.5)*h,epsilon, Num_coeff));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
