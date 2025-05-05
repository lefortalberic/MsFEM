function MM = matM_per(Nbpt_per)
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

MM = sparse(Nbpt_per,Nbpt_per); % matrice de masse

% boucle pour les matrices EF
% ------------------------
for l=1:(Nbpt_per-1)   

    MM(l,l)=h*(2/3);
    MM(l,l+1)=h*(1/6);
    MM(l+1,l)=h*(1/6);
        
end % for l
MM(1,1)=h*(2/3);
MM(1,2)=h*(1/6);
MM(2,1)=h*(1/6);
MM(1,Nbpt_per)=h*(1/6);
MM(Nbpt_per,1)=h*(1/6);
MM(Nbpt_per,Nbpt_per)=h*(2/3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
