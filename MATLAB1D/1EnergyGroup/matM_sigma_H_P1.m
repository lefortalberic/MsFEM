function MM_sigma = matM_sigma_H_P1(Nbpt_H,Nbpt_ref, epsilon)
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
H = 1/(Nbpt_H-1);
Nbpt_h = (Nbpt_ref-1)/(Nbpt_H-1) +1;
h_ref = 1/(Nbpt_ref-1);
% declarations
% ------------

MM_sigma = sparse(Nbpt_H,Nbpt_H); % matrice de masse

% boucle pour les matrices EF
% ------------------------
for l=2:(Nbpt_H-2)   %La première et dernière ligne doivent rester nulle

    Somme_M_sigma_1 = 0;
    Somme_M_sigma_2 = 0;

    for j=1:(2*Nbpt_h-2)
        if (j>=Nbpt_h)
            xj = 1-((j+0.5-Nbpt_h)*h_ref)/H;
        else 
            xj = (j+0.5-1)*h_ref/H;
        end
        Somme_M_sigma_1 = Somme_M_sigma_1 + Sigma((l-2)*H+(j-1+0.5)*h_ref,epsilon)*xj*xj;
    end
    
    for j=1:(Nbpt_h-1)
        xj = (j+0.5-1)*h_ref/H;
        Somme_M_sigma_2 = Somme_M_sigma_2 + Sigma((l-1)*H+(j+0.5-1)*h_ref,epsilon)*xj*(1-xj);
        %Somme_M_sigma_2 = Somme_M_sigma_2 + Sigma((l-1)*H+(j-1+0.5)*h,epsilon)*0.5*(Khi_g(Nbpt_h+j)+Khi_g(Nbpt_h-1+j))*0.5*(Khi_d(j+1)+Khi_d(j));
    end
    
    MM_sigma(l,l)= h_ref*Somme_M_sigma_1;
    MM_sigma(l,l+1)=h_ref*Somme_M_sigma_2;
    MM_sigma(l+1,l)=h_ref*Somme_M_sigma_2;

end % for l

%Dernière case des matrices
Somme_M_sigma_1 = 0;
l=Nbpt_H-1;

for j=1:(2*Nbpt_h-2)
    if (j>=Nbpt_h)
        xj = 1-((j+0.5-Nbpt_h)*h_ref)/H;
    else
        xj = (j+0.5-1)*h_ref/H;
    end
    Somme_M_sigma_1 = Somme_M_sigma_1 + Sigma((l-2)*H+(j+0.5-1)*h_ref,epsilon)*xj*xj;
end
MM_sigma(l,l)= h_ref*Somme_M_sigma_1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
