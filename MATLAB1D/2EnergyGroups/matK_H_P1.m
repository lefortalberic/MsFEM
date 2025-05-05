function KK = matK_H_P1(Nbpt_H,Nbpt_ref, epsilon)
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
KK = sparse(Nbpt_H,Nbpt_H);

H = 1/(Nbpt_H-1);
Nbpt_h = (Nbpt_ref-1)/(Nbpt_H-1) +1;
h_ref = 1/(Nbpt_ref-1);

for l=2:(Nbpt_H-2)   %La première et dernière ligne doivent rester nulle
      
    Somme_K_1 = 0;
    Somme_K_2 = 0;

    for j=1:(2*Nbpt_h-2)
        Somme_K_1 = Somme_K_1 + A((l-2)*H+(j-0.5)*h_ref,epsilon)/(H*H); %gradient de Khi_i au carré
    end

    for j=1:(Nbpt_h-1)
        Somme_K_2 = Somme_K_2 + A((l-1)*H+(j-0.5)*h_ref,epsilon)*(-1)/(H*H); %gradient de Khi_i*gradient de Khi_j
    end
        
    KK(l,l)= h_ref*Somme_K_1;
    KK(l,l+1)=h_ref*Somme_K_2;
    KK(l+1,l)=h_ref*Somme_K_2;
    
end % for l

%Dernière case des matrices
Somme_K_1 = 0;
l=Nbpt_H-1;

for j=1:(2*Nbpt_h-2)
    Somme_K_1 = Somme_K_1 + A((l-2)*H+(j-0.5)*h_ref,epsilon)/(H*H); %gradient de Khi_i au carré
end

KK(l,l)= h_ref*Somme_K_1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
