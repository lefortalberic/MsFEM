function KK = matK(Khi_global,Nbpt_H,Nbpt_h,epsilon)
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
H=1/(Nbpt_H-1);
h = 1/(Nbpt_h -1)*H;

%Khi_1
partie_gauche_fct_forme_gauche = Khi_global(:,1,1);
partie_droite_fct_forme_gauche = Khi_global(:,1,2);


for l=2:(Nbpt_H-2)   %La première et dernière ligne doivent rester nulle

    % Khi_(l)
    partie_gauche_fct_forme_droite = Khi_global(:,l,1);
    partie_droite_fct_forme_droite = Khi_global(:,l,2);

    Somme_K_1 = 0;
    Somme_K_2 = 0;

    for j=1:(Nbpt_h-1)
        Somme_K_1 = Somme_K_1 + A((l-2)*H+(j-0.5)*h,epsilon)*((partie_gauche_fct_forme_gauche(j+1)-partie_gauche_fct_forme_gauche(j))/h)^2; %gradient de Khi_i au carré
        Somme_K_1 = Somme_K_1 + A((l-1)*H+(j-0.5)*h,epsilon)*((partie_droite_fct_forme_gauche(j+1)-partie_droite_fct_forme_gauche(j))/h)^2; %gradient de Khi_i au carré
    end

    for j=1:(Nbpt_h-1)
        Somme_K_2 = Somme_K_2 + A((l-1)*H+(j-0.5)*h,epsilon)*((partie_droite_fct_forme_gauche(j+1)-partie_droite_fct_forme_gauche(j))/h)*((partie_gauche_fct_forme_droite(j+1)-partie_gauche_fct_forme_droite(j))/h); %gradient de Khi_i*gradient de Khi_j
    end

    KK(l,l)= h*Somme_K_1;
    KK(l,l+1)=h*Somme_K_2;
    KK(l+1,l)=h*Somme_K_2;

    % fct_forme_gauche=fct_forme_droite
    partie_gauche_fct_forme_gauche = partie_gauche_fct_forme_droite ;
    partie_droite_fct_forme_gauche = partie_droite_fct_forme_droite ;

end % for l

%Dernière case des matrices
Somme_K_1 = 0;
l=Nbpt_H-1;

for j=1:(Nbpt_h-1)
    Somme_K_1 = Somme_K_1 + A((l-2)*H+(j-0.5)*h,epsilon)*((partie_gauche_fct_forme_gauche(j+1)-partie_gauche_fct_forme_gauche(j))/h)^2; %gradient de Khi_i au carré
    Somme_K_1 = Somme_K_1 + A((l-1)*H+(j-0.5)*h,epsilon)*((partie_droite_fct_forme_gauche(j+1)-partie_droite_fct_forme_gauche(j))/h)^2; %gradient de Khi_i au carré
end

KK(l,l)= h*Somme_K_1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
