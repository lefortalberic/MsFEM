function MM_sigma = matM_sigma(Khi_global,Nbpt_H,Nbpt_h,epsilon)
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
MM_sigma = sparse(Nbpt_H,Nbpt_H);
H=1/(Nbpt_H-1);
h = 1/(Nbpt_h -1)*H;

%Khi_1
partie_gauche_fct_forme_gauche = Khi_global(:,1,1);
partie_droite_fct_forme_gauche = Khi_global(:,1,2);


for l=2:(Nbpt_H-2)   %La première et dernière ligne doivent rester nulle

    % Khi_(l)
    partie_gauche_fct_forme_droite = Khi_global(:,l,1);
    partie_droite_fct_forme_droite = Khi_global(:,l,2);

    Somme_M_sigma_1 = 0;
    Somme_M_sigma_2 = 0;

    for j=1:(Nbpt_h-1)
        demi_somme_gauche = 0.5*(partie_gauche_fct_forme_gauche(j+1) + partie_gauche_fct_forme_gauche(j));
        demi_somme_droite = 0.5*(partie_droite_fct_forme_gauche(j+1) + partie_droite_fct_forme_gauche(j));
        Somme_M_sigma_1 = Somme_M_sigma_1 + Sigma((l-2)*H+(j+0.5-1)*h,epsilon)*demi_somme_gauche^2 + Sigma((l-1)*H+(j+0.5-1)*h,epsilon)*demi_somme_droite^2;
    end

    for j=1:(Nbpt_h-1)
        demi_somme = Sigma((l-1)*H+(j-1+0.5)*h,epsilon)*0.5*(partie_droite_fct_forme_gauche(j+1)*partie_gauche_fct_forme_droite(j+1) + partie_droite_fct_forme_gauche(j)*partie_gauche_fct_forme_droite(j));
        Somme_M_sigma_2 = Somme_M_sigma_2 + demi_somme;
    end

    MM_sigma(l,l)= h*Somme_M_sigma_1;
    MM_sigma(l,l+1)=h*Somme_M_sigma_2;
    MM_sigma(l+1,l)=h*Somme_M_sigma_2;

    % fct_forme_gauche=fct_forme_droite
    partie_gauche_fct_forme_gauche = partie_gauche_fct_forme_droite ;
    partie_droite_fct_forme_gauche = partie_droite_fct_forme_droite ;


end % for l

%Dernière case des matrices
Somme_M_sigma_1 = 0;
l=Nbpt_H-1;

for j=1:(Nbpt_h-1)
    demi_somme_gauche = 0.5*(partie_gauche_fct_forme_gauche(j+1) + partie_gauche_fct_forme_gauche(j));
    demi_somme_droite = 0.5*(partie_droite_fct_forme_gauche(j+1) + partie_droite_fct_forme_gauche(j));
    Somme_M_sigma_1 = Somme_M_sigma_1 + Sigma((l-2)*H+(j+0.5-1)*h,epsilon)*demi_somme_gauche^2 + Sigma((l-1)*H+(j+0.5-1)*h,epsilon)*demi_somme_droite^2;
end

MM_sigma(l,l)= h*Somme_M_sigma_1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
