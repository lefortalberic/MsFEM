function MM = matM_vectoriel(Khi_global_1, Khi_global_2,Nbpt_H,Nbpt_h)
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
MM = sparse(Nbpt_H,Nbpt_H);
H=1/(Nbpt_H-1);
h = 1/(Nbpt_h -1)*H;

%Khi_1
partie_gauche_fct_forme_gauche_1 = Khi_global_1(:,1,1);
partie_droite_fct_forme_gauche_1 = Khi_global_1(:,1,2);
partie_gauche_fct_forme_gauche_2 = Khi_global_2(:,1,1);
partie_droite_fct_forme_gauche_2 = Khi_global_2(:,1,2);

for l=2:(Nbpt_H-2)   %La première et dernière ligne doivent rester nulle

    % Khi_(l)
    partie_gauche_fct_forme_droite_1 = Khi_global_1(:,l,1);
    partie_droite_fct_forme_droite_1 = Khi_global_1(:,l,2);
    partie_gauche_fct_forme_droite_2 = Khi_global_2(:,l,1);
    partie_droite_fct_forme_droite_2 = Khi_global_2(:,l,2);

    Somme_M_1 = 0;
    Somme_M_2 = 0;
    Somme_M_3 = 0;

    for j=1:(Nbpt_h-1)
        demi_somme_gauche_1 = 0.5*(partie_gauche_fct_forme_gauche_1(j+1) + partie_gauche_fct_forme_gauche_1(j));
        demi_somme_droite_1 = 0.5*(partie_droite_fct_forme_gauche_1(j+1) + partie_droite_fct_forme_gauche_1(j));
        demi_somme_gauche_2 = 0.5*(partie_gauche_fct_forme_gauche_2(j+1) + partie_gauche_fct_forme_gauche_2(j));
        demi_somme_droite_2 = 0.5*(partie_droite_fct_forme_gauche_2(j+1) + partie_droite_fct_forme_gauche_2(j));
        Somme_M_1 = Somme_M_1 + demi_somme_gauche_1*demi_somme_gauche_2 ...
            + demi_somme_droite_1*demi_somme_droite_2;
    end

    for j=1:(Nbpt_h-1)
        demi_somme = 0.5*(partie_droite_fct_forme_gauche_1(j+1)*partie_gauche_fct_forme_droite_2(j+1) ...
            + partie_droite_fct_forme_gauche_1(j)*partie_gauche_fct_forme_droite_2(j));
        Somme_M_2 = Somme_M_2 + demi_somme;
    end

    for j=1:(Nbpt_h-1)
        demi_somme = 0.5*(partie_droite_fct_forme_gauche_2(j+1)*partie_gauche_fct_forme_droite_1(j+1) ...
            + partie_droite_fct_forme_gauche_2(j)*partie_gauche_fct_forme_droite_1(j));
        Somme_M_3 = Somme_M_3 + demi_somme;
    end

    MM(l,l)= h*Somme_M_1;
    MM(l,l+1)=h*Somme_M_2;
    MM(l+1,l)=h*Somme_M_3;

    % fct_forme_gauche=fct_forme_droite
    partie_gauche_fct_forme_gauche_1 = partie_gauche_fct_forme_droite_1 ;
    partie_droite_fct_forme_gauche_1 = partie_droite_fct_forme_droite_1 ;
    partie_gauche_fct_forme_gauche_2 = partie_gauche_fct_forme_droite_2 ;
    partie_droite_fct_forme_gauche_2 = partie_droite_fct_forme_droite_2 ;

end % for l

%Dernière case des matrices
Somme_M_1 = 0;
l=Nbpt_H-1;

for j=1:(Nbpt_h-1)
    demi_somme_gauche_1 = 0.5*(partie_gauche_fct_forme_gauche_1(j+1) + partie_gauche_fct_forme_gauche_1(j));
    demi_somme_droite_1 = 0.5*(partie_droite_fct_forme_gauche_1(j+1) + partie_droite_fct_forme_gauche_1(j));
    demi_somme_gauche_2 = 0.5*(partie_gauche_fct_forme_gauche_2(j+1) + partie_gauche_fct_forme_gauche_2(j));
    demi_somme_droite_2 = 0.5*(partie_droite_fct_forme_gauche_2(j+1) + partie_droite_fct_forme_gauche_2(j));
    Somme_M_1 = Somme_M_1 + demi_somme_gauche_1*demi_somme_gauche_2 ...
        + demi_somme_droite_1*demi_somme_droite_2;
end

MM(l,l)= h*Somme_M_1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
