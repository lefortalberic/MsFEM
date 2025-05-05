function CC = matC(Khi_global,vecteur_J, Nbpt_H, Nbpt_h, epsilon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matM :
% calcul la matrice de masse
%
% INPUT *
% Num_coeff = 1 ou 2. Num_coeff = 1 Pour calculer la matrice de raideur
% correspondant au premier niveaux d'energie
% OUTPUT - Mel matrice de masse
%
% NOTE (1) calcul de l'intégrale par une quadrature
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = sparse(Nbpt_H,Nbpt_H);
H=1/(Nbpt_H-1);
h = 1/(Nbpt_h -1)*H;
Nbpt_ref = (Nbpt_h-1)*(Nbpt_H-1)+1;

vecteur_J_eps_g = zeros(Nbpt_h,1);
vecteur_J_eps_d = zeros(Nbpt_h,1);

%Khi_1
partie_gauche_fct_forme_gauche = Khi_global(:,1,1);
partie_droite_fct_forme_gauche = Khi_global(:,1,2);


for l=2:(Nbpt_H-2)   %La première et dernière ligne doivent rester nulle

    % Khi_(l)
    partie_gauche_fct_forme_droite = Khi_global(:,l,1);
    partie_droite_fct_forme_droite = Khi_global(:,l,2);

    for i=1:Nbpt_h
        xi_g = (l-1-1)*H+(i-1)*h ;
        xi_d = (l-1)*H+(i-1)*h ;
        yi_g = rem(xi_g/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk_g = fix(yi_g*(Nbpt_ref-1)) +1; %partie entière de y*(Nbpt-1) (indice pour les matrices)
        yi_d = rem(xi_d/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk_d = fix(yi_d*(Nbpt_ref-1)) +1;
        t_g = yi_g*(Nbpt_ref-1) - fix(yi_g*(Nbpt_ref-1)) ;
        t_d = yi_d*(Nbpt_ref-1) - fix(yi_d*(Nbpt_ref-1)) ;
        vecteur_J_eps_g(i) = vecteur_J(yk_g)*(1-t_g) +vecteur_J(yk_g+1)*t_g;
        vecteur_J_eps_d(i) = vecteur_J(yk_d)*(1-t_d) +vecteur_J(yk_d+1)*t_d;
    end

    Somme_K_1 = 0;
    Somme_K_2 = 0;
    Somme_K_3 = 0;

    for j=1:(Nbpt_h-1)
        vecteur_J_g = 0.5*(vecteur_J_eps_g(j)+vecteur_J_eps_g(j+1));
        vecteur_J_d = 0.5*(vecteur_J_eps_d(j)+vecteur_J_eps_d(j+1));

        Somme_K_1 = Somme_K_1 + vecteur_J_g*((partie_gauche_fct_forme_gauche(j+1)-partie_gauche_fct_forme_gauche(j))/h)*0.5*(partie_gauche_fct_forme_gauche(j+1)+partie_gauche_fct_forme_gauche(j)); %gradient de Khi_i au carré
        Somme_K_1 = Somme_K_1 + vecteur_J_d*((partie_droite_fct_forme_gauche(j+1)-partie_droite_fct_forme_gauche(j))/h)*0.5*(partie_droite_fct_forme_gauche(j+1)+partie_droite_fct_forme_gauche(j)); %gradient de Khi_i au carré
    end

    for j=1:(Nbpt_h-1)
        vecteur_J_d = 0.5*(vecteur_J_eps_d(j)+vecteur_J_eps_d(j+1));
        Somme_K_2 = Somme_K_2 + vecteur_J_d*0.5*(partie_droite_fct_forme_gauche(j+1)+partie_droite_fct_forme_gauche(j))*((partie_gauche_fct_forme_droite(j+1)-partie_gauche_fct_forme_droite(j))/h); %gradient de Khi_i*gradient de Khi_j
    end

    for j=1:(Nbpt_h-1)
        vecteur_J_d = 0.5*(vecteur_J_eps_d(j)+vecteur_J_eps_d(j+1));
        Somme_K_3 = Somme_K_3 + vecteur_J_d*0.5*(partie_gauche_fct_forme_droite(j+1)+partie_gauche_fct_forme_droite(j))*((partie_droite_fct_forme_gauche(j+1)-partie_droite_fct_forme_gauche(j))/h); %gradient de Khi_i*gradient de Khi_j
    end

    CC(l,l)= h*Somme_K_1;
    CC(l,l+1)=h*Somme_K_2;
    CC(l+1,l)=h*Somme_K_3;

    % fct_forme_gauche=fct_forme_droite
    partie_gauche_fct_forme_gauche = partie_gauche_fct_forme_droite ;
    partie_droite_fct_forme_gauche = partie_droite_fct_forme_droite ;

end % for l

%Dernière case des matrices
Somme_K_1 = 0;
l=Nbpt_H-1;

for i=1:Nbpt_h
    xi_g = (l-1-1)*H+(i-1)*h ;
    xi_d = (l-1)*H+(i-1)*h ;
    yi_g = rem(xi_g/epsilon,1) ; % y vaut x/epsilon modulo 1
    yk_g = fix(yi_g*(Nbpt_ref-1)) +1; %partie entière de y*(Nbpt-1) (indice pour les matrices)
    yi_d = rem(xi_d/epsilon,1) ; % y vaut x/epsilon modulo 1
    yk_d = fix(yi_d*(Nbpt_ref-1)) +1;
    t_g = yi_g*(Nbpt_ref-1) - fix(yi_g*(Nbpt_ref-1)) ;
    t_d = yi_d*(Nbpt_ref-1) - fix(yi_d*(Nbpt_ref-1)) ;
    vecteur_J_eps_g(i) = vecteur_J(yk_g)*(1-t_g) +vecteur_J(yk_g+1)*t_g;
    vecteur_J_eps_d(i) = vecteur_J(yk_d)*(1-t_d) +vecteur_J(yk_d+1)*t_d;
end

for j=1:(Nbpt_h-1)
    vecteur_J_g = 0.5*(vecteur_J_eps_g(j)+vecteur_J_eps_g(j+1));
    vecteur_J_d = 0.5*(vecteur_J_eps_d(j)+vecteur_J_eps_d(j+1));

    Somme_K_1 = Somme_K_1 + vecteur_J_g*((partie_gauche_fct_forme_gauche(j+1)-partie_gauche_fct_forme_gauche(j))/h)*0.5*(partie_gauche_fct_forme_gauche(j+1)+partie_gauche_fct_forme_gauche(j)); %gradient de Khi_i au carré
    Somme_K_1 = Somme_K_1 + vecteur_J_d*((partie_droite_fct_forme_gauche(j+1)-partie_droite_fct_forme_gauche(j))/h)*0.5*(partie_droite_fct_forme_gauche(j+1)+partie_droite_fct_forme_gauche(j)); %gradient de Khi_i au carré
end

CC(l,l)= h*Somme_K_1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
