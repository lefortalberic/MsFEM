function MM_sigma = matM_Q_tilde(Khi_global_1, Khi_global_2, Q_tilde, Nbpt_H,Nbpt_h,epsilon)
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
Nbpt_ref = (Nbpt_h-1)*(Nbpt_H-1)+1;

Q_tilde_eps_g = zeros(Nbpt_h,1);
Q_tilde_eps_d = zeros(Nbpt_h,1);

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


    for i=1:Nbpt_h
        xi_g = (l-1-1)*H+(i-1)*h ;
        xi_d = (l-1)*H+(i-1)*h ;
        yi_g = rem(xi_g/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk_g = fix(yi_g*(Nbpt_ref-1)) +1; %partie entière de y*(Nbpt-1) (indice pour les matrices)
        yi_d = rem(xi_d/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk_d = fix(yi_d*(Nbpt_ref-1)) +1;
        t_g = yi_g*(Nbpt_ref-1) - fix(yi_g*(Nbpt_ref-1)) ;
        t_d = yi_d*(Nbpt_ref-1) - fix(yi_d*(Nbpt_ref-1)) ;
        Q_tilde_eps_g(i) = Q_tilde(yk_g)*(1-t_g) +Q_tilde(yk_g+1)*t_g;
        Q_tilde_eps_d(i) = Q_tilde(yk_d)*(1-t_d) +Q_tilde(yk_d+1)*t_d;
    end

    Somme_M_sigma_1 = 0;
    Somme_M_sigma_2 = 0;
    Somme_M_sigma_3 = 0;

    for j=1:(Nbpt_h-1)
        Q_tilde_g = 0.5*(Q_tilde_eps_g(j)+Q_tilde_eps_g(j+1));
        Q_tilde_d = 0.5*(Q_tilde_eps_d(j)+Q_tilde_eps_d(j+1));
        demi_somme_gauche_1 = 0.5*(partie_gauche_fct_forme_gauche_1(j+1) + partie_gauche_fct_forme_gauche_1(j));
        demi_somme_droite_1 = 0.5*(partie_droite_fct_forme_gauche_1(j+1) + partie_droite_fct_forme_gauche_1(j));
        demi_somme_gauche_2 = 0.5*(partie_gauche_fct_forme_gauche_2(j+1) + partie_gauche_fct_forme_gauche_2(j));
        demi_somme_droite_2 = 0.5*(partie_droite_fct_forme_gauche_2(j+1) + partie_droite_fct_forme_gauche_2(j));
        Somme_M_sigma_1 = Somme_M_sigma_1 + Q_tilde_g*demi_somme_gauche_1*demi_somme_gauche_2 + Q_tilde_d*demi_somme_droite_1*demi_somme_droite_2;
    end

    for j=1:(Nbpt_h-1)
        Q_tilde_d = 0.5*(Q_tilde_eps_d(j)+Q_tilde_eps_d(j+1));
        demi_somme = Q_tilde_d*0.5*(partie_droite_fct_forme_gauche_1(j+1)*partie_gauche_fct_forme_droite_2(j+1) ...
            + partie_droite_fct_forme_gauche_1(j)*partie_gauche_fct_forme_droite_2(j));
        Somme_M_sigma_2 = Somme_M_sigma_2 + demi_somme;
    end

    for j=1:(Nbpt_h-1)
        Q_tilde_d = 0.5*(Q_tilde_eps_d(j)+Q_tilde_eps_d(j+1));
        demi_somme = Q_tilde_d*0.5*(partie_droite_fct_forme_gauche_2(j+1)*partie_gauche_fct_forme_droite_1(j+1) ...
            + partie_droite_fct_forme_gauche_2(j)*partie_gauche_fct_forme_droite_1(j));
        Somme_M_sigma_3 = Somme_M_sigma_3 + demi_somme;
    end

    MM_sigma(l,l)= h*Somme_M_sigma_1;
    MM_sigma(l,l+1)=h*Somme_M_sigma_2;
    MM_sigma(l+1,l)=h*Somme_M_sigma_3;

    % fct_forme_gauche=fct_forme_droite
    partie_gauche_fct_forme_gauche_1 = partie_gauche_fct_forme_droite_1 ;
    partie_droite_fct_forme_gauche_1 = partie_droite_fct_forme_droite_1 ;
    partie_gauche_fct_forme_gauche_2 = partie_gauche_fct_forme_droite_2 ;
    partie_droite_fct_forme_gauche_2 = partie_droite_fct_forme_droite_2 ;


end % for l

%Dernière case des matrices
Somme_M_sigma_1 = 0;
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
    Q_tilde_eps_g(i) = Q_tilde(yk_g)*(1-t_g) +Q_tilde(yk_g+1)*t_g;
    Q_tilde_eps_d(i) = Q_tilde(yk_d)*(1-t_d) +Q_tilde(yk_d+1)*t_d;
end

for j=1:(Nbpt_h-1)
    Q_tilde_g = 0.5*(Q_tilde_eps_g(j)+Q_tilde_eps_g(j+1));
    Q_tilde_d = 0.5*(Q_tilde_eps_d(j)+Q_tilde_eps_d(j+1));
    demi_somme_gauche_1 = 0.5*(partie_gauche_fct_forme_gauche_1(j+1) + partie_gauche_fct_forme_gauche_1(j));
    demi_somme_droite_1 = 0.5*(partie_droite_fct_forme_gauche_1(j+1) + partie_droite_fct_forme_gauche_1(j));
    demi_somme_gauche_2 = 0.5*(partie_gauche_fct_forme_gauche_2(j+1) + partie_gauche_fct_forme_gauche_2(j));
    demi_somme_droite_2 = 0.5*(partie_droite_fct_forme_gauche_2(j+1) + partie_droite_fct_forme_gauche_2(j));
    Somme_M_sigma_1 = Somme_M_sigma_1 + Q_tilde_g*demi_somme_gauche_1*demi_somme_gauche_2 + Q_tilde_d*demi_somme_droite_1*demi_somme_droite_2;
end

MM_sigma(l,l)= h*Somme_M_sigma_1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
