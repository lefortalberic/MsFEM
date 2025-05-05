function MM = matM(Khi_global,Psi , Psi_adjoint,Nbpt_H,Nbpt_h,epsilon)
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
Nbpt_ref = (Nbpt_h-1)*(Nbpt_H-1)+1;

Psi_eps_g = zeros(Nbpt_h,1);
Psi_eps_d = zeros(Nbpt_h,1);
Psi_eps_adjoint_g = zeros(Nbpt_h,1);
Psi_eps_adjoint_d = zeros(Nbpt_h,1);

%Khi_1
partie_gauche_fct_forme_gauche = Khi_global(:,1,1);
partie_droite_fct_forme_gauche = Khi_global(:,1,2);

for l=2:(Nbpt_H-2)   %La première et dernière ligne doivent rester nulle

    for i=1:Nbpt_h
        xi_g = (l-1-1)*H+(i-1)*h ;
        xi_d = (l-1)*H+(i-1)*h ;
        yi_g = rem(xi_g/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk_g = fix(yi_g*(Nbpt_ref-1)) +1; %partie entière de y*(Nbpt-1) (indice pour les matrices)
        yi_d = rem(xi_d/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk_d = fix(yi_d*(Nbpt_ref-1)) +1;
        t_g = yi_g*(Nbpt_ref-1) - fix(yi_g*(Nbpt_ref-1)) ;
        t_d = yi_d*(Nbpt_ref-1) - fix(yi_d*(Nbpt_ref-1)) ;
        Psi_eps_g(i) = Psi(yk_g)*(1-t_g) +Psi(yk_g+1)*t_g;
        Psi_eps_d(i) = Psi(yk_d)*(1-t_d) +Psi(yk_d+1)*t_d;
        Psi_eps_adjoint_g(i) = Psi_adjoint(yk_g)*(1-t_g) +Psi_adjoint(yk_g+1)*t_g;
        Psi_eps_adjoint_d(i) = Psi_adjoint(yk_d)*(1-t_d) +Psi_adjoint(yk_d+1)*t_d;
    end

    % Khi_(l)
    partie_gauche_fct_forme_droite = Khi_global(:,l,1);
    partie_droite_fct_forme_droite = Khi_global(:,l,2);

    Somme_M_1 = 0;
    Somme_M_2 = 0;

    for j=1:(Nbpt_h-1)
        psi_g = 0.5*(Psi_eps_g(j)+Psi_eps_g(j+1));
        psi_d = 0.5*(Psi_eps_d(j)+Psi_eps_d(j+1));
        psi_g_adjoint = 0.5*(Psi_eps_adjoint_g(j)+Psi_eps_adjoint_g(j+1));
        psi_d_adjoint = 0.5*(Psi_eps_adjoint_d(j)+Psi_eps_adjoint_d(j+1));

        demi_somme_gauche = 0.5*(partie_gauche_fct_forme_gauche(j+1) + partie_gauche_fct_forme_gauche(j));
        demi_somme_droite = 0.5*(partie_droite_fct_forme_gauche(j+1) + partie_droite_fct_forme_gauche(j));
        Somme_M_1 = Somme_M_1 + psi_g*psi_g_adjoint*demi_somme_gauche^2 + psi_d*psi_d_adjoint*demi_somme_droite^2;
    end

    for j=1:(Nbpt_h-1)
        psi_d = 0.5*(Psi_eps_d(j)+Psi_eps_d(j+1));
        psi_d_adjoint = 0.5*(Psi_eps_adjoint_d(j)+Psi_eps_adjoint_d(j+1));
        demi_somme = psi_d*psi_d_adjoint*0.5*(partie_droite_fct_forme_gauche(j+1)*partie_gauche_fct_forme_droite(j+1) ...
            + partie_droite_fct_forme_gauche(j)*partie_gauche_fct_forme_droite(j));

        Somme_M_2 = Somme_M_2 + demi_somme;
    end

    MM(l,l)= h*Somme_M_1;
    MM(l,l+1)=h*Somme_M_2;
    MM(l+1,l)=h*Somme_M_2;

    % fct_forme_gauche=fct_forme_droite
    partie_gauche_fct_forme_gauche = partie_gauche_fct_forme_droite ;
    partie_droite_fct_forme_gauche = partie_droite_fct_forme_droite ;


end % for l

%Dernière case des matrices
Somme_M_1 = 0;
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
    Psi_eps_g(i) = Psi(yk_g)*(1-t_g) +Psi(yk_g+1)*t_g;
    Psi_eps_d(i) = Psi(yk_d)*(1-t_d) +Psi(yk_d+1)*t_d;
    Psi_eps_adjoint_g(i) = Psi_adjoint(yk_g)*(1-t_g) +Psi_adjoint(yk_g+1)*t_g;
    Psi_eps_adjoint_d(i) = Psi_adjoint(yk_d)*(1-t_d) +Psi_adjoint(yk_d+1)*t_d;
end

for j=1:(Nbpt_h-1)
    psi_g = 0.5*(Psi_eps_g(j)+Psi_eps_g(j+1));
    psi_d = 0.5*(Psi_eps_d(j)+Psi_eps_d(j+1));
    psi_g_adjoint = 0.5*(Psi_eps_adjoint_g(j)+Psi_eps_adjoint_g(j+1));
    psi_d_adjoint = 0.5*(Psi_eps_adjoint_d(j)+Psi_eps_adjoint_d(j+1));

    demi_somme_gauche = 0.5*(partie_gauche_fct_forme_gauche(j+1) + partie_gauche_fct_forme_gauche(j));
    demi_somme_droite = 0.5*(partie_droite_fct_forme_gauche(j+1) + partie_droite_fct_forme_gauche(j));
    Somme_M_1 = Somme_M_1 + psi_g*psi_g_adjoint*demi_somme_gauche^2 + psi_d*psi_d_adjoint*demi_somme_droite^2;
end

MM(l,l)= h*Somme_M_1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
