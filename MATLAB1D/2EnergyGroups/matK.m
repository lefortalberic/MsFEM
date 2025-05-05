function KK = matK(Khi_global,Psi , Psi_adjoint, Nbpt_H, Nbpt_h, epsilon, Num_coeff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matM :
% calcul la matrice de masse
%
% INPUT *
% Num_coeff = 1 ou 2. Num_coeff = 1 Pour calculer la matrice de raideur
% correspondant au premier niveaux d'energie
% OUTPUT - Mel matrice de masse
%
% NOTE (1) calcul de l'int�grale par une quadrature
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KK = sparse(Nbpt_H,Nbpt_H);
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


for l=2:(Nbpt_H-2)   %La premi�re et derni�re ligne doivent rester nulle

    % Khi_(l)
    partie_gauche_fct_forme_droite = Khi_global(:,l,1);
    partie_droite_fct_forme_droite = Khi_global(:,l,2);

    for i=1:Nbpt_h
        xi_g = (l-1-1)*H+(i-1)*h ;
        xi_d = (l-1)*H+(i-1)*h ;
        yi_g = rem(xi_g/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk_g = fix(yi_g*(Nbpt_ref-1)) +1; %partie enti�re de y*(Nbpt-1) (indice pour les matrices)
        yi_d = rem(xi_d/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk_d = fix(yi_d*(Nbpt_ref-1)) +1;
        t_g = yi_g*(Nbpt_ref-1) - fix(yi_g*(Nbpt_ref-1)) ;
        t_d = yi_d*(Nbpt_ref-1) - fix(yi_d*(Nbpt_ref-1)) ;
        Psi_eps_g(i) = Psi(yk_g)*(1-t_g) +Psi(yk_g+1)*t_g;
        Psi_eps_d(i) = Psi(yk_d)*(1-t_d) +Psi(yk_d+1)*t_d;
        Psi_eps_adjoint_g(i) = Psi_adjoint(yk_g)*(1-t_g) +Psi_adjoint(yk_g+1)*t_g;
        Psi_eps_adjoint_d(i) = Psi_adjoint(yk_d)*(1-t_d) +Psi_adjoint(yk_d+1)*t_d;
    end

    Somme_K_1 = 0;
    Somme_K_2 = 0;

    for j=1:(Nbpt_h-1)
        psi_g = 0.5*(Psi_eps_g(j)+Psi_eps_g(j+1));
        psi_d = 0.5*(Psi_eps_d(j)+Psi_eps_d(j+1));
        psi_g_adjoint = 0.5*(Psi_eps_adjoint_g(j)+Psi_eps_adjoint_g(j+1));
        psi_d_adjoint = 0.5*(Psi_eps_adjoint_d(j)+Psi_eps_adjoint_d(j+1));

        Somme_K_1 = Somme_K_1 + psi_g*psi_g_adjoint*A((l-2)*H+(j-0.5)*h,epsilon, Num_coeff)*((partie_gauche_fct_forme_gauche(j+1)-partie_gauche_fct_forme_gauche(j))/h)^2; %gradient de Khi_i au carr�
        Somme_K_1 = Somme_K_1 + psi_d*psi_d_adjoint*A((l-1)*H+(j-0.5)*h,epsilon, Num_coeff)*((partie_droite_fct_forme_gauche(j+1)-partie_droite_fct_forme_gauche(j))/h)^2; %gradient de Khi_i au carr�
    end

    for j=1:(Nbpt_h-1)
        psi_d = 0.5*(Psi_eps_d(j)+Psi_eps_d(j+1));
        psi_d_adjoint = 0.5*(Psi_eps_adjoint_d(j)+Psi_eps_adjoint_d(j+1));
        Somme_K_2 = Somme_K_2 + psi_d*psi_d_adjoint*A((l-1)*H+(j-0.5)*h,epsilon, Num_coeff)*((partie_droite_fct_forme_gauche(j+1)-partie_droite_fct_forme_gauche(j))/h)*((partie_gauche_fct_forme_droite(j+1)-partie_gauche_fct_forme_droite(j))/h); %gradient de Khi_i*gradient de Khi_j
    end

    KK(l,l)= h*Somme_K_1;
    KK(l,l+1)=h*Somme_K_2;
    KK(l+1,l)=h*Somme_K_2;

    % fct_forme_gauche=fct_forme_droite
    partie_gauche_fct_forme_gauche = partie_gauche_fct_forme_droite ;
    partie_droite_fct_forme_gauche = partie_droite_fct_forme_droite ;

end % for l

%Derni�re case des matrices
Somme_K_1 = 0;
l=Nbpt_H-1;

for i=1:Nbpt_h
    xi_g = (l-1-1)*H+(i-1)*h ;
    xi_d = (l-1)*H+(i-1)*h ;
    yi_g = rem(xi_g/epsilon,1) ; % y vaut x/epsilon modulo 1
    yk_g = fix(yi_g*(Nbpt_ref-1)) +1; %partie enti�re de y*(Nbpt-1) (indice pour les matrices)
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

    Somme_K_1 = Somme_K_1 + psi_g*psi_g_adjoint*A((l-2)*H+(j-0.5)*h,epsilon, Num_coeff)*((partie_gauche_fct_forme_gauche(j+1)-partie_gauche_fct_forme_gauche(j))/h)^2; %gradient de Khi_i au carr�
    Somme_K_1 = Somme_K_1 + psi_d*psi_d_adjoint*A((l-1)*H+(j-0.5)*h,epsilon, Num_coeff)*((partie_droite_fct_forme_gauche(j+1)-partie_droite_fct_forme_gauche(j))/h)^2; %gradient de Khi_i au carr�
end

KK(l,l)= h*Somme_K_1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
