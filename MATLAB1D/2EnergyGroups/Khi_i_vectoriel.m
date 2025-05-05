function [ UU_g_cas1_1, UU_g_cas1_2 , UU_d_cas1_1, UU_d_cas1_2, UU_g_cas2_1, UU_g_cas2_2 , UU_d_cas2_1, UU_d_cas2_2] = Khi_i_vectoriel(Nbpt_H,Nbpt_h,epsilon,noeud_i,num_idee)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fonction de forme finale, combinaison de la solution du pb
% aux valeurs propres et des fonctions de base P1

% renvoie une matrice de taille Nbpt_h * 2.
%La première colonne (  (:,1)  ) est la partie gauche de la fonction de
%forme associé au noeud_i. La 2ieme est la partie droite.



% INPUT - %noeud_i : qui va de 1 à Nbpt_H-1. (Khi_0 est nulle)
%       - %num_idee : numéro de l'idée utilisée pour la génération des fonctions de forme
%                       num_idee = 6  --> Fonction de forme P1.
%                       num_idee = 8  --> Fonction de forme MsFEM triche, du problème en v_eps, avec Psi_eps^2*D(x/eps) comme
%                       matrice de diffusion. Fonction de forme multigroup(=2) pour résoudre le problème en v_eps
%                       num_idee = 9  --> Fonction de forme MsFEM triche, du problème en u_eps, avec Psi_eps^2*D(x/eps) comme
%                       matrice de diffusion. Je multiplie ensuite chaque fonction de forme par Psi. Fonction de forme
%                       multigroup(=2) pour résoudre le problème en u_eps.
%                       num_idee = 10  --> Fonction de forme MsFEM filtre, du problème en v_eps, avec Phi_eps^2*D(x/eps) comme
%                       matrice de diffusion. Je multiplie ensuite chaque fonction de forme par Phi. Fonction de forme
%                       multigroup(=2) pour résoudre le problème en v_eps.
%                       Je calcul les Phi_eps par la méthode du filtre.
%                       num_idee = 11  --> Fonction de forme MsFEM filtre, du problème en u_eps, avec Phi_eps^2*D(x/eps) comme
%                       matrice de diffusion. Je multiplie ensuite chaque fonction de forme par Phi. Fonction de forme
%                       multigroup(=2) pour résoudre le problème en u_eps.
%                       Je calcul les Phi_eps par la méthode du filtre.
%                       num_idee = 12  --> Fonction de forme MsFEM triche, du problème en u_eps, avec Psi_eps^2*D(x/eps) comme
%                       matrice de diffusion. Je multiplie ensuite chaque fonction de forme par Psi. Fonction de forme
%                       multigroup(=2) pour résoudre le problème en u_eps.
%                       La difference avec 8 est que j'impose les CL
%                       séparémment.

% OUTPUT - val:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H=1/(Nbpt_H -1);
h = 1/(Nbpt_h -1)*H;
Nbpt_ref = (Nbpt_h-1)*(Nbpt_H-1)+1;
h_ref = 1/(Nbpt_ref-1);
%Maillage
XX = zeros(Nbpt_h,1);     %vecteur des coordonnées des noeuds
for j=1:Nbpt_h
    XX(j)=(j-1)*h;
end
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if num_idee==12
    
    Psi_1_eps_g = zeros(Nbpt_h,1); Psi_2_eps_g = zeros(Nbpt_h,1);
    Psi_1_eps_adjoint_g = zeros(Nbpt_h,1); Psi_2_eps_adjoint_g = zeros(Nbpt_h,1);
    Q_tilde_11_eps_g = zeros(Nbpt_h,1);
    Q_tilde_22_eps_g = zeros(Nbpt_h,1);
    Q_tilde_21_eps_g = zeros(Nbpt_h,1);
    Q_tilde_12_eps_g = zeros(Nbpt_h,1);
    vecteur_J1_eps_g = zeros(Nbpt_h,1);
    vecteur_J2_eps_g = zeros(Nbpt_h,1);
    
    Psi_1_eps_d = zeros(Nbpt_h,1); Psi_2_eps_d = zeros(Nbpt_h,1);
    Psi_1_eps_adjoint_d = zeros(Nbpt_h,1); Psi_2_eps_adjoint_d = zeros(Nbpt_h,1);
    Q_tilde_11_eps_d = zeros(Nbpt_h,1);
    Q_tilde_22_eps_d = zeros(Nbpt_h,1);
    Q_tilde_21_eps_d = zeros(Nbpt_h,1);
    Q_tilde_12_eps_d = zeros(Nbpt_h,1);
    vecteur_J1_eps_d = zeros(Nbpt_h,1);
    vecteur_J2_eps_d = zeros(Nbpt_h,1);
    
    [~,Psi_1_adjoint,Psi_2_adjoint] = Solution_pb_spectral_adjoint(Nbpt_h);
    [~,Psi_1,Psi_2] = Solution_pb_spectral(Nbpt_h);
    
    Q_tilde_11 = Q_tilde(Psi_1,Psi_2,Psi_1_adjoint,Psi_2_adjoint,1);
    Q_tilde_22 = Q_tilde(Psi_1,Psi_2,Psi_1_adjoint,Psi_2_adjoint,4);
    Q_tilde_21 = - Q_tilde_22;
    Q_tilde_12 = - Q_tilde_11;
    
    vecteur_J1 = J(Psi_1,Psi_1_adjoint,1);
    vecteur_J2 = J(Psi_2,Psi_2_adjoint,2);
    
    for i=1:Nbpt_h
        xi_g = (noeud_i-1)*H+(i-1)*h ;
        xi_d = (noeud_i)*H+(i-1)*h ;
        yi_g = rem(xi_g/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk_g = fix(yi_g*(Nbpt_h-1)) +1; %partie entière de y*(Nbpt-1) (indice pour les matrices)
        yi_d = rem(xi_d/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk_d = fix(yi_d*(Nbpt_h-1)) +1;
        t_g = yi_g*(Nbpt_h-1) - fix(yi_g*(Nbpt_h-1)) ;
        Psi_1_eps_g(i) = Psi_1(yk_g)*(1-t_g) +Psi_1(yk_g+1)*t_g;
        Psi_2_eps_g(i) = Psi_2(yk_g)*(1-t_g) +Psi_2(yk_g+1)*t_g;
        Psi_1_eps_adjoint_g(i) = Psi_1_adjoint(yk_g)*(1-t_g) +Psi_1_adjoint(yk_g+1)*t_g;
        Psi_2_eps_adjoint_g(i) = Psi_2_adjoint(yk_g)*(1-t_g) +Psi_2_adjoint(yk_g+1)*t_g;
        Q_tilde_11_eps_g(i) = Q_tilde_11(yk_g)*(1-t_g) +Q_tilde_11(yk_g+1)*t_g;
        Q_tilde_22_eps_g(i) = Q_tilde_22(yk_g)*(1-t_g) +Q_tilde_22(yk_g+1)*t_g;
        Q_tilde_21_eps_g(i) = Q_tilde_21(yk_g)*(1-t_g) +Q_tilde_21(yk_g+1)*t_g;
        Q_tilde_12_eps_g(i) = Q_tilde_12(yk_g)*(1-t_g) +Q_tilde_12(yk_g+1)*t_g;
        vecteur_J1_eps_g(i) = vecteur_J1(yk_g)*(1-t_g) +vecteur_J1(yk_g+1)*t_g;
        vecteur_J2_eps_g(i) = vecteur_J2(yk_g)*(1-t_g) +vecteur_J2(yk_g+1)*t_g;
    
        t_d = yi_d*(Nbpt_h-1) - fix(yi_d*(Nbpt_h-1)) ;
        Psi_1_eps_d(i) = Psi_1(yk_d)*(1-t_d) +Psi_1(yk_d+1)*t_d;
        Psi_2_eps_d(i) = Psi_2(yk_d)*(1-t_d) +Psi_2(yk_d+1)*t_d;
        Psi_1_eps_adjoint_d(i) = Psi_1_adjoint(yk_d)*(1-t_d) +Psi_1_adjoint(yk_d+1)*t_d;
        Psi_2_eps_adjoint_d(i) = Psi_2_adjoint(yk_d)*(1-t_d) +Psi_2_adjoint(yk_d+1)*t_d;
        Q_tilde_11_eps_d(i) = Q_tilde_11(yk_d)*(1-t_d) +Q_tilde_11(yk_d+1)*t_d;
        Q_tilde_22_eps_d(i) = Q_tilde_22(yk_d)*(1-t_d) +Q_tilde_22(yk_d+1)*t_d;
        Q_tilde_21_eps_d(i) = Q_tilde_21(yk_d)*(1-t_d) +Q_tilde_21(yk_d+1)*t_d;
        Q_tilde_12_eps_d(i) = Q_tilde_12(yk_d)*(1-t_d) +Q_tilde_12(yk_d+1)*t_d;
        vecteur_J1_eps_d(i) = vecteur_J1(yk_d)*(1-t_d) +vecteur_J1(yk_d+1)*t_d;
        vecteur_J2_eps_d(i) = vecteur_J2(yk_d)*(1-t_d) +vecteur_J2(yk_d+1)*t_d;
    end

    Zeros = sparse(Nbpt_h,Nbpt_h); % matrice de masse
    
    KK_1_g = matK_ref_Psi(Nbpt_h,epsilon, Psi_1_eps_g , Psi_1_eps_adjoint_g, 1, (noeud_i-1)*H, h);
    KK_2_g = matK_ref_Psi(Nbpt_h,epsilon, Psi_2_eps_g , Psi_2_eps_adjoint_g, 2, (noeud_i-1)*H, h);
    
    KK_1_d = matK_ref_Psi(Nbpt_h,epsilon, Psi_1_eps_d , Psi_1_eps_adjoint_d, 1, (noeud_i)*H, h);
    KK_2_d = matK_ref_Psi(Nbpt_h,epsilon, Psi_2_eps_d , Psi_2_eps_adjoint_d, 2, (noeud_i)*H, h);
    
    MM_Q_tilde_11_g = matM_ref_Psi(Nbpt_h, Q_tilde_11_eps_g , ones(Nbpt_h,1),h);
    MM_Q_tilde_22_g = matM_ref_Psi(Nbpt_h, Q_tilde_22_eps_g , ones(Nbpt_h,1),h);
    MM_Q_tilde_21_g = matM_ref_Psi(Nbpt_h, Q_tilde_21_eps_g , ones(Nbpt_h,1),h);
    MM_Q_tilde_12_g = matM_ref_Psi(Nbpt_h, Q_tilde_12_eps_g , ones(Nbpt_h,1),h);
    
    MM_Q_tilde_11_d = matM_ref_Psi(Nbpt_h, Q_tilde_11_eps_d , ones(Nbpt_h,1),h);
    MM_Q_tilde_22_d = matM_ref_Psi(Nbpt_h, Q_tilde_22_eps_d , ones(Nbpt_h,1),h);
    MM_Q_tilde_21_d = matM_ref_Psi(Nbpt_h, Q_tilde_21_eps_d , ones(Nbpt_h,1),h);
    MM_Q_tilde_12_d = matM_ref_Psi(Nbpt_h, Q_tilde_12_eps_d , ones(Nbpt_h,1),h);
    
    CC_1_g = matC_ref(Nbpt_h, vecteur_J1_eps_g);
    CC_2_g = matC_ref(Nbpt_h, vecteur_J2_eps_g);
    
    CC_1_d = matC_ref(Nbpt_h, vecteur_J1_eps_d);
    CC_2_d = matC_ref(Nbpt_h, vecteur_J2_eps_d);

    % Calcul du second membre F
    % Cas 1 correspond à Chi_1 avec CL 0-1 et Chi_2 avec CL 0-0 (pour le
    % calcul des fonctions forme gauche)
    % Cas 0 correspond à Chi_2 avec CL 0-0 et Chi_2 avec CL 0-1 (pour le
    % calcul des fonctions forme gauche) 
    % -------------------------
    LL_1_div_g = (-1/H)*matF_div(Nbpt_h, epsilon, Psi_1_eps_g , Psi_1_eps_adjoint_g, 1, (noeud_i-1)*H, h);     % vecteur second membre
    LL_1_div_d = (1/H)*matF_div(Nbpt_h, epsilon, Psi_1_eps_d , Psi_1_eps_adjoint_d, 1, (noeud_i)*H, h);     % vecteur second membre
    LL_2_div_g = (-1/H)*matF_div(Nbpt_h, epsilon, Psi_2_eps_g , Psi_2_eps_adjoint_g, 2, (noeud_i-1)*H, h);     % vecteur second membre
    LL_2_div_d = (1/H)*matF_div(Nbpt_h, epsilon, Psi_2_eps_d , Psi_2_eps_adjoint_d, 2, (noeud_i)*H, h);     % vecteur second membre

    LL_1_J_g = (1/H)*matF_J(Nbpt_h, vecteur_J1_eps_g,h);
    LL_1_J_d = (-1/H)*matF_J(Nbpt_h, vecteur_J1_eps_d,h);
    LL_2_J_g = (1/H)*matF_J(Nbpt_h, vecteur_J2_eps_g,h);
    LL_2_J_d = (-1/H)*matF_J(Nbpt_h, vecteur_J2_eps_d,h);

    Phi_P1_g = zeros(Nbpt_h,1);
    Phi_P1_d = zeros(Nbpt_h,1);
    for l=1:Nbpt_h   %La première et dernière ligne doivent rester nulle
        Phi_P1_g(l) = (l-1)*h/H;
        Phi_P1_d(l) = 1-(l-1)*h/H;
    end % for l

%     LL_1_Q_tilde_g = matF_J(Nbpt_h, (Q_tilde_11_eps_g + Q_tilde_12_eps_g).*Phi_P1_g,h);
%     LL_1_Q_tilde_d = matF_J(Nbpt_h, (Q_tilde_11_eps_d + Q_tilde_12_eps_d).*Phi_P1_d,h);
%     LL_2_Q_tilde_g = matF_J(Nbpt_h, (Q_tilde_21_eps_g + Q_tilde_22_eps_g).*Phi_P1_g,h);
%     LL_2_Q_tilde_d = matF_J(Nbpt_h, (Q_tilde_21_eps_d + Q_tilde_22_eps_d).*Phi_P1_d,h);
    LL_1_Q_tilde_g_cas1 = matF_J(Nbpt_h, (Q_tilde_11_eps_g).*Phi_P1_g,h);
    LL_2_Q_tilde_g_cas1 = matF_J(Nbpt_h, (Q_tilde_21_eps_g).*Phi_P1_g,h);
    LL_1_Q_tilde_g_cas2 = matF_J(Nbpt_h, (Q_tilde_12_eps_g).*Phi_P1_g,h);
    LL_2_Q_tilde_g_cas2 = matF_J(Nbpt_h, (Q_tilde_22_eps_g).*Phi_P1_g,h);
    
    LL_1_Q_tilde_d_cas1 = matF_J(Nbpt_h, (Q_tilde_11_eps_d).*Phi_P1_d,h);
    LL_2_Q_tilde_d_cas1 = matF_J(Nbpt_h, (Q_tilde_21_eps_d).*Phi_P1_d,h);
    LL_1_Q_tilde_d_cas2 = matF_J(Nbpt_h, (Q_tilde_12_eps_d).*Phi_P1_d,h);
    LL_2_Q_tilde_d_cas2 = matF_J(Nbpt_h, (Q_tilde_22_eps_d).*Phi_P1_d,h);
    
    %Assemblage de toutes les matrices
    KK_g = [KK_1_g , Zeros ; Zeros ,KK_2_g ];
    MM_Q_tilde_g = [MM_Q_tilde_11_g , MM_Q_tilde_12_g ; MM_Q_tilde_21_g , MM_Q_tilde_22_g];
    CC_g = [CC_1_g , Zeros ; Zeros ,CC_2_g];

    LL_div_g_cas1 = [LL_1_div_g ; 0.*LL_2_div_g];
    LL_J_g_cas1 = [LL_1_J_g ; 0.*LL_2_J_g];
    LL_Q_tilde_g_cas1 = [LL_1_Q_tilde_g_cas1 ; LL_2_Q_tilde_g_cas1];

    LL_g_cas1 = LL_div_g_cas1 - LL_Q_tilde_g_cas1/(epsilon^2) - LL_J_g_cas1/epsilon;

    LL_div_g_cas2 = [0.*LL_1_div_g ; LL_2_div_g];
    LL_J_g_cas2 = [0.*LL_1_J_g ; LL_2_J_g];
    LL_Q_tilde_g_cas2 = [LL_1_Q_tilde_g_cas2 ; LL_2_Q_tilde_g_cas2];

    LL_g_cas2 = LL_div_g_cas2 - LL_Q_tilde_g_cas2/(epsilon^2) - LL_J_g_cas2/epsilon;


    LL_div_d_cas1 = [LL_1_div_d ; 0.*LL_2_div_d];
    LL_J_d_cas1 = [LL_1_J_g ; 0.*LL_2_J_g];
    LL_Q_tilde_d_cas1 = [LL_1_Q_tilde_d_cas1 ; LL_2_Q_tilde_d_cas1];

    LL_d_cas1 = LL_div_d_cas1 - LL_Q_tilde_d_cas1/(epsilon^2) - LL_J_d_cas1/epsilon;

    LL_div_d_cas2 = [0.*LL_1_div_d ; LL_2_div_d];
    LL_J_d_cas2 = [0.*LL_1_J_d ; LL_2_J_d];
    LL_Q_tilde_d_cas2 = [LL_1_Q_tilde_d_cas2 ; LL_2_Q_tilde_d_cas2];

    LL_d_cas2 = LL_div_d_cas2 - LL_Q_tilde_d_cas2/(epsilon^2) - LL_J_d_cas2/epsilon;

    AA_g = KK_g + MM_Q_tilde_g/(epsilon^2) + CC_g/epsilon;
    AA_g(1,1)=1; AA_g(2*Nbpt_h,2*Nbpt_h)=1;
    AA_g(Nbpt_h,Nbpt_h)=1 ; AA_g(Nbpt_h+1,Nbpt_h+1)=1;

    KK_d = [KK_1_d , Zeros ; Zeros ,KK_2_d ];
    MM_Q_tilde_d = [MM_Q_tilde_11_d , MM_Q_tilde_12_d ; MM_Q_tilde_21_d , MM_Q_tilde_22_d];
    CC_d = [CC_1_d , Zeros ; Zeros ,CC_2_d];

    AA_d = KK_d + MM_Q_tilde_d/(epsilon^2) + CC_d/epsilon;
    AA_d(1,1)=1; AA_d(2*Nbpt_h,2*Nbpt_h)=1;
    AA_d(Nbpt_h,Nbpt_h)=1 ; AA_d(Nbpt_h+1,Nbpt_h+1)=1;

    
    UU_g_cas1 = AA_g\LL_g_cas1;
    UU_g_cas2 = AA_g\LL_g_cas2;
    UU_d_cas1 = AA_d\LL_d_cas1;
    UU_d_cas2 = AA_d\LL_d_cas2;

    UU_g_cas1_1 = UU_g_cas1(1:Nbpt_h);
    UU_g_cas1_2 = UU_g_cas1((Nbpt_h+1):2*Nbpt_h);
    UU_g_cas2_1 = UU_g_cas2(1:Nbpt_h);
    UU_g_cas2_2 = UU_g_cas2((Nbpt_h+1):2*Nbpt_h);

    UU_d_cas1_1 = UU_d_cas1(1:Nbpt_h);
    UU_d_cas1_2 = UU_d_cas1((Nbpt_h+1):2*Nbpt_h);
    UU_d_cas2_1 = UU_d_cas2(1:Nbpt_h);
    UU_d_cas2_2 = UU_d_cas2((Nbpt_h+1):2*Nbpt_h);

    % réhaussement
    % ----------
    for l=1:Nbpt_h

        UU_g_cas1_1(l)=UU_g_cas1_1(l)+(l-1)*(h/H);
        %UU_g_cas1_2(l)=UU_g_cas1_2(l)+(l-1)*(h/H); %Lui n'est pas releve dans le cas1
        %UU_g_cas2_1(l)=UU_g_cas1_1(l)+(l-1)*(h/H); %Lui n'est pas releve dans le cas1
        UU_g_cas2_2(l)=UU_g_cas2_2(l)+(l-1)*(h/H);

        UU_d_cas1_1(l)=UU_d_cas1_1(l)-(l-1)*(h/H)+1;
        UU_d_cas2_2(l)=UU_d_cas2_2(l)-(l-1)*(h/H)+1;

    end % for l

    %Multiplication par Psi_eps
    for l=1:Nbpt_h

        UU_g_cas1_1(l)=UU_g_cas1_1(l)*Psi_1_eps_g(l);
        UU_g_cas1_2(l) = UU_g_cas1_2(l)*Psi_2_eps_g(l);
        UU_g_cas2_1(l) = UU_g_cas2_1(l)*Psi_1_eps_g(l);
        UU_g_cas2_2(l) = UU_g_cas2_2(l)*Psi_2_eps_g(l);
    
        UU_d_cas1_1(l) = UU_d_cas1_1(l)*Psi_1_eps_d(l);
        UU_d_cas1_2(l) = UU_d_cas1_2(l)*Psi_2_eps_d(l);
        UU_d_cas2_1(l) = UU_d_cas2_1(l)*Psi_1_eps_d(l);
        UU_d_cas2_2(l) = UU_d_cas2_2(l)*Psi_2_eps_d(l);

    end % for l
    %val = [ UU_g , UU_d];

elseif num_idee==9

    Psi_1_eps_g = zeros(Nbpt_h,1); Psi_2_eps_g = zeros(Nbpt_h,1);
    Psi_1_eps_adjoint_g = zeros(Nbpt_h,1); Psi_2_eps_adjoint_g = zeros(Nbpt_h,1);
    Q_tilde_11_eps_g = zeros(Nbpt_h,1);
    Q_tilde_22_eps_g = zeros(Nbpt_h,1);
    Q_tilde_21_eps_g = zeros(Nbpt_h,1);
    Q_tilde_12_eps_g = zeros(Nbpt_h,1);
    vecteur_J1_eps_g = zeros(Nbpt_h,1);
    vecteur_J2_eps_g = zeros(Nbpt_h,1);
    
    Psi_1_eps_d = zeros(Nbpt_h,1); Psi_2_eps_d = zeros(Nbpt_h,1);
    Psi_1_eps_adjoint_d = zeros(Nbpt_h,1); Psi_2_eps_adjoint_d = zeros(Nbpt_h,1);
    Q_tilde_11_eps_d = zeros(Nbpt_h,1);
    Q_tilde_22_eps_d = zeros(Nbpt_h,1);
    Q_tilde_21_eps_d = zeros(Nbpt_h,1);
    Q_tilde_12_eps_d = zeros(Nbpt_h,1);
    vecteur_J1_eps_d = zeros(Nbpt_h,1);
    vecteur_J2_eps_d = zeros(Nbpt_h,1);
    
    [~,Psi_1_adjoint,Psi_2_adjoint] = Solution_pb_spectral_adjoint(Nbpt_h);
    [~,Psi_1,Psi_2] = Solution_pb_spectral(Nbpt_h);
    
    Q_tilde_11 = Q_tilde(Psi_1,Psi_2,Psi_1_adjoint,Psi_2_adjoint,1);
    Q_tilde_22 = Q_tilde(Psi_1,Psi_2,Psi_1_adjoint,Psi_2_adjoint,4);
    Q_tilde_21 = - Q_tilde_22;
    Q_tilde_12 = - Q_tilde_11;
    
    vecteur_J1 = J(Psi_1,Psi_1_adjoint,1);
    vecteur_J2 = J(Psi_2,Psi_2_adjoint,2);
    
    for i=1:Nbpt_h
        xi_g = (noeud_i-1)*H+(i-1)*h ;
        xi_d = (noeud_i)*H+(i-1)*h ;
        yi_g = rem(xi_g/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk_g = fix(yi_g*(Nbpt_h-1)) +1; %partie entière de y*(Nbpt-1) (indice pour les matrices)
        yi_d = rem(xi_d/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk_d = fix(yi_d*(Nbpt_h-1)) +1;
        t_g = yi_g*(Nbpt_h-1) - fix(yi_g*(Nbpt_h-1)) ;
        Psi_1_eps_g(i) = Psi_1(yk_g)*(1-t_g) +Psi_1(yk_g+1)*t_g;
        Psi_2_eps_g(i) = Psi_2(yk_g)*(1-t_g) +Psi_2(yk_g+1)*t_g;
        Psi_1_eps_adjoint_g(i) = Psi_1_adjoint(yk_g)*(1-t_g) +Psi_1_adjoint(yk_g+1)*t_g;
        Psi_2_eps_adjoint_g(i) = Psi_2_adjoint(yk_g)*(1-t_g) +Psi_2_adjoint(yk_g+1)*t_g;
        Q_tilde_11_eps_g(i) = Q_tilde_11(yk_g)*(1-t_g) +Q_tilde_11(yk_g+1)*t_g;
        Q_tilde_22_eps_g(i) = Q_tilde_22(yk_g)*(1-t_g) +Q_tilde_22(yk_g+1)*t_g;
        Q_tilde_21_eps_g(i) = Q_tilde_21(yk_g)*(1-t_g) +Q_tilde_21(yk_g+1)*t_g;
        Q_tilde_12_eps_g(i) = Q_tilde_12(yk_g)*(1-t_g) +Q_tilde_12(yk_g+1)*t_g;
        vecteur_J1_eps_g(i) = vecteur_J1(yk_g)*(1-t_g) +vecteur_J1(yk_g+1)*t_g;
        vecteur_J2_eps_g(i) = vecteur_J2(yk_g)*(1-t_g) +vecteur_J2(yk_g+1)*t_g;
    
        t_d = yi_d*(Nbpt_h-1) - fix(yi_d*(Nbpt_h-1)) ;
        Psi_1_eps_d(i) = Psi_1(yk_d)*(1-t_d) +Psi_1(yk_d+1)*t_d;
        Psi_2_eps_d(i) = Psi_2(yk_d)*(1-t_d) +Psi_2(yk_d+1)*t_d;
        Psi_1_eps_adjoint_d(i) = Psi_1_adjoint(yk_d)*(1-t_d) +Psi_1_adjoint(yk_d+1)*t_d;
        Psi_2_eps_adjoint_d(i) = Psi_2_adjoint(yk_d)*(1-t_d) +Psi_2_adjoint(yk_d+1)*t_d;
        Q_tilde_11_eps_d(i) = Q_tilde_11(yk_d)*(1-t_d) +Q_tilde_11(yk_d+1)*t_d;
        Q_tilde_22_eps_d(i) = Q_tilde_22(yk_d)*(1-t_d) +Q_tilde_22(yk_d+1)*t_d;
        Q_tilde_21_eps_d(i) = Q_tilde_21(yk_d)*(1-t_d) +Q_tilde_21(yk_d+1)*t_d;
        Q_tilde_12_eps_d(i) = Q_tilde_12(yk_d)*(1-t_d) +Q_tilde_12(yk_d+1)*t_d;
        vecteur_J1_eps_d(i) = vecteur_J1(yk_d)*(1-t_d) +vecteur_J1(yk_d+1)*t_d;
        vecteur_J2_eps_d(i) = vecteur_J2(yk_d)*(1-t_d) +vecteur_J2(yk_d+1)*t_d;
    end
    
    Zeros = sparse(Nbpt_h,Nbpt_h); % matrice de masse
    
    KK_1_g = matK_ref_Psi(Nbpt_h,epsilon, Psi_1_eps_g , Psi_1_eps_adjoint_g, 1, (noeud_i-1)*H, h);
    KK_2_g = matK_ref_Psi(Nbpt_h,epsilon, Psi_2_eps_g , Psi_2_eps_adjoint_g, 2, (noeud_i-1)*H, h);
    
    KK_1_d = matK_ref_Psi(Nbpt_h,epsilon, Psi_1_eps_d , Psi_1_eps_adjoint_d, 1, (noeud_i)*H, h);
    KK_2_d = matK_ref_Psi(Nbpt_h,epsilon, Psi_2_eps_d , Psi_2_eps_adjoint_d, 2, (noeud_i)*H, h);
    
    MM_Q_tilde_11_g = matM_ref_Psi(Nbpt_h, Q_tilde_11_eps_g , ones(Nbpt_h,1),h);
    MM_Q_tilde_22_g = matM_ref_Psi(Nbpt_h, Q_tilde_22_eps_g , ones(Nbpt_h,1),h);
    MM_Q_tilde_21_g = matM_ref_Psi(Nbpt_h, Q_tilde_21_eps_g , ones(Nbpt_h,1),h);
    MM_Q_tilde_12_g = matM_ref_Psi(Nbpt_h, Q_tilde_12_eps_g , ones(Nbpt_h,1),h);
    
    MM_Q_tilde_11_d = matM_ref_Psi(Nbpt_h, Q_tilde_11_eps_d , ones(Nbpt_h,1),h);
    MM_Q_tilde_22_d = matM_ref_Psi(Nbpt_h, Q_tilde_22_eps_d , ones(Nbpt_h,1),h);
    MM_Q_tilde_21_d = matM_ref_Psi(Nbpt_h, Q_tilde_21_eps_d , ones(Nbpt_h,1),h);
    MM_Q_tilde_12_d = matM_ref_Psi(Nbpt_h, Q_tilde_12_eps_d , ones(Nbpt_h,1),h);
    
    CC_1_g = matC_ref(Nbpt_h, vecteur_J1_eps_g);
    CC_2_g = matC_ref(Nbpt_h, vecteur_J2_eps_g);
    
    CC_1_d = matC_ref(Nbpt_h, vecteur_J1_eps_d);
    CC_2_d = matC_ref(Nbpt_h, vecteur_J2_eps_d);

    % Calcul du second membre F
    % -------------------------
    LL_1_div_g = (-1/H)*matF_div(Nbpt_h, epsilon, Psi_1_eps_g , Psi_1_eps_adjoint_g, 1, (noeud_i-1)*H, h);     % vecteur second membre
    LL_1_div_d = (1/H)*matF_div(Nbpt_h, epsilon, Psi_1_eps_d , Psi_1_eps_adjoint_d, 1, (noeud_i)*H, h);     % vecteur second membre
    LL_2_div_g = (-1/H)*matF_div(Nbpt_h, epsilon, Psi_2_eps_g , Psi_2_eps_adjoint_g, 2, (noeud_i-1)*H, h);     % vecteur second membre
    LL_2_div_d = (1/H)*matF_div(Nbpt_h, epsilon, Psi_2_eps_d , Psi_2_eps_adjoint_d, 2, (noeud_i)*H, h);     % vecteur second membre

    LL_1_J_g = (1/H)*matF_J(Nbpt_h, vecteur_J1_eps_g,h);
    LL_1_J_d = (-1/H)*matF_J(Nbpt_h, vecteur_J1_eps_d,h);
    LL_2_J_g = (1/H)*matF_J(Nbpt_h, vecteur_J2_eps_g,h);
    LL_2_J_d = (-1/H)*matF_J(Nbpt_h, vecteur_J2_eps_d,h);

    Phi_P1_g = zeros(Nbpt_h,1);
    Phi_P1_d = zeros(Nbpt_h,1);
    for l=1:Nbpt_h   %La première et dernière ligne doivent rester nulle
        Phi_P1_g(l) = (l-1)*h/H;
        Phi_P1_d(l) = 1-(l-1)*h/H;
    end % for l

    LL_1_Q_tilde_g = matF_J(Nbpt_h, (Q_tilde_11_eps_g + Q_tilde_12_eps_g).*Phi_P1_g,h);
    LL_1_Q_tilde_d = matF_J(Nbpt_h, (Q_tilde_11_eps_d + Q_tilde_12_eps_d).*Phi_P1_d,h);
    LL_2_Q_tilde_g = matF_J(Nbpt_h, (Q_tilde_21_eps_g + Q_tilde_22_eps_g).*Phi_P1_g,h);
    LL_2_Q_tilde_d = matF_J(Nbpt_h, (Q_tilde_21_eps_d + Q_tilde_22_eps_d).*Phi_P1_d,h);
    
    %Assemblage de toutes les matrices
    KK_g = [KK_1_g , Zeros ; Zeros ,KK_2_g ];
    MM_Q_tilde_g = [MM_Q_tilde_11_g , MM_Q_tilde_12_g ; MM_Q_tilde_21_g , MM_Q_tilde_22_g];
    CC_g = [CC_1_g , Zeros ; Zeros ,CC_2_g];

    LL_div_g = [LL_1_div_g ; LL_2_div_g];
    LL_J_g = [LL_1_J_g ; LL_2_J_g];
    LL_Q_tilde_g = [LL_1_Q_tilde_g ; LL_2_Q_tilde_g];

    LL_g = LL_div_g - LL_Q_tilde_g/(epsilon^2) - LL_J_g/epsilon;

    AA_g = KK_g + MM_Q_tilde_g/(epsilon^2) + CC_g/epsilon;
    AA_g(1,1)=1; AA_g(2*Nbpt_h,2*Nbpt_h)=1;
    AA_g(Nbpt_h,Nbpt_h)=1 ; AA_g(Nbpt_h+1,Nbpt_h+1)=1;

    KK_d = [KK_1_d , Zeros ; Zeros ,KK_2_d ];
    MM_Q_tilde_d = [MM_Q_tilde_11_d , MM_Q_tilde_12_d ; MM_Q_tilde_21_d , MM_Q_tilde_22_d];
    CC_d = [CC_1_d , Zeros ; Zeros ,CC_2_d];

    LL_div_d = [LL_1_div_d ; LL_2_div_d];
    LL_J_d = [LL_1_J_d ; LL_2_J_d];
    LL_Q_tilde_d = [LL_1_Q_tilde_d ; LL_2_Q_tilde_d];

    LL_d = LL_div_d - LL_Q_tilde_d/(epsilon^2) - LL_J_d/epsilon;

    AA_d = KK_d + MM_Q_tilde_d/(epsilon^2) + CC_d/epsilon;
    AA_d(1,1)=1; AA_d(2*Nbpt_h,2*Nbpt_h)=1;
    AA_d(Nbpt_h,Nbpt_h)=1 ; AA_d(Nbpt_h+1,Nbpt_h+1)=1;

    
    UU_g = AA_g\LL_g;
    UU_d = AA_d\LL_d;

    UU_g_1 = UU_g(1:Nbpt_h);
    UU_g_2 = UU_g((Nbpt_h+1):2*Nbpt_h);
    UU_d_1 = UU_d(1:Nbpt_h);
    UU_d_2 = UU_d((Nbpt_h+1):2*Nbpt_h);

    % réhaussement
    % ----------
    for l=1:Nbpt_h

        UU_g_1(l)=UU_g_1(l)+(l-1)*(h/H);
        UU_g_2(l)=UU_g_2(l)+(l-1)*(h/H);
        UU_d_1(l)=UU_d_1(l)-(l-1)*(h/H)+1;
        UU_d_2(l)=UU_d_2(l)-(l-1)*(h/H)+1;

    end % for l

    %Multiplication par Psi_eps
    for l=1:Nbpt_h

        UU_g_1(l)=UU_g_1(l)*Psi_1_eps_g(l);
        UU_d_1(l)=UU_d_1(l)*Psi_1_eps_d(l);
        UU_g_2(l)=UU_g_2(l)*Psi_2_eps_g(l);
        UU_d_2(l)=UU_d_2(l)*Psi_2_eps_d(l);

    end % for l

    UU_g_cas1_1 = UU_g_1;
    UU_g_cas1_2 = UU_g_1*0;
    UU_g_cas2_1 = UU_g_1*0;
    UU_g_cas2_2 = UU_g_2 ;

    UU_d_cas1_1 = UU_d_1;
    UU_d_cas1_2 = UU_d_1*0;
    UU_d_cas2_1 = UU_d_1*0;
    UU_d_cas2_2 = UU_d_2 ;
    %val = [ UU_g , UU_d];
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
