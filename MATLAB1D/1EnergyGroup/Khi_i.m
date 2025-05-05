function [ UU_g , UU_d] = Khi_i(Nbpt_H,Nbpt_h,epsilon,noeud_i,num_idee)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fonction de forme finale, combinaison de la solution du pb
% aux valeurs propres et des fonctions de base P1

% renvoie une matrice de taille Nbpt_h * 2.
%La première colonne (  (:,1)  ) est la partie gauche de la fonction de
%forme associé au noeud_i. La 2ieme est la partie droite.



% INPUT - %noeud_i : qui va de 1 à Nbpt_H-1. (Khi_0 est nulle)
%       - %num_idee : numéro de l'idée utilisée pour la génération des fonctions de forme
%                       num_idee = 0  --> Fonction de forme MsFEM_0, du problème purement diffusif
%                       num_idee = 1  --> Fonction de forme MsFEM_1, du problème diffusif avec réaction
%                       num_idee = 2  --> Fonction de forme MsFEM_2, du problème aux valeurs propres (avec filtre PH).
%                       Fonction continue par addition de fonction affine.
%                       num_idee = 3  --> Fonction de forme MsFEM_2, du problème aux valeurs propres (avec filtre PH).
%                       Fonction continue par multiplication de fonction affine.
%                       num_idee = 4  --> Fonction de forme MsFEM_2, du problème aux valeurs propres (sans filtre PH).
%                       Fonction continue par multiplication des fonctions P1.
%                       num_idee = 5  --> Fonction de forme MsFEM_2, du problème aux valeurs propres (sans filtre PH).
%                       Fonction discontinue par multiplication des fonctions P1.
%                       num_idee = 6  --> Fonction de forme P1 classique, qui sert de cas test.
%                       num_idee = 7  --> Fonction de forme MsFEM, du problème purement diffusif, avec Phi_eps^2*D(x/eps) comme
%                       matrice de diffusion, où Phi_eps solution du problème stationnaire avec CL Neumann homogène (avec filtre PH).
%                       On multiplie ensuite ces fonctions par Phi_eps (avec filtre PH aussi).
%                       num_idee = 8  --> Fonction de forme MsFEM, du problème purement diffusif, avec Psi_eps^2*D(x/eps) comme
%                       matrice de diffusion. On multiplie ensuite ces fonctions par Psi_eps
%                       num_idee = 9  --> Fonction de forme MsFEM, du problème purement diffusif, avec Psi_eps^2*D(x/eps) comme
%                       matrice de diffusion. On multiplie ensuite ces fonctions par Psi_eps. La différence avec l'idée 8
%                       est que l'on prend intentionnellement un grand pas de discrétisation, pour obtenir une mauvaise
%                       approximation de la fonction Psi.
%                       num_idee = 10  --> Fonction de forme MsFEM_3, du problème aux valeurs propres, avec conditions aux
%                       limites de Dirichlet. On résout le pb aux vp sur [0,2H].
%                       num_idee = 11  --> Fonction de forme MsFEM, du problème aux valeurs propres, avec Phi_eps^2*D(x/eps) comme
%                       matrice de diffusion, où Phi_eps restriction sur [H,2H] de la solution du problème stationnaire avec CL Dirichlet homogène sur [0,3H].
%                       On multiplie ensuite ces fonctions par Phi_eps.
%                       num_idee = 12  --> Fonction de forme MsFEM, du problème aux valeurs propres, avec Phi_eps^2*D(x/eps) comme
%                       matrice de diffusion, où Phi_eps restriction sur [H,2H] de la solution du problème stationnaire avec CL Neumann homogène sur [0,3H].
%                       On multiplie ensuite ces fonctions par Phi_eps.
%                       num_idee = 13  --> Fonction de forme MsFEM, du problème aux valeurs propres, avec Phi_eps^2*D(x/eps) comme
%                       matrice de diffusion, où Phi_eps restriction sur [H,2H] puis sur [2H,3H] de la solution du problème stationnaire avec CL Dirichlet homogène sur [0,4H].
%                       On multiplie ensuite ces fonctions par Phi_eps.
%                       num_idee = 14  --> Fonction de forme MsFEM, du problème aux valeurs propres, avec Phi_eps^2*D(x/eps) comme
%                       matrice de diffusion, où Phi_eps restriction sur [0,H], puis sur [H,2H] de la solution du problème stationnaire avec CL Dirichlet homogène sur [0,2H].
%                       On multiplie ensuite ces fonctions par Phi_eps.
%                       num_idee = 15  --> Fonction de forme MsFEM, du problème aux valeurs propres, avec Phi_eps^2*D(x/eps) comme
%                       matrice de diffusion, où Phi_eps restriction sur [H,2H] de la solution du problème stationnaire avec CL Dirichlet homogène sur [0,3H].
%                       On multiplie ensuite ces fonctions par Phi_eps. On calcul ensuite u_eps à l'aide de Phi_eps, et on divise la fonction de forme finale par u_eps
%                       num_idee = 16  --> Fonction de forme MsFEM, du problème purement diffusif, avec Phi_eps^2*D(x/eps) comme
%                       matrice de diffusion, où Phi_eps restriction sur [4H,5H] de la solution du problème stationnaire avec CL Dirichlet homogène sur [0,9H].
%                       On multiplie ensuite ces fonctions par Phi_eps.
%                       num_idee = 17  --> Fonction de forme MsFEM, du problème aux valeurs propres, avec Phi_eps^2*D(x/eps) comme
%                       matrice de diffusion, où Phi_eps restriction sur [H,2H] de la solution du problème stationnaire avec CL Dirichlet homogène sur [0,3H].
%                       On enlève ensuite la moyenne locale afin d'essayer d'enlever le sinus macroscopique. On multiplie ensuite ces fonctions par Phi_eps.
%                       num_idee = 18  --> Fonction de forme MsFEM, du problème aux valeurs propres, avec Phi_eps^2*D(x/eps) comme
%                       matrice de diffusion, où Phi_eps solution du problème stationnaire avec CL Dirichlet homogène.
%                       On enlève ensuite la moyenne locale, mais cette fois si en trichant, car on divise directement par un sinus, afin d'essayer d'enlever le sinus macroscopique.
%                       On multiplie ensuite ces fonctions par Phi_eps.
%                       num_idee = 19  --> Fonction de forme MsFEM, du problème aux valeurs propres, avec Phi_eps^2*D(x/eps) comme
%                       matrice de diffusion, où Phi_eps restriction sur [H,2H] de la solution du problème stationnaire avec CL periodique sur [0,3H].
%                       On multiplie ensuite ces fonctions par Phi_eps.
%                       num_idee = 20  --> Fonction de forme P1 classique, auxquelles on multiplie Psi(x/eps).
%                       num_idee = 21  --> Fonction de forme MsFEM, du problème aux valeurs propres, avec Phi_eps^2*D(x/eps) comme
%                       matrice de diffusion, où Phi_eps restriction sur [4H,5H] de la solution du problème stationnaire avec CL periodique sur [0,9H].
%                       On multiplie ensuite ces fonctions par Phi_eps.
%                       num_idee = 22  --> Fonction de forme MsFEM, du problème aux valeurs propres, avec Phi_eps^2*D(x/eps) comme
%                       matrice de diffusion, où Phi_eps restriction sur
%                       [H,2H] de la solution du problème stationnaire avec CL periodique sur [0,3H] domaine sur lequel on a
%                       modifié sigma_eps afin d'avoir la même valeur moyenne de sigma_eps sur chaque element. Pour cela
%                       on ajoute une sorte de dirac. On multiplie ensuite ces fonctions par Phi_eps.
%                       num_idee = 23  --> Fonction de forme MsFEM, du problème aux valeurs propres(ou on a ajoute un filtre), avec Phi_eps^2*D(x/eps) comme
%                       matrice de diffusion, où Phi_eps restriction sur [H,2H] de la solution du problème stationnaire avec
%                       CL du filtre voir article CLB+XB sur [0,3H]. On multiplie ensuite ces fonctions par Phi_eps.
%
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
if num_idee==0

    KK_g = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_g = sparse(Nbpt_h,1);     % vecteur second membre
    KK_d = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_d = sparse(Nbpt_h,1);     % vecteur second membre

    % boucle pour les matrices EF
    % ------------------------
    for l=2:(Nbpt_h-2)   %La première et dernière ligne doivent rester nulle

        KK_g(l,l)=(1/h)*(A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_g(l,l+1)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK_g(l+1,l)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        KK_d(l,l)=(1/h)*(A((noeud_i)*H+(l-1-0.5)*h,epsilon)+A((noeud_i)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_d(l,l+1)=(-1/h)*A((noeud_i)*H+(l-1+0.5)*h,epsilon);
        KK_d(l+1,l)=(-1/h)*A((noeud_i)*H+(l-1+0.5)*h,epsilon);

    end % for l

    KK_g(Nbpt_h-1,Nbpt_h-1)=(1/h)*(A((noeud_i-1)*H+(Nbpt_h-2-0.5)*h,epsilon)+A((noeud_i-1)*H+(Nbpt_h-2+0.5)*h,epsilon));
    KK_d(Nbpt_h-1,Nbpt_h-1)=(1/h)*(A((noeud_i)*H+(Nbpt_h-2-0.5)*h,epsilon)+A((noeud_i)*H+(Nbpt_h-2+0.5)*h,epsilon));

    % Calcul du second membre F
    % -------------------------
    for l=2:(Nbpt_h-1)   %La première et dernière ligne doivent rester nulle

        LL_g(l)=(1/H)*(A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)-A((noeud_i-1)*H+(l-1-0.5)*h,epsilon));
        LL_d(l)=(-1/H)*(A((noeud_i)*H+(l-1+0.5)*h,epsilon)-A((noeud_i)*H+(l-1-0.5)*h,epsilon));

    end % for l

    AA_g = KK_g ;
    AA_g(1,1)=1; AA_g(Nbpt_h,Nbpt_h)=1;

    AA_d = KK_d ;
    AA_d(1,1)=1; AA_d(Nbpt_h,Nbpt_h)=1;
    % inversion
    % ----------

    UU_g = AA_g\LL_g;
    UU_d = AA_d\LL_d;

    % réhaussement
    % ----------
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)+(l-1)*(h/H);
        UU_d(l)=UU_d(l)-(l-1)*(h/H)+1;

    end % for l

    %val = [ UU_g , UU_d];
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif num_idee==1

    % declarations
    % ------------
    KK_g = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_g = sparse(Nbpt_h,1);     % vecteur second membre
    KK_d = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_d = sparse(Nbpt_h,1);     % vecteur second membre
    MM_sigma = sparse(Nbpt_h,Nbpt_h); % matrice de masse

    % boucle pour les matrices EF
    % ------------------------
    for l=2:(Nbpt_h-2)   %La première et dernière ligne doivent rester nulle

        MM_sigma(l,l)=Sigma((noeud_i-1)*H+(l-1)*h,epsilon)*h*(2/3);
        MM_sigma(l,l+1)=Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon)*h*(1/6);
        MM_sigma(l+1,l)=Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon)*h*(1/6);

        KK_g(l,l)=(1/h)*(A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_g(l,l+1)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK_g(l+1,l)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        KK_d(l,l)=(1/h)*(A((noeud_i)*H+(l-1-0.5)*h,epsilon)+A((noeud_i)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_d(l,l+1)=(-1/h)*A((noeud_i)*H+(l-1+0.5)*h,epsilon);
        KK_d(l+1,l)=(-1/h)*A((noeud_i)*H+(l-1+0.5)*h,epsilon);

    end % for l

    KK_g(Nbpt_h-1,Nbpt_h-1)=(1/h)*(A((noeud_i-1)*H+(Nbpt_h-2-0.5)*h,epsilon)+A((noeud_i-1)*H+(Nbpt_h-2+0.5)*h,epsilon));
    KK_d(Nbpt_h-1,Nbpt_h-1)=(1/h)*(A((noeud_i)*H+(Nbpt_h-2-0.5)*h,epsilon)+A((noeud_i)*H+(Nbpt_h-2+0.5)*h,epsilon));
    MM_sigma(Nbpt_h-1,Nbpt_h-1)=Sigma((noeud_i-1)*H+(Nbpt_h-2)*h,epsilon)*h*(2/3);

    % Calcul du second membre
    % -------------------------
    for l=2:(Nbpt_h-1)   %La première et dernière ligne doivent rester nulle

        LL_g(l)=(1/(2*H))*(A((noeud_i-1)*H+l*h,epsilon)-A((noeud_i-1)*H+(l-2)*h,epsilon)) - 1/(epsilon^2)*Sigma((noeud_i-1)*H+(l-1)*h,epsilon)*h*(l/Nbpt_h) ;

        LL_d(l)=(-1/(2*H))*(A((noeud_i)*H+l*h,epsilon)-A((noeud_i)*H+(l-2)*h,epsilon)) - 1/(epsilon^2)*Sigma((noeud_i-1)*H+(l-1)*h,epsilon)*h*((Nbpt_h-l)/Nbpt_h);
    end % for l


    AA_g = KK_g + 1/(epsilon^2)*MM_sigma ;
    AA_g(1,1)=1; AA_g(Nbpt_h,Nbpt_h)=1;

    AA_d = KK_d + 1/(epsilon^2)*MM_sigma;
    AA_d(1,1)=1; AA_d(Nbpt_h,Nbpt_h)=1;
    % inversion
    % ----------

    UU_g = AA_g\LL_g;
    UU_d = AA_d\LL_d;

    % réhaussement
    % ----------
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)+(l-1)*(h/H);
        UU_d(l)=UU_d(l)-(l-1)*(h/H)+1;

    end % for l

    %val = [ UU_g ; UU_d(2:Nbpt_h)];
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif num_idee==2

    [~,vect_propre_eps]=fonction_propre_eps(Nbpt_H,Nbpt_h,epsilon,noeud_i);
    [~,vect_propre_droite]=fonction_propre_eps(Nbpt_H,Nbpt_h,epsilon,noeud_i+1);

    Phi_eps_gauche = zeros(Nbpt_h-1,1);
    Phi_eps_droite = zeros(Nbpt_h-1,1);
    for i=1:(Nbpt_h)
        xi = (i-1)*h ;
        Phi_eps_gauche(i) = vect_propre_eps(i) - (vect_propre_eps(1) + (xi/H)*(vect_propre_eps(Nbpt_h) - vect_propre_eps(1) ) );
        Phi_eps_droite(i) = vect_propre_droite(i) - (vect_propre_droite(1) + (xi/H)*(vect_propre_droite(Nbpt_h) - vect_propre_droite(1) ) );

    end
    Phi_eps_droite = Phi_eps_droite - mean(Phi_eps_droite);
    Phi_eps_gauche = Phi_eps_gauche - mean(Phi_eps_gauche);
    Y_Phi_droite=fft(Phi_eps_droite);
    Y_Phi_gauche=fft(Phi_eps_gauche);

    %Frequence de coupure, potentiellement à ajuster
    fc = 1/(1*epsilon);

    Y_Phi_droite(1:round(fc*H)) = 0;
    Y_Phi_droite(end-round(fc*H)+1:end) = 0;
    Y_Phi_gauche(1:round(fc*H)) = 0;
    Y_Phi_gauche(end-round(fc*H)+1:end) = 0;

    Phi_eps_filtre_gauche = real(ifft(Y_Phi_gauche));
    Phi_eps_filtre_gauche = Phi_eps_filtre_gauche +1;
    Phi_eps_filtre_droite = real(ifft(Y_Phi_droite));
    Phi_eps_filtre_droite = Phi_eps_filtre_droite +1;
    UU_g = zeros(Nbpt_h,1);
    UU_d = zeros(Nbpt_h,1);

    % boucle pour les matrices EF
    % ------------------------
    for i=1:(Nbpt_h)

        UU_g(i) = Phi_eps_filtre_gauche(i)- Phi_eps_filtre_gauche(1) + (1-(Phi_eps_filtre_gauche(Nbpt_h)- Phi_eps_filtre_gauche(1)))*(i-1)*(h/H);
        UU_d(i) = Phi_eps_filtre_droite(i)- Phi_eps_filtre_droite(Nbpt_h) + (1-(Phi_eps_filtre_droite(1)- Phi_eps_filtre_droite(Nbpt_h)))*(1-((i-1)*(h/H)));

    end % for l

    %val = [ UU_g ; UU_d(2:Nbpt_h)];
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif num_idee==3

    [~,vect_propre_eps]=fonction_propre_eps(Nbpt_H,Nbpt_h,epsilon,noeud_i);
    [~,vect_propre_droite]=fonction_propre_eps(Nbpt_H,Nbpt_h,epsilon,noeud_i+1);

    Phi_eps_gauche = zeros(Nbpt_h-1,1);
    Phi_eps_droite = zeros(Nbpt_h-1,1);
    for i=1:(Nbpt_h)
        xi = (i-1)*h ;
        Phi_eps_gauche(i) = vect_propre_eps(i) - (vect_propre_eps(1) + (xi/H)*(vect_propre_eps(Nbpt_h) - vect_propre_eps(1) ) );
        Phi_eps_droite(i) = vect_propre_droite(i) - (vect_propre_droite(1) + (xi/H)*(vect_propre_droite(Nbpt_h) - vect_propre_droite(1) ) );

    end
    Phi_eps_droite = Phi_eps_droite - mean(Phi_eps_droite);
    Phi_eps_gauche = Phi_eps_gauche - mean(Phi_eps_gauche);
    Y_Phi_droite=fft(Phi_eps_droite);
    Y_Phi_gauche=fft(Phi_eps_gauche);

    %Frequence de coupure, potentiellement à ajuster
    fc = 1/(1*epsilon);

    Y_Phi_droite(1:round(fc*H)) = 0;
    Y_Phi_droite(end-round(fc*H)+1:end) = 0;
    Y_Phi_gauche(1:round(fc*H)) = 0;
    Y_Phi_gauche(end-round(fc*H)+1:end) = 0;

    Phi_eps_filtre_gauche = real(ifft(Y_Phi_gauche));
    Phi_eps_filtre_gauche = Phi_eps_filtre_gauche +1;
    Phi_eps_filtre_droite = real(ifft(Y_Phi_droite));
    Phi_eps_filtre_droite = Phi_eps_filtre_droite +1;
    UU_g = zeros(Nbpt_h,1);
    UU_d = zeros(Nbpt_h,1);

    % boucle pour les matrices EF
    % ------------------------
    for i=1:(Nbpt_h)

        UU_g(i) = Phi_eps_filtre_gauche(i)*(i-1)*(h/H);
        UU_d(i) = Phi_eps_filtre_droite(i)*(1-((i-1)*(h/H)));

    end % for l
    Kg = Phi_eps_filtre_gauche(Nbpt_h);
    UU_g = UU_g./ Kg ;
    Kd = Phi_eps_filtre_droite(1);
    UU_d = UU_d./Kd ;

    %val = [ UU_g ; UU_d(2:Nbpt_h)];
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif num_idee==4

    [~,vect_propre_eps]=fonction_propre_eps(Nbpt_H,Nbpt_h,epsilon,noeud_i);
    [~,vect_propre_droite]=fonction_propre_eps(Nbpt_H,Nbpt_h,epsilon,noeud_i+1);


    UU_g = zeros(Nbpt_h,1);
    UU_d = zeros(Nbpt_h,1);

    % boucle pour les matrices EF
    % ------------------------
    for i=1:(Nbpt_h)

        UU_g(i) = vect_propre_eps(i)*(i-1)*(h/H);
        UU_d(i) = vect_propre_droite(i)*(1-((i-1)*(h/H)));

    end % for l
    Kg = vect_propre_eps(Nbpt_h);
    UU_g = UU_g./ Kg ;
    Kd = vect_propre_droite(1);
    UU_d = UU_d./Kd ;

    %val = [ UU_g ; UU_d(2:Nbpt_h)];
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif num_idee==5


    [~,vect_propre_eps]=fonction_propre_eps(Nbpt_H,Nbpt_h,epsilon,noeud_i);
    [~,vect_propre_droite]=fonction_propre_eps(Nbpt_H,Nbpt_h,epsilon,noeud_i+1);

    UU_g = zeros(Nbpt_h,1);
    UU_d = zeros(Nbpt_h,1);

    % boucle pour les matrices EF
    % ------------------------
    for i=1:(Nbpt_h)
        UU_g(i) = vect_propre_eps(i)*(i-1)*(h/H);
        UU_d(i) = vect_propre_droite(i)*(1-((i-1)*(h/H)));
    end % for l

    %val = [ UU_g ; UU_d(2:Nbpt_h)];
    %val = val*(sqrt(H)); %On "renormalise" pour que le vecteur soit toujours aux alentours de 1
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif num_idee==6

    UU_g = zeros(Nbpt_h,1);
    UU_d = zeros(Nbpt_h,1);

    % boucle pour les matrices EF
    % ------------------------
    for i=1:(Nbpt_h)
        UU_g(i) = (i-1)*(h/H);
        UU_d(i) = (1-((i-1)*(h/H)));
    end % for l

    %val = [ UU_g ; UU_d(2:Nbpt_h)];
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif num_idee==7

    [~,vect_propre_eps]=fonction_propre_eps(Nbpt_H,Nbpt_h,epsilon,noeud_i);
    [~,vect_propre_droite]=fonction_propre_eps(Nbpt_H,Nbpt_h,epsilon,noeud_i+1);

    Phi_eps_gauche = zeros(Nbpt_h-1,1);
    Phi_eps_droite = zeros(Nbpt_h-1,1);
    for i=1:(Nbpt_h)
        xi = (i-1)*h ;
        Phi_eps_gauche(i) = vect_propre_eps(i) - (vect_propre_eps(1) + (xi/H)*(vect_propre_eps(Nbpt_h) - vect_propre_eps(1) ) );
        Phi_eps_droite(i) = vect_propre_droite(i) - (vect_propre_droite(1) + (xi/H)*(vect_propre_droite(Nbpt_h) - vect_propre_droite(1) ) );

    end
    Phi_eps_droite = Phi_eps_droite - mean(Phi_eps_droite);
    Phi_eps_gauche = Phi_eps_gauche - mean(Phi_eps_gauche);
    Y_Phi_droite=fft(Phi_eps_droite);
    Y_Phi_gauche=fft(Phi_eps_gauche);

    %Frequence de coupure, potentiellement à ajuster
    fc = 1/(1*epsilon);

    Y_Phi_droite(1:round(fc*H)) = 0;
    Y_Phi_droite(end-round(fc*H)+1:end) = 0;
    Y_Phi_gauche(1:round(fc*H)) = 0;
    Y_Phi_gauche(end-round(fc*H)+1:end) = 0;

    Phi_eps_filtre_gauche = real(ifft(Y_Phi_gauche));
    Phi_eps_filtre_gauche = Phi_eps_filtre_gauche +1;
    Phi_eps_filtre_droite = real(ifft(Y_Phi_droite));
    Phi_eps_filtre_droite = Phi_eps_filtre_droite +1;

    KK_g = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_g = zeros(Nbpt_h,1);     % vecteur second membre
    KK_d = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_d = zeros(Nbpt_h,1);     % vecteur second membre

    % boucle pour les matrices EF
    % ------------------------
    for l=2:(Nbpt_h-2)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(Phi_eps_filtre_gauche(l)+Phi_eps_filtre_gauche(l-1));
        phi_apres_g = 0.5*(Phi_eps_filtre_gauche(l)+Phi_eps_filtre_gauche(l+1));

        KK_g(l,l)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_g(l,l+1)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK_g(l+1,l)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        phi_avant_d = 0.5*(Phi_eps_filtre_droite(l)+Phi_eps_filtre_droite(l-1));
        phi_apres_d = 0.5*(Phi_eps_filtre_droite(l)+Phi_eps_filtre_droite(l+1));

        KK_d(l,l)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_d(l,l+1)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);
        KK_d(l+1,l)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);

    end % for l
    phi_avant_d = 0.5*(Phi_eps_filtre_droite(Nbpt_h-1)+Phi_eps_filtre_droite(Nbpt_h-1-1));
    phi_apres_d = 0.5*(Phi_eps_filtre_droite(Nbpt_h-1)+Phi_eps_filtre_droite(Nbpt_h-1+1));
    phi_avant_g = 0.5*(Phi_eps_filtre_gauche(Nbpt_h-1)+Phi_eps_filtre_gauche(Nbpt_h-1-1));
    phi_apres_g = 0.5*(Phi_eps_filtre_gauche(Nbpt_h-1)+Phi_eps_filtre_gauche(Nbpt_h-1+1));

    KK_g(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
    KK_d(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)

    % Calcul du second membre F
    % -------------------------
    for l=2:(Nbpt_h-1)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(Phi_eps_filtre_gauche(l)+Phi_eps_filtre_gauche(l-1));
        phi_apres_g = 0.5*(Phi_eps_filtre_gauche(l)+Phi_eps_filtre_gauche(l+1));
        phi_avant_d = 0.5*(Phi_eps_filtre_droite(l)+Phi_eps_filtre_droite(l-1));
        phi_apres_d = 0.5*(Phi_eps_filtre_droite(l)+Phi_eps_filtre_droite(l+1));

        LL_g(l)=(-1/H)*(phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)-phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon));

        LL_d(l)=(1/H)*(phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)-phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon));

    end % for l

    AA_g = KK_g ;
    AA_g(1,1)=1; AA_g(Nbpt_h,Nbpt_h)=1;

    AA_d = KK_d ;
    AA_d(1,1)=1; AA_d(Nbpt_h,Nbpt_h)=1;
    % inversion
    % ----------

    UU_g = AA_g\LL_g;
    UU_d = AA_d\LL_d;

    % réhaussement
    % ----------
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)+(l-1)*(h/H);
        UU_d(l)=UU_d(l)-(l-1)*(h/H)+1;

    end % for l

    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)*Phi_eps_filtre_gauche(l);
        UU_d(l)=UU_d(l)*Phi_eps_filtre_droite(l);

    end % for l


    %val = [ UU_g ; UU_d(2:Nbpt_h)];
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif num_idee==8

    [~,Psi] = Solution_pb_spectral(Nbpt_h);
    vect_propre_eps = zeros(Nbpt_h,1);
    vect_propre_droite = zeros(Nbpt_h,1);
    for i=1:Nbpt_h
        xi_g = (noeud_i-1)*H+(i-1)*h ;
        xi_d = (noeud_i)*H+(i-1)*h ;
        yi_g = rem(xi_g/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk_g = fix(yi_g*(Nbpt_h-1)) +1; %partie entière de y*(Nbpt-1) (indice pour les matrices)
        yi_d = rem(xi_d/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk_d = fix(yi_d*(Nbpt_h-1)) +1;

        %         vect_propre_eps(i) = Psi(yk_g);
        %         vect_propre_droite(i) = Psi(yk_d);

        t_g = yi_g*(Nbpt_h-1) - fix(yi_g*(Nbpt_h-1)) ;
        vect_propre_eps(i) = Psi(yk_g)*(1-t_g) +Psi(yk_g+1)*t_g;
        t_d = yi_d*(Nbpt_h-1) - fix(yi_d*(Nbpt_h-1)) ;
        vect_propre_droite(i) = Psi(yk_d)*(1-t_d) +Psi(yk_d+1)*t_d;
    end

    KK_g = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_g = zeros(Nbpt_h,1);     % vecteur second membre
    KK_d = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_d = zeros(Nbpt_h,1);     % vecteur second membre

    % boucle pour les matrices EF
    % ------------------------
    for l=2:(Nbpt_h-2)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));

        KK_g(l,l)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_g(l,l+1)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK_g(l+1,l)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        KK_d(l,l)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_d(l,l+1)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);
        KK_d(l+1,l)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);

    end % for l
    phi_avant_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1-1));
    phi_apres_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1+1));
    phi_avant_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1-1));
    phi_apres_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1+1));

    KK_g(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
    KK_d(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)

    % Calcul du second membre F
    % -------------------------
    for l=2:(Nbpt_h-1)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));
        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        LL_g(l)=(-1/H)*(-1*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)+phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon));

        LL_d(l)=(1/H)*(-phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)+phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon));

    end % for l

    AA_g = KK_g ;
    AA_g(1,1)=1; AA_g(Nbpt_h,Nbpt_h)=1;

    AA_d = KK_d ;
    AA_d(1,1)=1; AA_d(Nbpt_h,Nbpt_h)=1;
    % inversion
    % ----------

    UU_g = AA_g\LL_g;
    UU_d = AA_d\LL_d;

    % réhaussement
    % ----------
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)+(l-1)*(h/H);
        UU_d(l)=UU_d(l)-(l-1)*(h/H)+1;

    end % for l

    %Multiplication par Psi_eps
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)*vect_propre_eps(l);
        UU_d(l)=UU_d(l)*vect_propre_droite(l);

    end % for l
    %val = [ UU_g , UU_d];
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif num_idee==9

    Nbpt_pour_mauvaise_approximation = 11 ; %Valeur à modifier en fonction de l'approximation plus ou moins mauvaise que l'on veut

    [~,Psi_mauvaise_approximation] = Solution_pb_spectral(Nbpt_pour_mauvaise_approximation);
    Psi = decompose(Psi_mauvaise_approximation,Nbpt_h); %Interpollation de sur le maillage fin de l'approximation grossière de Psi

    vect_propre_eps = zeros(Nbpt_h,1);
    vect_propre_droite = zeros(Nbpt_h,1);
    for i=1:Nbpt_h
        xi_g = (noeud_i-1)*H+(i-1)*h ;
        xi_d = (noeud_i)*H+(i-1)*h ;
        yi_g = rem(xi_g/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk_g = fix(yi_g*(Nbpt_h-1)) +1; %partie entière de y*(Nbpt-1) (indice pour les matrices)

        t_g = yi_g*(Nbpt_h-1) - fix(yi_g*(Nbpt_h-1)) ;
        vect_propre_eps(i) = Psi(yk_g)*(1-t_g) +Psi(yk_g+1)*t_g;


        yi_d = rem(xi_d/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk_d = fix(yi_d*(Nbpt_h-1)) +1;

        t_d = yi_d*(Nbpt_h-1) - fix(yi_d*(Nbpt_h-1)) ;

        vect_propre_droite(i) = Psi(yk_d)*(1-t_d) +Psi(yk_d+1)*t_d;
    end

    KK_g = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_g = zeros(Nbpt_h,1);     % vecteur second membre
    KK_d = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_d = zeros(Nbpt_h,1);     % vecteur second membre

    % boucle pour les matrices EF
    % ------------------------
    for l=2:(Nbpt_h-2)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));

        KK_g(l,l)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_g(l,l+1)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK_g(l+1,l)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        KK_d(l,l)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_d(l,l+1)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);
        KK_d(l+1,l)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);

    end % for l
    phi_avant_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1-1));
    phi_apres_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1+1));
    phi_avant_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1-1));
    phi_apres_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1+1));

    KK_g(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
    KK_d(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)

    % Calcul du second membre F
    % -------------------------
    for l=2:(Nbpt_h-1)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));
        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        LL_g(l)=(-1/H)*(-1*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)+phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon));

        LL_d(l)=(1/H)*(-phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)+phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon));

    end % for l

    AA_g = KK_g ;
    AA_g(1,1)=1; AA_g(Nbpt_h,Nbpt_h)=1;

    AA_d = KK_d ;
    AA_d(1,1)=1; AA_d(Nbpt_h,Nbpt_h)=1;
    % inversion
    % ----------

    UU_g = AA_g\LL_g;
    UU_d = AA_d\LL_d;

    % réhaussement
    % ----------
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)+(l-1)*(h/H);
        UU_d(l)=UU_d(l)-(l-1)*(h/H)+1;

    end % for l

    %Multiplication par Psi_eps
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)*vect_propre_eps(l);
        UU_d(l)=UU_d(l)*vect_propre_droite(l);

    end % for l
    %val = [ UU_g ; UU_d(2:Nbpt_h)];

    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif num_idee==10

    % declarations
    % ------------
    MM = matM_ref(2*Nbpt_h-1); % matrice de masse
    KK = sparse(2*Nbpt_h-1,2*Nbpt_h-1); % matrice de rigidite
    MM_sigma = sparse(2*Nbpt_h-1,2*Nbpt_h-1); % matrice de masse

    for l=2:(2*Nbpt_h-1-2)   %La première et dernière ligne doivent rester nulle

        KK(l,l)=(1/h)*(A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK(l+1,l)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        MM_sigma(l,l)=(h/3)*(Sigma((noeud_i-1)*H+(l-1-0.5)*h,epsilon)+Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon));
        MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

    end % for l

    MM_sigma(2*Nbpt_h-1-1,2*Nbpt_h-1-1)=(h/3)*(Sigma((noeud_i-1)*H+(2*Nbpt_h-1-2-0.5)*h,epsilon)+Sigma((noeud_i-1)*H+(2*Nbpt_h-1-2+0.5)*h,epsilon));
    KK(2*Nbpt_h-1-1,2*Nbpt_h-1-1)=(1/h)*(A((noeud_i-1)*H+(2*Nbpt_h-1-2-0.5)*h,epsilon)+A((noeud_i-1)*H+(2*Nbpt_h-1-2+0.5)*h,epsilon));

    AA = MM_sigma + KK*(epsilon^2) ;
    AA(1,1)=1; AA(2*Nbpt_h-1,2*Nbpt_h-1)=1;

    %[EV,DV] = eigs(AA,MM,5,'smallestabs');
    [EV,DV] = eigs(AA,MM,1,'smallestabs');


    [~,arg_min] = min(diag(DV));
    Vecteur_propre = EV(:,arg_min);
    Vecteur_propre = Vecteur_propre*(sqrt(H)); %On "renormalise" pour que le vecteur soit toujours aux alentours de 1

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre(fix(Nbpt_h/2))<0)
        Vecteur_propre = -Vecteur_propre;
    end

    %Vecteur_propre = Vecteur_propre/Vecteur_propre(Nbpt_h);
    % réhaussement
    % ----------
    %     for l=1:Nbpt_h
    %
    %         UU_g(l)=Vecteur_propre(l)+(1-Vecteur_propre(Nbpt_h))*(l-1)*(h/H);
    %         UU_d(l)=Vecteur_propre(Nbpt_h-1+l)+(1-Vecteur_propre(Nbpt_h))*(1-((l-1)*(h/H)));
    %
    %     end % for l

    %val = Vecteur_propre;
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif num_idee==11

    % declarations
    % ------------
    MM = matM_ref(3*(Nbpt_h-1)+1); % matrice de masse
    KK = sparse(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1); % matrice de rigidite
    MM_sigma = sparse(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1); % matrice de masse

    for l=2:(3*(Nbpt_h-1)+1-2)   %La première et dernière ligne doivent rester nulle

        KK(l,l)=(1/h)*(A((noeud_i-2)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-2)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-2)*H+(l-1+0.5)*h,epsilon);
        KK(l+1,l)=(-1/h)*A((noeud_i-2)*H+(l-1+0.5)*h,epsilon);

        MM_sigma(l,l)=(h/3)*(Sigma((noeud_i-2)*H+(l-1-0.5)*h,epsilon)+Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon));
        MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon);
        MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon);

    end % for l

    MM_sigma(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=(h/3)*(Sigma((noeud_i-2)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+Sigma((noeud_i-2)*H+(3*(Nbpt_h-1)+1-2+0.5)*h,epsilon));
    KK(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=(1/h)*(A((noeud_i-2)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+A((noeud_i-2)*H+(3*(Nbpt_h-1)+1-2+0.5)*h,epsilon));

    AA = MM_sigma + KK*(epsilon^2) ;
    AA(1,1)=1; AA(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1)=1;

    %[EV,DV] = eigs(AA,MM,5,'smallestabs');
    [EV,DV] = eigs(AA,MM,1,'smallestabs');


    [~,arg_min] = min(diag(DV));
    Vecteur_propre = EV(:,arg_min);
    Vecteur_propre = Vecteur_propre*(sqrt(H)); %On "renormalise" pour que le vecteur soit toujours aux alentours de 1

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre(fix(Nbpt_h/2))<0)
        Vecteur_propre = -Vecteur_propre;
    end

    vect_propre_eps = Vecteur_propre(Nbpt_h:2*Nbpt_h-1);
    % -------------------------------------------------------------------------

    % declarations
    % ------------
    MM = matM_ref(3*(Nbpt_h-1)+1); % matrice de masse
    KK = sparse(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1); % matrice de rigidite
    MM_sigma = sparse(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1); % matrice de masse

    for l=2:(3*(Nbpt_h-1)+1-2)   %La première et dernière ligne doivent rester nulle

        KK(l,l)=(1/h)*(A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK(l+1,l)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        MM_sigma(l,l)=(h/3)*(Sigma((noeud_i-1)*H+(l-1-0.5)*h,epsilon)+Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon));
        MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

    end % for l

    MM_sigma(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=(h/3)*(Sigma((noeud_i-1)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+Sigma((noeud_i-1)*H+(3*(Nbpt_h-1)+1-2+0.5)*h,epsilon));
    KK(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=(1/h)*(A((noeud_i-1)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+A((noeud_i-1)*H+(3*(Nbpt_h-1)+1-2+0.5)*h,epsilon));

    AA = MM_sigma + KK*(epsilon^2) ;
    AA(1,1)=1; AA(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1)=1;

    %[EV,DV] = eigs(AA,MM,5,'smallestabs');
    [EV,DV] = eigs(AA,MM,1,'smallestabs');


    [~,arg_min] = min(diag(DV));
    Vecteur_propre = EV(:,arg_min);
    Vecteur_propre = Vecteur_propre*(sqrt(H)); %On "renormalise" pour que le vecteur soit toujours aux alentours de 1

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre(fix(Nbpt_h/2))<0)
        Vecteur_propre = -Vecteur_propre;
    end

    vect_propre_droite = Vecteur_propre(Nbpt_h:2*Nbpt_h-1);
    % -------------------------------------------------------------------------

    KK_g = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_g = zeros(Nbpt_h,1);     % vecteur second membre
    KK_d = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_d = zeros(Nbpt_h,1);     % vecteur second membre

    % boucle pour les matrices EF
    % ------------------------
    for l=2:(Nbpt_h-2)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));

        KK_g(l,l)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_g(l,l+1)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK_g(l+1,l)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        KK_d(l,l)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_d(l,l+1)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);
        KK_d(l+1,l)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);

    end % for l
    phi_avant_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1-1));
    phi_apres_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1+1));
    phi_avant_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1-1));
    phi_apres_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1+1));

    KK_g(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
    KK_d(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)

    % Calcul du second membre F
    % -------------------------
    for l=2:(Nbpt_h-1)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));
        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        LL_g(l)=(-1/H)*(-1*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)+phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon));

        LL_d(l)=(1/H)*(-phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)+phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon));

    end % for l

    AA_g = KK_g ;
    AA_g(1,1)=1; AA_g(Nbpt_h,Nbpt_h)=1;

    AA_d = KK_d ;
    AA_d(1,1)=1; AA_d(Nbpt_h,Nbpt_h)=1;
    % inversion
    % ----------

    UU_g = AA_g\LL_g;
    UU_d = AA_d\LL_d;

    % réhaussement
    % ----------
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)+(l-1)*(h/H);
        UU_d(l)=UU_d(l)-(l-1)*(h/H)+1;

    end % for l

    %Multiplication par Psi_eps
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)*vect_propre_eps(l)/vect_propre_eps(Nbpt_h);
        UU_d(l)=UU_d(l)*vect_propre_droite(l)/vect_propre_droite(1);
        %         UU_g(l)=UU_g(l)*vect_propre_eps(l);
        %         UU_d(l)=UU_d(l)*vect_propre_droite(l);

    end % for l
    %val = [ UU_g ; UU_d(2:Nbpt_h)];
    %     coeff = UU_g(Nbpt_h);
    %     val = val./coeff ;
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif num_idee==12

    % declarations
    % ------------
    MM = sparse(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1); % matrice de masse
    KK = sparse(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1); % matrice de rigidite
    MM_sigma = sparse(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1); % matrice de masse
    
    % boucle pour les matrices EF
    % ------------------------
    for l=2:(3*(Nbpt_h-1)+1-2)   %La première et dernière ligne doivent rester nulle
        
        MM(l,l)=h*(2/3);
        MM(l,l+1)=h*(1/6);
        MM(l+1,l)=h*(1/6);
        
        MM_sigma(l,l)=(h/3)*(Sigma((noeud_i-2)*H+(l-1-0.5)*h,epsilon)+Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon));
        MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon);
        MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon);
    
        KK(l,l)=(1/h)*(A((noeud_i-2)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-2)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-2)*H+(l-1+0.5)*h,epsilon);
        KK(l+1,l)=(-1/h)*A((noeud_i-2)*H+(l-1+0.5)*h,epsilon);
    
    end % for l
    
    MM(1,1)=h*(1/3);
    MM(1,2)=h*(1/6);
    MM(2,1)=h*(1/6);
    MM(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1)=h*(1/3);
    MM(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1)=h*(1/6);
    MM(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1-1)=h*(1/6);
    MM(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=h*(2/3);
    
    KK(1,1)=(1/h)*A((noeud_i-2)*H+(0.5)*h,epsilon);
    KK(1,2)=(-1/h)*A((noeud_i-2)*H+(0.5)*h,epsilon);
    KK(2,1)=(-1/h)*A((noeud_i-2)*H+(0.5)*h,epsilon);
    KK(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1)=(1/h)*A((noeud_i-2)*H+(3*(Nbpt_h-1)+1-1-0.5)*h,epsilon);
    KK(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1)=(-1/h)*A((noeud_i-2)*H+(3*(Nbpt_h-1)+1-1-1+0.5)*h,epsilon);
    KK(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1-1)=(-1/h)*A((noeud_i-2)*H+(3*(Nbpt_h-1)+1-1-1+0.5)*h,epsilon);
    KK(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=(1/h)*(A((noeud_i-2)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+A((noeud_i-2)*H+(3*(Nbpt_h-1)+1-2+0.5)*h,epsilon));
    
    MM_sigma(1,1)=(h/3)*Sigma((noeud_i-2)*H+(0.5)*h,epsilon);
    MM_sigma(1,2)=h*(1/6)*Sigma((noeud_i-2)*H+(0.5)*h,epsilon);
    MM_sigma(2,1)=h*(1/6)*Sigma((noeud_i-2)*H+(0.5)*h,epsilon);
    MM_sigma(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1)=(h/3)*Sigma((noeud_i-2)*H+(Nbpt_h-1-0.5)*h,epsilon);
    MM_sigma(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1)=h*(1/6)*Sigma((noeud_i-2)*H+(3*(Nbpt_h-1)+1-1-1+0.5)*h,epsilon);
    MM_sigma(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1-1)=h*(1/6)*Sigma((noeud_i-2)*H+(3*(Nbpt_h-1)+1-1-1+0.5)*h,epsilon);
    MM_sigma(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=(h/3)*(Sigma((noeud_i-2)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+Sigma((noeud_i-2)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon));
    
    
    AA = MM_sigma + KK*(epsilon*epsilon) ;
    
    %[EV,DV] = eigs(AA,MM,5,'smallestabs');
    [EV,DV] = eigs(AA,MM,1,'smallestabs');
    
    
    [~,arg_min] = min(diag(DV));
    Vecteur_propre = EV(:,arg_min);
    Vecteur_propre = Vecteur_propre*(sqrt(3*H)); %On "renormalise" pour que le vecteur soit toujours aux alentours de 1
    
    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre(fix(Nbpt_h/2))<0)
        Vecteur_propre = -Vecteur_propre;
    end
    
    vect_propre_eps_gauche = Vecteur_propre(Nbpt_h:2*Nbpt_h-1);
 % ------------------------ % ------------------------ % ------------------------ % ------------------------

    KK = sparse(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1); % matrice de rigidite
    MM_sigma = sparse(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1); % matrice de masse
    
    % boucle pour les matrices EF
    % ------------------------
    for l=2:(3*(Nbpt_h-1)+1-2)   %La première et dernière ligne doivent rester nulle
        
        MM_sigma(l,l)=(h/3)*(Sigma((noeud_i-1)*H+(l-1-0.5)*h,epsilon)+Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon));
        MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
    
        KK(l,l)=(1/h)*(A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK(l+1,l)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
    
    end % for l
    
    MM(1,1)=h*(1/3);
    MM(1,2)=h*(1/6);
    MM(2,1)=h*(1/6);
    MM(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1)=h*(1/3);
    MM(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1)=h*(1/6);
    MM(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1-1)=h*(1/6);
    MM(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=h*(2/3);
    
    KK(1,1)=(1/h)*A((noeud_i-1)*H+(0.5)*h,epsilon);
    KK(1,2)=(-1/h)*A((noeud_i-1)*H+(0.5)*h,epsilon);
    KK(2,1)=(-1/h)*A((noeud_i-1)*H+(0.5)*h,epsilon);
    KK(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1)=(1/h)*A((noeud_i-1)*H+(3*(Nbpt_h-1)+1-1-0.5)*h,epsilon);
    KK(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1)=(-1/h)*A((noeud_i-1)*H+(3*(Nbpt_h-1)+1-1-1+0.5)*h,epsilon);
    KK(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1-1)=(-1/h)*A((noeud_i-1)*H+(3*(Nbpt_h-1)+1-1-1+0.5)*h,epsilon);
    KK(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=(1/h)*(A((noeud_i-1)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+A((noeud_i-1)*H+(3*(Nbpt_h-1)+1-2+0.5)*h,epsilon));
    
    MM_sigma(1,1)=(h/3)*Sigma((noeud_i-1)*H+(0.5)*h,epsilon);
    MM_sigma(1,2)=h*(1/6)*Sigma((noeud_i-1)*H+(0.5)*h,epsilon);
    MM_sigma(2,1)=h*(1/6)*Sigma((noeud_i-1)*H+(0.5)*h,epsilon);
    MM_sigma(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1)=(h/3)*Sigma((noeud_i-1)*H+(3*(Nbpt_h-1)+1-1-0.5)*h,epsilon);
    MM_sigma(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1)=h*(1/6)*Sigma((noeud_i-1)*H+(3*(Nbpt_h-1)+1-1-1+0.5)*h,epsilon);
    MM_sigma(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1-1)=h*(1/6)*Sigma((noeud_i-1)*H+(3*(Nbpt_h-1)+1-1-1+0.5)*h,epsilon);
    MM_sigma(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=(h/3)*(Sigma((noeud_i-1)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+Sigma((noeud_i-1)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon));
    
    
    AA = MM_sigma + KK*(epsilon*epsilon) ;
    
    %[EV,DV] = eigs(AA,MM,5,'smallestabs');
    [EV,DV] = eigs(AA,MM,1,'smallestabs');
    
    
    [~,arg_min] = min(diag(DV));
    Vecteur_propre = EV(:,arg_min);
    Vecteur_propre = Vecteur_propre*(sqrt(3*H)); %On "renormalise" pour que le vecteur soit toujours aux alentours de 1
    
    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre(fix(Nbpt_h/2))<0)
        Vecteur_propre = -Vecteur_propre;
    end
    
    vect_propre_droite = Vecteur_propre(Nbpt_h:2*Nbpt_h-1);
 % ------------------------ % ------------------------ % ------------------------ % ------------------------

    KK_g = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_g = zeros(Nbpt_h,1);     % vecteur second membre
    KK_d = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_d = zeros(Nbpt_h,1);     % vecteur second membre

    % boucle pour les matrices EF
    % ------------------------
    for l=2:(Nbpt_h-2)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps_gauche(l)+vect_propre_eps_gauche(l-1));
        phi_apres_g = 0.5*(vect_propre_eps_gauche(l)+vect_propre_eps_gauche(l+1));

        KK_g(l,l)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_g(l,l+1)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK_g(l+1,l)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        KK_d(l,l)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_d(l,l+1)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);
        KK_d(l+1,l)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);

    end % for l
    phi_avant_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1-1));
    phi_apres_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1+1));
    phi_avant_g = 0.5*(vect_propre_eps_gauche(Nbpt_h-1)+vect_propre_eps_gauche(Nbpt_h-1-1));
    phi_apres_g = 0.5*(vect_propre_eps_gauche(Nbpt_h-1)+vect_propre_eps_gauche(Nbpt_h-1+1));

    KK_g(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
    KK_d(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)

    % Calcul du second membre F
    % -------------------------
    for l=2:(Nbpt_h-1)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps_gauche(l)+vect_propre_eps_gauche(l-1));
        phi_apres_g = 0.5*(vect_propre_eps_gauche(l)+vect_propre_eps_gauche(l+1));
        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        LL_g(l)=(-1/H)*(phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)-phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon));

        LL_d(l)=(1/H)*(phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)-phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon));

    end % for l

    AA_g = KK_g ;
    AA_g(1,1)=1; AA_g(Nbpt_h,Nbpt_h)=1;

    AA_d = KK_d ;
    AA_d(1,1)=1; AA_d(Nbpt_h,Nbpt_h)=1;
    % inversion
    % ----------

    UU_g = AA_g\LL_g;
    UU_d = AA_d\LL_d;

    % réhaussement
    % ----------
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)+(l-1)*(h/H);
        UU_d(l)=UU_d(l)-(l-1)*(h/H)+1;

    end % for l

    % Multiplication par Psi_eps
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)*vect_propre_eps_gauche(l)/vect_propre_eps_gauche(Nbpt_h);
        UU_d(l)=UU_d(l)*vect_propre_droite(l)/vect_propre_droite(1);

    end % for l

    %val = [ UU_g ; UU_d(2:Nbpt_h)];

elseif num_idee==13

    % declarations
    % ------------
    MM = matM_ref(4*(Nbpt_h-1)+1); % matrice de masse
    KK = sparse(4*(Nbpt_h-1)+1,4*(Nbpt_h-1)+1); % matrice de rigidite
    MM_sigma = sparse(4*(Nbpt_h-1)+1,4*(Nbpt_h-1)+1); % matrice de masse

    for l=2:(4*(Nbpt_h-1)+1-2)   %La première et dernière ligne doivent rester nulle

        KK(l,l)=(1/h)*(A((noeud_i-2)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-2)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-2)*H+(l-1+0.5)*h,epsilon);
        KK(l+1,l)=(-1/h)*A((noeud_i-2)*H+(l-1+0.5)*h,epsilon);

        MM_sigma(l,l)=(h/3)*(Sigma((noeud_i-2)*H+(l-1-0.5)*h,epsilon)+Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon));
        MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon);
        MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon);

    end % for l

    MM_sigma(4*(Nbpt_h-1)+1-1,4*(Nbpt_h-1)+1-1)=(h/3)*(Sigma((noeud_i-2)*H+(4*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+Sigma((noeud_i-2)*H+(4*(Nbpt_h-1)+1-2+0.5)*h,epsilon));
    KK(4*(Nbpt_h-1)+1-1,4*(Nbpt_h-1)+1-1)=(1/h)*(A((noeud_i-2)*H+(4*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+A((noeud_i-2)*H+(4*(Nbpt_h-1)+1-2+0.5)*h,epsilon));

    AA = MM_sigma + KK*(epsilon^2) ;
    AA(1,1)=1; AA(4*(Nbpt_h-1)+1,4*(Nbpt_h-1)+1)=1;

    %[EV,DV] = eigs(AA,MM,5,'smallestabs');
    [EV,DV] = eigs(AA,MM,1,'smallestabs');


    [~,arg_min] = min(diag(DV));
    Vecteur_propre = EV(:,arg_min);
    Vecteur_propre = Vecteur_propre*(sqrt(H)); %On "renormalise" pour que le vecteur soit toujours aux alentours de 1

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre(fix(Nbpt_h/2))<0)
        Vecteur_propre = -Vecteur_propre;
    end

    vect_propre_eps = Vecteur_propre(Nbpt_h:2*Nbpt_h-1);
    % -------------------------------------------------------------------------

    vect_propre_droite = Vecteur_propre(2*Nbpt_h:3*Nbpt_h-1);
    % -------------------------------------------------------------------------


    KK_g = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_g = zeros(Nbpt_h,1);     % vecteur second membre
    KK_d = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_d = zeros(Nbpt_h,1);     % vecteur second membre

    % boucle pour les matrices EF
    % ------------------------
    for l=2:(Nbpt_h-2)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));

        KK_g(l,l)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_g(l,l+1)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK_g(l+1,l)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        KK_d(l,l)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_d(l,l+1)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);
        KK_d(l+1,l)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);

    end % for l
    phi_avant_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1-1));
    phi_apres_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1+1));
    phi_avant_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1-1));
    phi_apres_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1+1));

    KK_g(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
    KK_d(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)

    % Calcul du second membre F
    % -------------------------
    for l=2:(Nbpt_h-1)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));
        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        LL_g(l)=(-1/H)*(-1*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)+phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon));

        LL_d(l)=(1/H)*(-phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)+phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon));

    end % for l

    AA_g = KK_g ;
    AA_g(1,1)=1; AA_g(Nbpt_h,Nbpt_h)=1;

    AA_d = KK_d ;
    AA_d(1,1)=1; AA_d(Nbpt_h,Nbpt_h)=1;
    % inversion
    % ----------

    UU_g = AA_g\LL_g;
    UU_d = AA_d\LL_d;

    % réhaussement
    % ----------
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)+(l-1)*(h/H);
        UU_d(l)=UU_d(l)-(l-1)*(h/H)+1;

    end % for l

    %Multiplication par Psi_eps
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)*vect_propre_eps(l);
        UU_d(l)=UU_d(l)*vect_propre_droite(l);

    end % for l
    %val = [ UU_g ; UU_d(2:Nbpt_h)];
    %     coeff = UU_g(Nbpt_h);
    %     val = val./coeff ;
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif num_idee==14

    % declarations
    % ------------
    MM = matM_ref(2*(Nbpt_h-1)+1); % matrice de masse
    KK = sparse(2*(Nbpt_h-1)+1,2*(Nbpt_h-1)+1); % matrice de rigidite
    MM_sigma = sparse(2*(Nbpt_h-1)+1,2*(Nbpt_h-1)+1); % matrice de masse

    for l=2:(2*(Nbpt_h-1)+1-2)   %La première et dernière ligne doivent rester nulle

        KK(l,l)=(1/h)*(A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK(l+1,l)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        MM_sigma(l,l)=(h/3)*(Sigma((noeud_i-1)*H+(l-1-0.5)*h,epsilon)+Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon));
        MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

    end % for l

    MM_sigma(2*(Nbpt_h-1)+1-1,2*(Nbpt_h-1)+1-1)=(h/3)*(Sigma((noeud_i-1)*H+(2*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+Sigma((noeud_i-1)*H+(2*(Nbpt_h-1)+1-2+0.5)*h,epsilon));
    KK(2*(Nbpt_h-1)+1-1,2*(Nbpt_h-1)+1-1)=(1/h)*(A((noeud_i-1)*H+(2*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+A((noeud_i-1)*H+(2*(Nbpt_h-1)+1-2+0.5)*h,epsilon));

    AA = MM_sigma + KK*(epsilon^2) ;
    AA(1,1)=1; AA(2*(Nbpt_h-1)+1,2*(Nbpt_h-1)+1)=1;

    %[EV,DV] = eigs(AA,MM,5,'smallestabs');
    [EV,DV] = eigs(AA,MM,1,'smallestabs');


    [~,arg_min] = min(diag(DV));
    Vecteur_propre = EV(:,arg_min);
    Vecteur_propre = Vecteur_propre*(sqrt(H)); %On "renormalise" pour que le vecteur soit toujours aux alentours de 1

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre(fix(Nbpt_h))<0)
        Vecteur_propre = -Vecteur_propre;
    end

    vect_propre_eps = Vecteur_propre(1:Nbpt_h);
    % -------------------------------------------------------------------------

    vect_propre_droite = Vecteur_propre(Nbpt_h:2*Nbpt_h-1);
    % -------------------------------------------------------------------------


    KK_g = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_g = zeros(Nbpt_h,1);     % vecteur second membre
    KK_d = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_d = zeros(Nbpt_h,1);     % vecteur second membre

    % boucle pour les matrices EF
    % ------------------------
    for l=2:(Nbpt_h-2)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));

        KK_g(l,l)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_g(l,l+1)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK_g(l+1,l)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        KK_d(l,l)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_d(l,l+1)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);
        KK_d(l+1,l)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);

    end % for l
    phi_avant_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1-1));
    phi_apres_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1+1));
    phi_avant_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1-1));
    phi_apres_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1+1));

    KK_g(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
    KK_d(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)

    % Calcul du second membre F
    % -------------------------
    for l=2:(Nbpt_h-1)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));
        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        LL_g(l)=(-1/H)*(-1*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)+phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon));

        LL_d(l)=(1/H)*(-phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)+phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon));

    end % for l

    AA_g = KK_g ;
    AA_g(1,1)=1; AA_g(Nbpt_h,Nbpt_h)=1;

    AA_d = KK_d ;
    AA_d(1,1)=1; AA_d(Nbpt_h,Nbpt_h)=1;
    % inversion
    % ----------

    UU_g = AA_g\LL_g;
    UU_d = AA_d\LL_d;

    % réhaussement
    % ----------
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)+(l-1)*(h/H);
        UU_d(l)=UU_d(l)-(l-1)*(h/H)+1;

    end % for l

    %Multiplication par Psi_eps
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)*vect_propre_eps(l);
        UU_d(l)=UU_d(l)*vect_propre_droite(l);

    end % for l
    %val = [ UU_g ; UU_d(2:Nbpt_h)];
    %     coeff = UU_g(Nbpt_h);
    %     val = val./coeff ;
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif num_idee==15

    %Calcul de Phi_eps sur [0;3H], puis on sélectionne la partie sur [H,2H]
    % declarations
    % ------------
    MM = matM_ref(3*(Nbpt_h-1)+1); % matrice de masse
    KK = sparse(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1); % matrice de rigidite
    MM_sigma = sparse(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1); % matrice de masse

    for l=2:(3*(Nbpt_h-1)+1-2)   %La première et dernière ligne doivent rester nulle

        KK(l,l)=(1/h)*(A((noeud_i-2)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-2)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-2)*H+(l-1+0.5)*h,epsilon);
        KK(l+1,l)=(-1/h)*A((noeud_i-2)*H+(l-1+0.5)*h,epsilon);

        MM_sigma(l,l)=(h/3)*(Sigma((noeud_i-2)*H+(l-1-0.5)*h,epsilon)+Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon));
        MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon);
        MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon);

    end % for l

    MM_sigma(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=(h/3)*(Sigma((noeud_i-2)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+Sigma((noeud_i-2)*H+(3*(Nbpt_h-1)+1-2+0.5)*h,epsilon));
    KK(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=(1/h)*(A((noeud_i-2)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+A((noeud_i-2)*H+(3*(Nbpt_h-1)+1-2+0.5)*h,epsilon));

    AA = MM_sigma + KK*(epsilon^2) ;
    AA(1,1)=1; AA(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1)=1;

    [EV,DV] = eigs(AA,MM,1,'smallestabs');


    [~,arg_min] = min(diag(DV));
    Vecteur_propre_0_3H = EV(:,arg_min);

    a = norm(Vecteur_propre_0_3H,2)/sqrt(length(Vecteur_propre_0_3H)-1); % Rescale pour une norme L^2 = 1
    Vecteur_propre_0_3H = Vecteur_propre_0_3H./a;

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre_0_3H(fix(Nbpt_h/2))<0)
        Vecteur_propre_0_3H = -Vecteur_propre_0_3H;
    end

    vect_propre_eps_mid = Vecteur_propre_0_3H(Nbpt_h:2*Nbpt_h-1);
    % ------------------------------------------------------------------------------------------------------------------------------------
    %Calcul de U_eps à l'aide du Phi_eps calculé ci-dessus

    KK = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    MM = sparse(Nbpt_h,Nbpt_h); % matrice de masse
    % boucle pour les matrices EF
    % ------------------------
    for l=2:(Nbpt_h-2)   %La première et dernière ligne doivent rester nulle

        phi_avant = 0.5*(vect_propre_eps_mid(l)+vect_propre_eps_mid(l-1));
        phi_apres = 0.5*(vect_propre_eps_mid(l)+vect_propre_eps_mid(l+1));

        KK(l,l)=(1/h)*(phi_avant*phi_avant*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)  +  phi_apres*phi_apres*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*phi_apres*phi_apres*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK(l+1,l)=(-1/h)*phi_apres*phi_apres*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        MM(l,l)=h*(2/3)*(vect_propre_eps_mid(l)*vect_propre_eps_mid(l));
        MM(l,l+1)=h*(1/6)*phi_apres*phi_apres;
        MM(l+1,l)=h*(1/6)*phi_apres*phi_apres;

    end % for l
    phi_avant = 0.5*(vect_propre_eps_mid(Nbpt_h-1)+vect_propre_eps_mid(Nbpt_h-1-1));
    phi_apres = 0.5*(vect_propre_eps_mid(Nbpt_h-1)+vect_propre_eps_mid(Nbpt_h-1+1));

    KK(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant*phi_avant*A((noeud_i-1)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres*phi_apres*A((noeud_i-1)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
    MM(Nbpt_h-1,Nbpt_h-1)=h*(2/3)*(vect_propre_eps_mid(Nbpt_h-1)*vect_propre_eps_mid(Nbpt_h-1));


    AA =  KK*(epsilon^2) ;
    AA(1,1)=1; AA(Nbpt_h,Nbpt_h)=1;


    [EV,DV] = eigs(AA,MM,1,'smallestabs');

    [~,arg_min] = min(diag(DV));
    Vecteur_propre_u_eps = EV(:,arg_min);

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre_u_eps(fix(Nbpt_h/2))<0)
        Vecteur_propre_u_eps = -Vecteur_propre_u_eps;
    end

    a = norm(Vecteur_propre_u_eps,2)/sqrt(length(Vecteur_propre_u_eps)-1); % Rescale pour une norme L^2 = 1
    Vecteur_propre_u_eps = Vecteur_propre_u_eps./a;

    % ------------------------------------------------------------------------------------------------------------------------------------
    %Calcul de Phi_eps sur [0;H]

    [~,Vecteur_propre_0_H] = fonction_propre_eps_Dirichlet(Nbpt_H,Nbpt_h,epsilon,noeud_i);


    Approximation_Psi_eps_g = Vecteur_propre_0_H./Vecteur_propre_u_eps;
    Approximation_Psi_eps_g(1) = Approximation_Psi_eps_g(2);
    Approximation_Psi_eps_g(Nbpt_h) = Approximation_Psi_eps_g(Nbpt_h-1);

    % ----------------------------------------------------------------------
    %Calcul de Phi_eps sur [0;3H], puis on sélectionne la partie sur [H,2H]
    % declarations
    % ------------
    MM = matM_ref(3*(Nbpt_h-1)+1); % matrice de masse
    KK = sparse(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1); % matrice de rigidite
    MM_sigma = sparse(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1); % matrice de masse

    for l=2:(3*(Nbpt_h-1)+1-2)   %La première et dernière ligne doivent rester nulle

        KK(l,l)=(1/h)*(A((noeud_i+1-2)*H+(l-1-0.5)*h,epsilon)+A((noeud_i+1-2)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i+1-2)*H+(l-1+0.5)*h,epsilon);
        KK(l+1,l)=(-1/h)*A((noeud_i+1-2)*H+(l-1+0.5)*h,epsilon);

        MM_sigma(l,l)=(h/3)*(Sigma((noeud_i+1-2)*H+(l-1-0.5)*h,epsilon)+Sigma((noeud_i+1-2)*H+(l-1+0.5)*h,epsilon));
        MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i+1-2)*H+(l-1+0.5)*h,epsilon);
        MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i+1-2)*H+(l-1+0.5)*h,epsilon);

    end % for l

    MM_sigma(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=(h/3)*(Sigma((noeud_i+1-2)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+Sigma((noeud_i+1-2)*H+(3*(Nbpt_h-1)+1-2+0.5)*h,epsilon));
    KK(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=(1/h)*(A((noeud_i+1-2)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+A((noeud_i+1-2)*H+(3*(Nbpt_h-1)+1-2+0.5)*h,epsilon));

    AA = MM_sigma + KK*(epsilon^2) ;
    AA(1,1)=1; AA(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1)=1;

    [EV,DV] = eigs(AA,MM,1,'smallestabs');


    [~,arg_min] = min(diag(DV));
    Vecteur_propre_0_3H = EV(:,arg_min);

    a = norm(Vecteur_propre_0_3H,2)/sqrt(length(Vecteur_propre_0_3H)-1); % Rescale pour une norme L^2 = 1
    Vecteur_propre_0_3H = Vecteur_propre_0_3H./a;

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre_0_3H(fix(Nbpt_h/2))<0)
        Vecteur_propre_0_3H = -Vecteur_propre_0_3H;
    end

    vect_propre_eps_mid = Vecteur_propre_0_3H(Nbpt_h:2*Nbpt_h-1);
    % ------------------------------------------------------------------------------------------------------------------------------------
    %Calcul de U_eps à l'aide du Phi_eps calculé ci-dessus

    KK = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    MM = sparse(Nbpt_h,Nbpt_h); % matrice de masse
    % boucle pour les matrices EF
    % ------------------------
    for l=2:(Nbpt_h-2)   %La première et dernière ligne doivent rester nulle

        phi_avant = 0.5*(vect_propre_eps_mid(l)+vect_propre_eps_mid(l-1));
        phi_apres = 0.5*(vect_propre_eps_mid(l)+vect_propre_eps_mid(l+1));

        KK(l,l)=(1/h)*(phi_avant*phi_avant*A((noeud_i+1-1)*H+(l-1-0.5)*h,epsilon)  +  phi_apres*phi_apres*A((noeud_i+1-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*phi_apres*phi_apres*A((noeud_i+1-1)*H+(l-1+0.5)*h,epsilon);
        KK(l+1,l)=(-1/h)*phi_apres*phi_apres*A((noeud_i+1-1)*H+(l-1+0.5)*h,epsilon);

        MM(l,l)=h*(2/3)*(vect_propre_eps_mid(l)*vect_propre_eps_mid(l));
        MM(l,l+1)=h*(1/6)*phi_apres*phi_apres;
        MM(l+1,l)=h*(1/6)*phi_apres*phi_apres;

    end % for l
    phi_avant = 0.5*(vect_propre_eps_mid(Nbpt_h-1)+vect_propre_eps_mid(Nbpt_h-1-1));
    phi_apres = 0.5*(vect_propre_eps_mid(Nbpt_h-1)+vect_propre_eps_mid(Nbpt_h-1+1));

    KK(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant*phi_avant*A((noeud_i+1-1)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres*phi_apres*A((noeud_i+1-1)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
    MM(Nbpt_h-1,Nbpt_h-1)=h*(2/3)*(vect_propre_eps_mid(Nbpt_h-1)*vect_propre_eps_mid(Nbpt_h-1));


    AA =  KK*(epsilon^2) ;
    AA(1,1)=1; AA(Nbpt_h,Nbpt_h)=1;


    [EV,DV] = eigs(AA,MM,1,'smallestabs');

    [~,arg_min] = min(diag(DV));
    Vecteur_propre_u_eps = EV(:,arg_min);

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre_u_eps(fix(Nbpt_h/2))<0)
        Vecteur_propre_u_eps = -Vecteur_propre_u_eps;
    end

    a = norm(Vecteur_propre_u_eps,2)/sqrt(length(Vecteur_propre_u_eps)-1); % Rescale pour une norme L^2 = 1
    Vecteur_propre_u_eps = Vecteur_propre_u_eps./a;

    % ------------------------------------------------------------------------------------------------------------------------------------
    %Calcul de Phi_eps sur [0;H]

    [~,Vecteur_propre_0_H] = fonction_propre_eps_Dirichlet(Nbpt_H,Nbpt_h,epsilon,noeud_i+1);


    Approximation_Psi_eps_d = Vecteur_propre_0_H./Vecteur_propre_u_eps;
    Approximation_Psi_eps_d(1) = Approximation_Psi_eps_d(2);
    Approximation_Psi_eps_d(Nbpt_h) = Approximation_Psi_eps_d(Nbpt_h-1);



    vect_propre_eps = Approximation_Psi_eps_g;
    % -------------------------------------------------------------------------

    vect_propre_droite = Approximation_Psi_eps_d;
    % -------------------------------------------------------------------------




    KK_g = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_g = zeros(Nbpt_h,1);     % vecteur second membre
    KK_d = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_d = zeros(Nbpt_h,1);     % vecteur second membre

    % boucle pour les matrices EF
    % ------------------------
    for l=2:(Nbpt_h-2)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));

        KK_g(l,l)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_g(l,l+1)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK_g(l+1,l)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        KK_d(l,l)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_d(l,l+1)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);
        KK_d(l+1,l)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);

    end % for l
    phi_avant_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1-1));
    phi_apres_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1+1));
    phi_avant_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1-1));
    phi_apres_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1+1));

    KK_g(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
    KK_d(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)

    % Calcul du second membre F
    % -------------------------
    for l=2:(Nbpt_h-1)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));
        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        LL_g(l)=(-1/H)*(-1*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)+phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon));

        LL_d(l)=(1/H)*(-phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)+phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon));

    end % for l

    AA_g = KK_g ;
    AA_g(1,1)=1; AA_g(Nbpt_h,Nbpt_h)=1;

    AA_d = KK_d ;
    AA_d(1,1)=1; AA_d(Nbpt_h,Nbpt_h)=1;
    % inversion
    % ----------

    UU_g = AA_g\LL_g;
    UU_d = AA_d\LL_d;

    % réhaussement
    % ----------
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)+(l-1)*(h/H);
        UU_d(l)=UU_d(l)-(l-1)*(h/H)+1;

    end % for l

    %Multiplication par Psi_eps
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)*vect_propre_eps(l);
        UU_d(l)=UU_d(l)*vect_propre_droite(l);

    end % for l
    %val = [ UU_g ; UU_d(2:Nbpt_h)];
    %     coeff = UU_g(Nbpt_h);
    %     val = val./coeff ;
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

elseif num_idee==16

    % declarations
    % ------------
    MM = matM_ref(9*(Nbpt_h-1)+1); % matrice de masse
    KK = sparse(9*(Nbpt_h-1)+1,9*(Nbpt_h-1)+1); % matrice de rigidite
    MM_sigma = sparse(9*(Nbpt_h-1)+1,9*(Nbpt_h-1)+1); % matrice de masse

    for l=2:(9*(Nbpt_h-1)+1-2)   %La première et dernière ligne doivent rester nulle

        KK(l,l)=(1/h)*(A((noeud_i-5)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-5)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-5)*H+(l-1+0.5)*h,epsilon);
        KK(l+1,l)=(-1/h)*A((noeud_i-5)*H+(l-1+0.5)*h,epsilon);

        MM_sigma(l,l)=(h/3)*(Sigma((noeud_i-5)*H+(l-1-0.5)*h,epsilon)+Sigma((noeud_i-5)*H+(l-1+0.5)*h,epsilon));
        MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i-5)*H+(l-1+0.5)*h,epsilon);
        MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i-5)*H+(l-1+0.5)*h,epsilon);

    end % for l

    MM_sigma(9*(Nbpt_h-1)+1-1,9*(Nbpt_h-1)+1-1)=(h/3)*(Sigma((noeud_i-5)*H+(9*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+Sigma((noeud_i-5)*H+(9*(Nbpt_h-1)+1-2+0.5)*h,epsilon));
    KK(9*(Nbpt_h-1)+1-1,9*(Nbpt_h-1)+1-1)=(1/h)*(A((noeud_i-5)*H+(9*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+A((noeud_i-5)*H+(9*(Nbpt_h-1)+1-2+0.5)*h,epsilon));

    AA = MM_sigma + KK*(epsilon^2) ;
    AA(1,1)=1; AA(9*(Nbpt_h-1)+1,9*(Nbpt_h-1)+1)=1;

    %[EV,DV] = eigs(AA,MM,5,'smallestabs');
    [EV,DV] = eigs(AA,MM,1,'smallestabs');


    [~,arg_min] = min(diag(DV));
    Vecteur_propre = EV(:,arg_min);
    Vecteur_propre = Vecteur_propre*(sqrt(H)); %On "renormalise" pour que le vecteur soit toujours aux alentours de 1

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre(fix(9*Nbpt_h/2))<0)
        Vecteur_propre = -Vecteur_propre;
    end

    vect_propre_eps = Vecteur_propre(4*(Nbpt_h-1)+1:5*(Nbpt_h-1)+1);
    % -------------------------------------------------------------------------

    % declarations
    % ------------
    MM = matM_ref(9*(Nbpt_h-1)+1); % matrice de masse
    KK = sparse(9*(Nbpt_h-1)+1,9*(Nbpt_h-1)+1); % matrice de rigidite
    MM_sigma = sparse(9*(Nbpt_h-1)+1,9*(Nbpt_h-1)+1); % matrice de masse

    for l=2:(9*(Nbpt_h-1)+1-2)   %La première et dernière ligne doivent rester nulle

        KK(l,l)=(1/h)*(A((noeud_i-4)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-4)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-4)*H+(l-1+0.5)*h,epsilon);
        KK(l+1,l)=(-1/h)*A((noeud_i-4)*H+(l-1+0.5)*h,epsilon);

        MM_sigma(l,l)=(h/3)*(Sigma((noeud_i-4)*H+(l-1-0.5)*h,epsilon)+Sigma((noeud_i-4)*H+(l-1+0.5)*h,epsilon));
        MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i-4)*H+(l-1+0.5)*h,epsilon);
        MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i-4)*H+(l-1+0.5)*h,epsilon);

    end % for l

    MM_sigma(9*(Nbpt_h-1)+1-1,9*(Nbpt_h-1)+1-1)=(h/3)*(Sigma((noeud_i-4)*H+(9*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+Sigma((noeud_i-4)*H+(9*(Nbpt_h-1)+1-2+0.5)*h,epsilon));
    KK(9*(Nbpt_h-1)+1-1,9*(Nbpt_h-1)+1-1)=(1/h)*(A((noeud_i-4)*H+(9*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+A((noeud_i-4)*H+(9*(Nbpt_h-1)+1-2+0.5)*h,epsilon));

    AA = MM_sigma + KK*(epsilon^2) ;
    AA(1,1)=1; AA(9*(Nbpt_h-1)+1,9*(Nbpt_h-1)+1)=1;

    %[EV,DV] = eigs(AA,MM,5,'smallestabs');
    [EV,DV] = eigs(AA,MM,1,'smallestabs');


    [~,arg_min] = min(diag(DV));
    Vecteur_propre = EV(:,arg_min);
    Vecteur_propre = Vecteur_propre*(sqrt(H)); %On "renormalise" pour que le vecteur soit toujours aux alentours de 1

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre(fix(9*Nbpt_h/2))<0)
        Vecteur_propre = -Vecteur_propre;
    end

    vect_propre_droite = Vecteur_propre(4*(Nbpt_h-1)+1:5*(Nbpt_h-1)+1);
    % -------------------------------------------------------------------------
    KK_g = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_g = zeros(Nbpt_h,1);     % vecteur second membre
    KK_d = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_d = zeros(Nbpt_h,1);     % vecteur second membre

    % boucle pour les matrices EF
    % ------------------------
    for l=2:(Nbpt_h-2)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));

        KK_g(l,l)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_g(l,l+1)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK_g(l+1,l)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        KK_d(l,l)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_d(l,l+1)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);
        KK_d(l+1,l)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);

    end % for l
    phi_avant_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1-1));
    phi_apres_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1+1));
    phi_avant_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1-1));
    phi_apres_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1+1));

    KK_g(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
    KK_d(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)

    % Calcul du second membre F
    % -------------------------
    for l=2:(Nbpt_h-1)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));
        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        LL_g(l)=(-1/H)*(-1*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)+phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon));

        LL_d(l)=(1/H)*(-phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)+phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon));

    end % for l

    AA_g = KK_g ;
    AA_g(1,1)=1; AA_g(Nbpt_h,Nbpt_h)=1;

    AA_d = KK_d ;
    AA_d(1,1)=1; AA_d(Nbpt_h,Nbpt_h)=1;
    % inversion
    % ----------

    UU_g = AA_g\LL_g;
    UU_d = AA_d\LL_d;

    % réhaussement
    % ----------
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)+(l-1)*(h/H);
        UU_d(l)=UU_d(l)-(l-1)*(h/H)+1;

    end % for l

    %Multiplication par Psi_eps
    for l=1:Nbpt_h

%         UU_g(l)=UU_g(l)*vect_propre_eps(l)/vect_propre_eps(Nbpt_h);
%         UU_d(l)=UU_d(l)*vect_propre_droite(l)/vect_propre_droite(1);
                UU_g(l)=UU_g(l)*vect_propre_eps(l);
                UU_d(l)=UU_d(l)*vect_propre_droite(l);

    end % for l
    %val = [ UU_g ; UU_d(2:Nbpt_h)];
    %     coeff = UU_g(Nbpt_h);
    %     val = val./coeff ;
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

elseif num_idee==17

    % declarations
    % ------------
    MM = matM_ref(3*(Nbpt_h-1)+1); % matrice de masse
    KK = sparse(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1); % matrice de rigidite
    MM_sigma = sparse(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1); % matrice de masse

    for l=2:(3*(Nbpt_h-1)+1-2)   %La première et dernière ligne doivent rester nulle

        KK(l,l)=(1/h)*(A((noeud_i-2)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-2)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-2)*H+(l-1+0.5)*h,epsilon);
        KK(l+1,l)=(-1/h)*A((noeud_i-2)*H+(l-1+0.5)*h,epsilon);

        MM_sigma(l,l)=(h/3)*(Sigma((noeud_i-2)*H+(l-1-0.5)*h,epsilon)+Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon));
        MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon);
        MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon);

    end % for l

    MM_sigma(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=(h/3)*(Sigma((noeud_i-2)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+Sigma((noeud_i-2)*H+(3*(Nbpt_h-1)+1-2+0.5)*h,epsilon));
    KK(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=(1/h)*(A((noeud_i-2)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+A((noeud_i-2)*H+(3*(Nbpt_h-1)+1-2+0.5)*h,epsilon));

    AA = MM_sigma + KK*(epsilon^2) ;
    AA(1,1)=1; AA(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1)=1;

    %[EV,DV] = eigs(AA,MM,5,'smallestabs');
    [EV,DV] = eigs(AA,MM,1,'smallestabs');


    [~,arg_min] = min(diag(DV));
    Vecteur_propre = EV(:,arg_min);
    Vecteur_propre = Vecteur_propre*(sqrt(H)); %On "renormalise" pour que le vecteur soit toujours aux alentours de 1

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre(fix(Nbpt_h/2))<0)
        Vecteur_propre = -Vecteur_propre;
    end

    vect_propre_gauche = zeros(Nbpt_h,1);

    %Calcul de la moyenne locale, afin d'enlever la sinus macroscopique
    for l=Nbpt_h:(2*Nbpt_h-1)

        yl_min = l - (fix(epsilon*((Nbpt_ref-1)))+1);
        yl_max = l + (fix(epsilon*((Nbpt_ref-1))));

        Moyenne_local = 0;
        for j=yl_min:yl_max
            %Moyenne_local = Moyenne_local + Vecteur_propre(j);
            Moyenne_local = Moyenne_local + Vecteur_propre(j)*Vecteur_propre(j);
        end
        Moyenne_local = Moyenne_local/(yl_max-yl_min+1);
        vect_propre_gauche(l-Nbpt_h+1) = Vecteur_propre(l)/sqrt(Moyenne_local);

    end


    %vect_propre_gauche = Vecteur_propre(Nbpt_h:2*Nbpt_h-1);
    % -------------------------------------------------------------------------




    % declarations
    % ------------
    MM = matM_ref(3*(Nbpt_h-1)+1); % matrice de masse
    KK = sparse(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1); % matrice de rigidite
    MM_sigma = sparse(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1); % matrice de masse

    for l=2:(3*(Nbpt_h-1)+1-2)   %La première et dernière ligne doivent rester nulle

        KK(l,l)=(1/h)*(A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK(l+1,l)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        MM_sigma(l,l)=(h/3)*(Sigma((noeud_i-1)*H+(l-1-0.5)*h,epsilon)+Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon));
        MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

    end % for l

    MM_sigma(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=(h/3)*(Sigma((noeud_i-1)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+Sigma((noeud_i-1)*H+(3*(Nbpt_h-1)+1-2+0.5)*h,epsilon));
    KK(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=(1/h)*(A((noeud_i-1)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+A((noeud_i-1)*H+(3*(Nbpt_h-1)+1-2+0.5)*h,epsilon));

    AA = MM_sigma + KK*(epsilon^2) ;
    AA(1,1)=1; AA(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1)=1;

    %[EV,DV] = eigs(AA,MM,5,'smallestabs');
    [EV,DV] = eigs(AA,MM,1,'smallestabs');


    [~,arg_min] = min(diag(DV));
    Vecteur_propre = EV(:,arg_min);
    Vecteur_propre = Vecteur_propre*(sqrt(H)); %On "renormalise" pour que le vecteur soit toujours aux alentours de 1

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre(fix(Nbpt_h/2))<0)
        Vecteur_propre = -Vecteur_propre;
    end

    vect_propre_droite = zeros(Nbpt_h,1);

    %Calcul de la moyenne locale, afin d'enlever la sinus macroscopique
    for l=Nbpt_h:(2*Nbpt_h-1)

        yl_min = l - (fix(epsilon*((Nbpt_ref-1)))+1);
        yl_max = l + (fix(epsilon*((Nbpt_ref-1))));

        Moyenne_local = 0;
        for j=yl_min:yl_max
            %Moyenne_local = Moyenne_local + Vecteur_propre(j);
            Moyenne_local = Moyenne_local + Vecteur_propre(j)*Vecteur_propre(j);
        end
        Moyenne_local = Moyenne_local/(yl_max-yl_min+1);
        vect_propre_droite(l-Nbpt_h+1) = Vecteur_propre(l)/sqrt(Moyenne_local);

    end

    %vect_propre_droite = Vecteur_propre(Nbpt_h:2*Nbpt_h-1);
    % -------------------------------------------------------------------------

    KK_g = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_g = zeros(Nbpt_h,1);     % vecteur second membre
    KK_d = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_d = zeros(Nbpt_h,1);     % vecteur second membre

    % boucle pour les matrices EF
    % ------------------------
    for l=2:(Nbpt_h-2)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_gauche(l)+vect_propre_gauche(l-1));
        phi_apres_g = 0.5*(vect_propre_gauche(l)+vect_propre_gauche(l+1));

        KK_g(l,l)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_g(l,l+1)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK_g(l+1,l)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        KK_d(l,l)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_d(l,l+1)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);
        KK_d(l+1,l)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);

    end % for l
    phi_avant_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1-1));
    phi_apres_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1+1));
    phi_avant_g = 0.5*(vect_propre_gauche(Nbpt_h-1)+vect_propre_gauche(Nbpt_h-1-1));
    phi_apres_g = 0.5*(vect_propre_gauche(Nbpt_h-1)+vect_propre_gauche(Nbpt_h-1+1));

    KK_g(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
    KK_d(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)

    % Calcul du second membre F
    % -------------------------
    for l=2:(Nbpt_h-1)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_gauche(l)+vect_propre_gauche(l-1));
        phi_apres_g = 0.5*(vect_propre_gauche(l)+vect_propre_gauche(l+1));
        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        LL_g(l)=(-1/H)*(-1*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)+phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon));

        LL_d(l)=(1/H)*(-phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)+phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon));

    end % for l

    AA_g = KK_g ;
    AA_g(1,1)=1; AA_g(Nbpt_h,Nbpt_h)=1;

    AA_d = KK_d ;
    AA_d(1,1)=1; AA_d(Nbpt_h,Nbpt_h)=1;
    % inversion
    % ----------

    UU_g = AA_g\LL_g;
    UU_d = AA_d\LL_d;

    % réhaussement
    % ----------
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)+(l-1)*(h/H);
        UU_d(l)=UU_d(l)-(l-1)*(h/H)+1;

    end % for l

    %Multiplication par Psi_eps
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)*vect_propre_gauche(l)/vect_propre_gauche(Nbpt_h);
        UU_d(l)=UU_d(l)*vect_propre_droite(l)/vect_propre_droite(1);
%                 UU_g(l)=UU_g(l)*vect_propre_gauche(l);
%                 UU_d(l)=UU_d(l)*vect_propre_droite(l);

    end % for l
    %val = [ UU_g , UU_d];
    
    %     coeff = UU_g(Nbpt_h);
    %     val = val./coeff ;
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

elseif num_idee==18

    % declarations
    % ------------
    MM = matM_ref(3*(Nbpt_h-1)+1); % matrice de masse
    KK = sparse(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1); % matrice de rigidite
    MM_sigma = sparse(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1); % matrice de masse

    for l=2:(3*(Nbpt_h-1)+1-2)   %La première et dernière ligne doivent rester nulle

        KK(l,l)=(1/h)*(A((noeud_i-2)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-2)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-2)*H+(l-1+0.5)*h,epsilon);
        KK(l+1,l)=(-1/h)*A((noeud_i-2)*H+(l-1+0.5)*h,epsilon);

        MM_sigma(l,l)=(h/3)*(Sigma((noeud_i-2)*H+(l-1-0.5)*h,epsilon)+Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon));
        MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon);
        MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon);

    end % for l

    MM_sigma(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=(h/3)*(Sigma((noeud_i-2)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+Sigma((noeud_i-2)*H+(3*(Nbpt_h-1)+1-2+0.5)*h,epsilon));
    KK(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=(1/h)*(A((noeud_i-2)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+A((noeud_i-2)*H+(3*(Nbpt_h-1)+1-2+0.5)*h,epsilon));

    AA = MM_sigma + KK*(epsilon^2) ;
    AA(1,1)=1; AA(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1)=1;

    %[EV,DV] = eigs(AA,MM,5,'smallestabs');
    [EV,DV] = eigs(AA,MM,1,'smallestabs');


    [~,arg_min] = min(diag(DV));
    Vecteur_propre = EV(:,arg_min);
    Vecteur_propre = Vecteur_propre*(sqrt(H)); %On "renormalise" pour que le vecteur soit toujours aux alentours de 1

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre(fix(Nbpt_h/2))<0)
        Vecteur_propre = -Vecteur_propre;
    end

    a = sin(pi*(0:(3*(Nbpt_h-1)))/(3*(Nbpt_h-1)));

    a = a* norm(Vecteur_propre,2)/norm(a,2);
    Vect_propre_sinus_triche = Vecteur_propre./a';

    vect_propre_gauche = Vect_propre_sinus_triche(Nbpt_h:2*Nbpt_h-1);


    % declarations
    % ------------
    MM = matM_ref(3*(Nbpt_h-1)+1); % matrice de masse
    KK = sparse(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1); % matrice de rigidite
    MM_sigma = sparse(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1); % matrice de masse

    for l=2:(3*(Nbpt_h-1)+1-2)   %La première et dernière ligne doivent rester nulle

        KK(l,l)=(1/h)*(A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK(l+1,l)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        MM_sigma(l,l)=(h/3)*(Sigma((noeud_i-1)*H+(l-1-0.5)*h,epsilon)+Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon));
        MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

    end % for l

    MM_sigma(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=(h/3)*(Sigma((noeud_i-1)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+Sigma((noeud_i-1)*H+(3*(Nbpt_h-1)+1-2+0.5)*h,epsilon));
    KK(3*(Nbpt_h-1)+1-1,3*(Nbpt_h-1)+1-1)=(1/h)*(A((noeud_i-1)*H+(3*(Nbpt_h-1)+1-2-0.5)*h,epsilon)+A((noeud_i-1)*H+(3*(Nbpt_h-1)+1-2+0.5)*h,epsilon));

    AA = MM_sigma + KK*(epsilon^2) ;
    AA(1,1)=1; AA(3*(Nbpt_h-1)+1,3*(Nbpt_h-1)+1)=1;

    %[EV,DV] = eigs(AA,MM,5,'smallestabs');
    [EV,DV] = eigs(AA,MM,1,'smallestabs');


    [~,arg_min] = min(diag(DV));
    Vecteur_propre = EV(:,arg_min);
    Vecteur_propre = Vecteur_propre*(sqrt(H)); %On "renormalise" pour que le vecteur soit toujours aux alentours de 1

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre(fix(Nbpt_h/2))<0)
        Vecteur_propre = -Vecteur_propre;
    end

    a = sin(pi*(0:(3*(Nbpt_h-1)))/(3*(Nbpt_h-1)));

    a = a* norm(Vecteur_propre,2)/norm(a,2);
    Vect_propre_sinus_triche = Vecteur_propre./a';

    vect_propre_droite = Vect_propre_sinus_triche(Nbpt_h:2*Nbpt_h-1);

    KK_g = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_g = zeros(Nbpt_h,1);     % vecteur second membre
    KK_d = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_d = zeros(Nbpt_h,1);     % vecteur second membre

    % boucle pour les matrices EF
    % ------------------------
    for l=2:(Nbpt_h-2)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_gauche(l)+vect_propre_gauche(l-1));
        phi_apres_g = 0.5*(vect_propre_gauche(l)+vect_propre_gauche(l+1));

        KK_g(l,l)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_g(l,l+1)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK_g(l+1,l)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        KK_d(l,l)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_d(l,l+1)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);
        KK_d(l+1,l)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);

    end % for l
    phi_avant_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1-1));
    phi_apres_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1+1));
    phi_avant_g = 0.5*(vect_propre_gauche(Nbpt_h-1)+vect_propre_gauche(Nbpt_h-1-1));
    phi_apres_g = 0.5*(vect_propre_gauche(Nbpt_h-1)+vect_propre_gauche(Nbpt_h-1+1));

    KK_g(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
    KK_d(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)

    % Calcul du second membre F
    % -------------------------
    for l=2:(Nbpt_h-1)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_gauche(l)+vect_propre_gauche(l-1));
        phi_apres_g = 0.5*(vect_propre_gauche(l)+vect_propre_gauche(l+1));
        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        LL_g(l)=(-1/H)*(-1*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)+phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon));

        LL_d(l)=(1/H)*(-phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)+phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon));

    end % for l

    AA_g = KK_g ;
    AA_g(1,1)=1; AA_g(Nbpt_h,Nbpt_h)=1;

    AA_d = KK_d ;
    AA_d(1,1)=1; AA_d(Nbpt_h,Nbpt_h)=1;
    % inversion
    % ----------

    UU_g = AA_g\LL_g;
    UU_d = AA_d\LL_d;

    % réhaussement
    % ----------
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)+(l-1)*(h/H);
        UU_d(l)=UU_d(l)-(l-1)*(h/H)+1;

    end % for l

    %Multiplication par Psi_eps
    for l=1:Nbpt_h

        % UU_g(l)=UU_g(l)*vect_propre_gauche(l)/vect_propre_gauche(Nbpt_h);
        % UU_d(l)=UU_d(l)*vect_propre_droite(l)/vect_propre_droite(1);
        UU_g(l)=UU_g(l)*vect_propre_gauche(l);
        UU_d(l)=UU_d(l)*vect_propre_droite(l);

    end % for l
    %val = [ UU_g ; UU_d(2:Nbpt_h)];
    %     coeff = UU_g(Nbpt_h);
    %     val = val./coeff ;
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif num_idee==19

    MM_ref = matM_ref(Nbpt_h);
    Nbpt_per = 3*(Nbpt_h-1)+1-1;
    %Matrice EF
    MM = matM_per(Nbpt_per);
    MM_sigma = sparse(Nbpt_per,Nbpt_per); % matrice de masse avec sigma
    KK = sparse(Nbpt_per,Nbpt_per); % matrice de rigidite

    for l=2:(Nbpt_per-1)   %La première et dernière ligne doivent rester nulle

        KK(l,l)=(1/h)*(A((noeud_i-2)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-2)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-2)*H+(l-1+0.5)*h,epsilon); %#ok<*SPRIX>
        KK(l+1,l)=(-1/h)*A((noeud_i-2)*H+(l-1+0.5)*h,epsilon);

        MM_sigma(l,l)=(h/3)*(Sigma((noeud_i-2)*H+(l-1-0.5)*h,epsilon)+Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon));
        MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon);
        MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon);

    end % for l

    KK(1,1)=(1/h)*(A((noeud_i-2)*H+0.5*h,epsilon)+A((noeud_i+1)*H-0.5*h,epsilon));
    KK(1,2)=(-1/h)*A((noeud_i-2)*H+0.5*h,epsilon);
    KK(2,1)=(-1/h)*A((noeud_i-2)*H+0.5*h,epsilon);
    KK(1,Nbpt_per)=(-1/h)*A((noeud_i+1)*H-0.5*h,epsilon);
    KK(Nbpt_per,1)=(-1/h)*A((noeud_i+1)*H-0.5*h,epsilon);
    KK(Nbpt_per,Nbpt_per)=(1/h)*(A((noeud_i-2)*H+(Nbpt_per-1-0.5)*h,epsilon)+A((noeud_i-2)*H+(Nbpt_per-1+0.5)*h,epsilon));

    MM_sigma(1,1)=(h/3)*(Sigma((noeud_i-2)*H+0.5*h,epsilon)+Sigma((noeud_i+1)*H-0.5*h,epsilon));
    MM_sigma(1,2)=h*(1/6)*Sigma((noeud_i-2)*H+0.5*h,epsilon);
    MM_sigma(2,1)=h*(1/6)*Sigma((noeud_i-2)*H+0.5*h,epsilon);
    MM_sigma(1,Nbpt_per)=h*(1/6)*Sigma((noeud_i+1)*H-0.5*h,epsilon);
    MM_sigma(Nbpt_per,1)=h*(1/6)*Sigma((noeud_i+1)*H-0.5*h,epsilon);
    MM_sigma(Nbpt_per,Nbpt_per)=(h/3)*(Sigma((noeud_i-2)*H+(Nbpt_per-1-0.5)*h,epsilon)+Sigma((noeud_i-2)*H+(Nbpt_per-1+0.5)*h,epsilon));


    AA = MM_sigma + KK*(epsilon^2) ;

    [EV,DV] = eigs(AA,MM,1,'smallestabs');

    [~,arg_min] = min(diag(DV));
    Vecteur_propre_H_periodic = EV(:,arg_min);

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre_H_periodic(fix(Nbpt_h/2))<0)
        Vecteur_propre_H_periodic = -Vecteur_propre_H_periodic;
    end

    Vecteur_propre_H_periodic = [Vecteur_propre_H_periodic;Vecteur_propre_H_periodic(1)];

    Vecteur_propre_H_periodic_mid = Vecteur_propre_H_periodic(Nbpt_h:2*Nbpt_h-1);

    % % Rescale pour une norme L^2 = 1
    Norme_L2 = sqrt(Vecteur_propre_H_periodic_mid'*MM_ref*Vecteur_propre_H_periodic_mid);
    vect_propre_eps = Vecteur_propre_H_periodic_mid/Norme_L2;
    % -------------------------------------------------------------------------

    MM_sigma = sparse(Nbpt_per,Nbpt_per); % matrice de masse avec sigma
    KK = sparse(Nbpt_per,Nbpt_per); % matrice de rigidite

    for l=2:(Nbpt_per-1)   %La première et dernière ligne doivent rester nulle

        KK(l,l)=(1/h)*(A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon); %#ok<*SPRIX>
        KK(l+1,l)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        MM_sigma(l,l)=(h/3)*(Sigma((noeud_i-1)*H+(l-1-0.5)*h,epsilon)+Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon));
        MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

    end % for l

    KK(1,1)=(1/h)*(A((noeud_i-1)*H+0.5*h,epsilon)+A((noeud_i+2)*H-0.5*h,epsilon));
    KK(1,2)=(-1/h)*A((noeud_i-1)*H+0.5*h,epsilon);
    KK(2,1)=(-1/h)*A((noeud_i-1)*H+0.5*h,epsilon);
    KK(1,Nbpt_per)=(-1/h)*A((noeud_i+2)*H-0.5*h,epsilon);
    KK(Nbpt_per,1)=(-1/h)*A((noeud_i+2)*H-0.5*h,epsilon);
    KK(Nbpt_per,Nbpt_per)=(1/h)*(A((noeud_i-1)*H+(Nbpt_per-1-0.5)*h,epsilon)+A((noeud_i-1)*H+(Nbpt_per-1+0.5)*h,epsilon));

    MM_sigma(1,1)=(h/3)*(Sigma((noeud_i-1)*H+0.5*h,epsilon)+Sigma((noeud_i+2)*H-0.5*h,epsilon));
    MM_sigma(1,2)=h*(1/6)*Sigma((noeud_i-1)*H+0.5*h,epsilon);
    MM_sigma(2,1)=h*(1/6)*Sigma((noeud_i-1)*H+0.5*h,epsilon);
    MM_sigma(1,Nbpt_per)=h*(1/6)*Sigma((noeud_i+2)*H-0.5*h,epsilon);
    MM_sigma(Nbpt_per,1)=h*(1/6)*Sigma((noeud_i+2)*H-0.5*h,epsilon);
    MM_sigma(Nbpt_per,Nbpt_per)=(h/3)*(Sigma((noeud_i-1)*H+(Nbpt_per-1-0.5)*h,epsilon)+Sigma((noeud_i-1)*H+(Nbpt_per-1+0.5)*h,epsilon));


    AA = MM_sigma + KK*(epsilon^2) ;

    [EV,DV] = eigs(AA,MM,1,'smallestabs');

    [~,arg_min] = min(diag(DV));
    Vecteur_propre_H_periodic = EV(:,arg_min);

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre_H_periodic(fix(Nbpt_h/2))<0)
        Vecteur_propre_H_periodic = -Vecteur_propre_H_periodic;
    end

    Vecteur_propre_H_periodic = [Vecteur_propre_H_periodic;Vecteur_propre_H_periodic(1)];

    Vecteur_propre_H_periodic_mid = Vecteur_propre_H_periodic(Nbpt_h:2*Nbpt_h-1);

    % % Rescale pour une norme L^2 = 1
    Norme_L2 = sqrt(Vecteur_propre_H_periodic_mid'*MM_ref*Vecteur_propre_H_periodic_mid);
    vect_propre_droite = Vecteur_propre_H_periodic_mid/Norme_L2;
    % -------------------------------------------------------------------------

    % -------------------------------------------------------------------------

    KK_g = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_g = zeros(Nbpt_h,1);     % vecteur second membre
    KK_d = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_d = zeros(Nbpt_h,1);     % vecteur second membre

    % boucle pour les matrices EF
    % ------------------------
    for l=2:(Nbpt_h-2)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));

        KK_g(l,l)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_g(l,l+1)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK_g(l+1,l)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        KK_d(l,l)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_d(l,l+1)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);
        KK_d(l+1,l)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);

    end % for l
    phi_avant_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1-1));
    phi_apres_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1+1));
    phi_avant_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1-1));
    phi_apres_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1+1));

    KK_g(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
    KK_d(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)

    % Calcul du second membre F
    % -------------------------
    for l=2:(Nbpt_h-1)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));
        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        LL_g(l)=(-1/H)*(-1*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)+phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon));

        LL_d(l)=(1/H)*(-phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)+phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon));

    end % for l

    AA_g = KK_g ;
    AA_g(1,1)=1; AA_g(Nbpt_h,Nbpt_h)=1;

    AA_d = KK_d ;
    AA_d(1,1)=1; AA_d(Nbpt_h,Nbpt_h)=1;
    % inversion
    % ----------

    UU_g = AA_g\LL_g;
    UU_d = AA_d\LL_d;

    % réhaussement
    % ----------
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)+(l-1)*(h/H);
        UU_d(l)=UU_d(l)-(l-1)*(h/H)+1;

    end % for l

    %Multiplication par Psi_eps
    for l=1:Nbpt_h

%         UU_g(l)=UU_g(l)*vect_propre_eps(l)/vect_propre_eps(Nbpt_h);
%         UU_d(l)=UU_d(l)*vect_propre_droite(l)/vect_propre_droite(1);
                UU_g(l)=UU_g(l)*vect_propre_eps(l);
                UU_d(l)=UU_d(l)*vect_propre_droite(l);

    end % for l
    %val = [ UU_g ; UU_d(2:Nbpt_h)];
    %     coeff = UU_g(Nbpt_h);
    %     val = val./coeff ;
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif num_idee==20

    [~,Psi] = Solution_pb_spectral(Nbpt_h);
    vect_propre_eps = zeros(Nbpt_h,1);
    vect_propre_droite = zeros(Nbpt_h,1);
    for i=1:Nbpt_h
        xi_g = (noeud_i-1)*H+(i-1)*h ;
        xi_d = (noeud_i)*H+(i-1)*h ;
        yi_g = rem(xi_g/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk_g = fix(yi_g*(Nbpt_h-1)) +1; %partie entière de y*(Nbpt-1) (indice pour les matrices)
        yi_d = rem(xi_d/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk_d = fix(yi_d*(Nbpt_h-1)) +1;

        %         vect_propre_eps(i) = Psi(yk_g);
        %         vect_propre_droite(i) = Psi(yk_d);

        t_g = yi_g*(Nbpt_h-1) - fix(yi_g*(Nbpt_h-1)) ;
        vect_propre_eps(i) = Psi(yk_g)*(1-t_g) +Psi(yk_g+1)*t_g;
        t_d = yi_d*(Nbpt_h-1) - fix(yi_d*(Nbpt_h-1)) ;
        vect_propre_droite(i) = Psi(yk_d)*(1-t_d) +Psi(yk_d+1)*t_d;
    end

    UU_g = zeros(Nbpt_h,1);
    UU_d = zeros(Nbpt_h,1);

    % boucle pour les matrices EF
    % ------------------------
    for i=1:(Nbpt_h)
        UU_g(i) = (i-1)*(h/H);
        UU_d(i) = (1-((i-1)*(h/H)));
    end % for l

        %Multiplication par Psi_eps
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)*vect_propre_eps(l);
        UU_d(l)=UU_d(l)*vect_propre_droite(l);

    end % for l

    %val = [ UU_g ; UU_d(2:Nbpt_h)];
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif num_idee==21

    MM_ref = matM_ref(Nbpt_h);
    Nbpt_per = 9*(Nbpt_h-1)+1-1;
    %Matrice EF
    MM = matM_per(Nbpt_per);
    MM_sigma = sparse(Nbpt_per,Nbpt_per); % matrice de masse avec sigma
    KK = sparse(Nbpt_per,Nbpt_per); % matrice de rigidite

    for l=2:(Nbpt_per-1)   %La première et dernière ligne doivent rester nulle

        KK(l,l)=(1/h)*(A((noeud_i-5)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-5)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-5)*H+(l-1+0.5)*h,epsilon); %#ok<*SPRIX>
        KK(l+1,l)=(-1/h)*A((noeud_i-5)*H+(l-1+0.5)*h,epsilon);

        MM_sigma(l,l)=(h/3)*(Sigma((noeud_i-5)*H+(l-1-0.5)*h,epsilon)+Sigma((noeud_i-5)*H+(l-1+0.5)*h,epsilon));
        MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i-5)*H+(l-1+0.5)*h,epsilon);
        MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i-5)*H+(l-1+0.5)*h,epsilon);

    end % for l

    KK(1,1)=(1/h)*(A((noeud_i-5)*H+0.5*h,epsilon)+A((noeud_i+4)*H-0.5*h,epsilon));
    KK(1,2)=(-1/h)*A((noeud_i-5)*H+0.5*h,epsilon);
    KK(2,1)=(-1/h)*A((noeud_i-5)*H+0.5*h,epsilon);
    KK(1,Nbpt_per)=(-1/h)*A((noeud_i+4)*H-0.5*h,epsilon);
    KK(Nbpt_per,1)=(-1/h)*A((noeud_i+4)*H-0.5*h,epsilon);
    KK(Nbpt_per,Nbpt_per)=(1/h)*(A((noeud_i-5)*H+(Nbpt_per-1-0.5)*h,epsilon)+A((noeud_i-5)*H+(Nbpt_per-1+0.5)*h,epsilon));

    MM_sigma(1,1)=(h/3)*(Sigma((noeud_i-5)*H+0.5*h,epsilon)+Sigma((noeud_i+4)*H-0.5*h,epsilon));
    MM_sigma(1,2)=h*(1/6)*Sigma((noeud_i-5)*H+0.5*h,epsilon);
    MM_sigma(2,1)=h*(1/6)*Sigma((noeud_i-5)*H+0.5*h,epsilon);
    MM_sigma(1,Nbpt_per)=h*(1/6)*Sigma((noeud_i+4)*H-0.5*h,epsilon);
    MM_sigma(Nbpt_per,1)=h*(1/6)*Sigma((noeud_i+4)*H-0.5*h,epsilon);
    MM_sigma(Nbpt_per,Nbpt_per)=(h/3)*(Sigma((noeud_i-5)*H+(Nbpt_per-1-0.5)*h,epsilon)+Sigma((noeud_i-5)*H+(Nbpt_per-1+0.5)*h,epsilon));


    AA = MM_sigma + KK*(epsilon^2) ;

    [EV,DV] = eigs(AA,MM,1,'smallestabs');

    [~,arg_min] = min(diag(DV));
    Vecteur_propre_H_periodic = EV(:,arg_min);

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre_H_periodic(fix(Nbpt_h/2))<0)
        Vecteur_propre_H_periodic = -Vecteur_propre_H_periodic;
    end

    Vecteur_propre_H_periodic = [Vecteur_propre_H_periodic;Vecteur_propre_H_periodic(1)];

    Vecteur_propre_H_periodic_mid = Vecteur_propre_H_periodic(4*(Nbpt_h-1)+1:5*(Nbpt_h-1)+1);

    % % Rescale pour une norme L^2 = 1
    Norme_L2 = sqrt(Vecteur_propre_H_periodic_mid'*MM_ref*Vecteur_propre_H_periodic_mid);
    vect_propre_eps = Vecteur_propre_H_periodic_mid/Norme_L2;
    % -------------------------------------------------------------------------

    MM_sigma = sparse(Nbpt_per,Nbpt_per); % matrice de masse avec sigma
    KK = sparse(Nbpt_per,Nbpt_per); % matrice de rigidite

    for l=2:(Nbpt_per-1)   %La première et dernière ligne doivent rester nulle

        KK(l,l)=(1/h)*(A((noeud_i-4)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-4)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-4)*H+(l-1+0.5)*h,epsilon); %#ok<*SPRIX>
        KK(l+1,l)=(-1/h)*A((noeud_i-4)*H+(l-1+0.5)*h,epsilon);

        MM_sigma(l,l)=(h/3)*(Sigma((noeud_i-4)*H+(l-1-0.5)*h,epsilon)+Sigma((noeud_i-4)*H+(l-1+0.5)*h,epsilon));
        MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i-4)*H+(l-1+0.5)*h,epsilon);
        MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i-4)*H+(l-1+0.5)*h,epsilon);

    end % for l

    KK(1,1)=(1/h)*(A((noeud_i-4)*H+0.5*h,epsilon)+A((noeud_i+5)*H-0.5*h,epsilon));
    KK(1,2)=(-1/h)*A((noeud_i-4)*H+0.5*h,epsilon);
    KK(2,1)=(-1/h)*A((noeud_i-4)*H+0.5*h,epsilon);
    KK(1,Nbpt_per)=(-1/h)*A((noeud_i+5)*H-0.5*h,epsilon);
    KK(Nbpt_per,1)=(-1/h)*A((noeud_i+5)*H-0.5*h,epsilon);
    KK(Nbpt_per,Nbpt_per)=(1/h)*(A((noeud_i-4)*H+(Nbpt_per-1-0.5)*h,epsilon)+A((noeud_i-4)*H+(Nbpt_per-1+0.5)*h,epsilon));

    MM_sigma(1,1)=(h/3)*(Sigma((noeud_i-4)*H+0.5*h,epsilon)+Sigma((noeud_i+5)*H-0.5*h,epsilon));
    MM_sigma(1,2)=h*(1/6)*Sigma((noeud_i-4)*H+0.5*h,epsilon);
    MM_sigma(2,1)=h*(1/6)*Sigma((noeud_i-4)*H+0.5*h,epsilon);
    MM_sigma(1,Nbpt_per)=h*(1/6)*Sigma((noeud_i+5)*H-0.5*h,epsilon);
    MM_sigma(Nbpt_per,1)=h*(1/6)*Sigma((noeud_i+5)*H-0.5*h,epsilon);
    MM_sigma(Nbpt_per,Nbpt_per)=(h/3)*(Sigma((noeud_i-4)*H+(Nbpt_per-1-0.5)*h,epsilon)+Sigma((noeud_i-4)*H+(Nbpt_per-1+0.5)*h,epsilon));


    AA = MM_sigma + KK*(epsilon^2) ;

    [EV,DV] = eigs(AA,MM,1,'smallestabs');

    [~,arg_min] = min(diag(DV));
    Vecteur_propre_H_periodic = EV(:,arg_min);

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre_H_periodic(fix(Nbpt_h/2))<0)
        Vecteur_propre_H_periodic = -Vecteur_propre_H_periodic;
    end

    Vecteur_propre_H_periodic = [Vecteur_propre_H_periodic;Vecteur_propre_H_periodic(1)];

    Vecteur_propre_H_periodic_mid = Vecteur_propre_H_periodic(4*(Nbpt_h-1)+1:5*(Nbpt_h-1)+1);

    % % Rescale pour une norme L^2 = 1
    Norme_L2 = sqrt(Vecteur_propre_H_periodic_mid'*MM_ref*Vecteur_propre_H_periodic_mid);
    vect_propre_droite = Vecteur_propre_H_periodic_mid/Norme_L2;
    % -------------------------------------------------------------------------

    % -------------------------------------------------------------------------

    KK_g = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_g = zeros(Nbpt_h,1);     % vecteur second membre
    KK_d = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_d = zeros(Nbpt_h,1);     % vecteur second membre

    % boucle pour les matrices EF
    % ------------------------
    for l=2:(Nbpt_h-2)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));

        KK_g(l,l)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_g(l,l+1)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK_g(l+1,l)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        KK_d(l,l)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_d(l,l+1)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);
        KK_d(l+1,l)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);

    end % for l
    phi_avant_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1-1));
    phi_apres_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1+1));
    phi_avant_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1-1));
    phi_apres_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1+1));

    KK_g(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
    KK_d(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)

    % Calcul du second membre F
    % -------------------------
    for l=2:(Nbpt_h-1)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));
        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        LL_g(l)=(-1/H)*(-1*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)+phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon));

        LL_d(l)=(1/H)*(-phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)+phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon));

    end % for l

    AA_g = KK_g ;
    AA_g(1,1)=1; AA_g(Nbpt_h,Nbpt_h)=1;

    AA_d = KK_d ;
    AA_d(1,1)=1; AA_d(Nbpt_h,Nbpt_h)=1;
    % inversion
    % ----------

    UU_g = AA_g\LL_g;
    UU_d = AA_d\LL_d;

    % réhaussement
    % ----------
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)+(l-1)*(h/H);
        UU_d(l)=UU_d(l)-(l-1)*(h/H)+1;

    end % for l

    %Multiplication par Psi_eps
    for l=1:Nbpt_h

%         UU_g(l)=UU_g(l)*vect_propre_eps(l)/vect_propre_eps(Nbpt_h);
%         UU_d(l)=UU_d(l)*vect_propre_droite(l)/vect_propre_droite(1);
                UU_g(l)=UU_g(l)*vect_propre_eps(l);
                UU_d(l)=UU_d(l)*vect_propre_droite(l);

    end % for l
    %val = [ UU_g ; UU_d(2:Nbpt_h)];
    %     coeff = UU_g(Nbpt_h);
    %     val = val./coeff ;
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif num_idee==22
    bool_moyennisation = 1;
    Largeur_dirac = 0.05; %% C'est un pourcentage de epsilon (0.1 := le dirac sera 1/10 de la taille de epsilon)
    X_0 = (noeud_i-2)*H;
    taille_domaine = 3*H;
    Sigma_epsXX = zeros((3*(Nbpt_h-1)+1),1);     %vecteur des coordonnées des noeuds
    Sigma_yXX = zeros(Nbpt_ref,1);     %vecteur des coordonnées des noeuds
    for j=1:(3*(Nbpt_h-1)+1)
        Sigma_epsXX(j)=Sigma((noeud_i-2)*H+(j-0.5)*h,epsilon);
    end
    Moy_sigma_eps = mean(Sigma_epsXX);
    %Sigma_yXXmean = mean(Sigma_yXX);

    for j=1:Nbpt_ref
    Sigma_yXX(j)=Sigma((j-1)*h,epsilon/1); % On calcul la moyenne de sigma comme la moyenne sur Omega
    end
    Moy_sigma_y = mean(Sigma_yXX);

    Sigma_moyXX = zeros((3*(Nbpt_h-1)+1),1);     %vecteur des coordonnées des noeuds
    for j=1:(3*(Nbpt_h-1)+1)
        Sigma_moyXX(j)=Sigma_moy_nulle((noeud_i-2)*H+(j-0.5)*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation);
    end
    Sigma_moyXXmean = mean(Sigma_moyXX);

    Pourcentage_ecart_moyenne_sigma = ((Moy_sigma_y - Moy_sigma_eps)/Moy_sigma_y)*100

    MM_ref = matM_ref(Nbpt_h);
    Nbpt_per = 3*(Nbpt_h-1)+1-1;
    %Matrice EF
    MM = matM_per(Nbpt_per);
    MM_sigma = sparse(Nbpt_per,Nbpt_per); % matrice de masse avec sigma
    KK = sparse(Nbpt_per,Nbpt_per); % matrice de rigidite

    for l=2:(Nbpt_per-1)   %La première et dernière ligne doivent rester nulle

        KK(l,l)=(1/h)*(A((noeud_i-2)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-2)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-2)*H+(l-1+0.5)*h,epsilon); %#ok<*SPRIX>
        KK(l+1,l)=(-1/h)*A((noeud_i-2)*H+(l-1+0.5)*h,epsilon);

        MM_sigma(l,l)=(h/3)*(Sigma_moy_nulle((noeud_i-2)*H+(l-1-0.5)*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation) ...
            +Sigma_moy_nulle((noeud_i-2)*H+(l-1+0.5)*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation));
        MM_sigma(l,l+1)=h*(1/6)*Sigma_moy_nulle((noeud_i-2)*H+(l-1+0.5)*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation);
        MM_sigma(l+1,l)=h*(1/6)*Sigma_moy_nulle((noeud_i-2)*H+(l-1+0.5)*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation);

    end % for l

    KK(1,1)=(1/h)*(A((noeud_i-2)*H+0.5*h,epsilon)+A((noeud_i+1)*H-0.5*h,epsilon));
    KK(1,2)=(-1/h)*A((noeud_i-2)*H+0.5*h,epsilon);
    KK(2,1)=(-1/h)*A((noeud_i-2)*H+0.5*h,epsilon);
    KK(1,Nbpt_per)=(-1/h)*A((noeud_i+1)*H-0.5*h,epsilon);
    KK(Nbpt_per,1)=(-1/h)*A((noeud_i+1)*H-0.5*h,epsilon);
    KK(Nbpt_per,Nbpt_per)=(1/h)*(A((noeud_i-2)*H+(Nbpt_per-1-0.5)*h,epsilon)+A((noeud_i-2)*H+(Nbpt_per-1+0.5)*h,epsilon));

    MM_sigma(1,1)=(h/3)*(Sigma_moy_nulle((noeud_i-2)*H+0.5*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation) ...
        +Sigma_moy_nulle((noeud_i+1)*H-0.5*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation));
    MM_sigma(1,2)=h*(1/6)*Sigma_moy_nulle((noeud_i-2)*H+0.5*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation);
    MM_sigma(2,1)=h*(1/6)*Sigma_moy_nulle((noeud_i-2)*H+0.5*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation);
    MM_sigma(1,Nbpt_per)=h*(1/6)*Sigma_moy_nulle((noeud_i+1)*H-0.5*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation);
    MM_sigma(Nbpt_per,1)=h*(1/6)*Sigma_moy_nulle((noeud_i+1)*H-0.5*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation);
    MM_sigma(Nbpt_per,Nbpt_per)=(h/3)*(Sigma_moy_nulle((noeud_i-2)*H+(Nbpt_per-1-0.5)*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation) ...
        +Sigma_moy_nulle((noeud_i-2)*H+(Nbpt_per-1+0.5)*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation));


    AA = MM_sigma + KK*(epsilon^2) ;

    [EV,DV] = eigs(AA,MM,1,'smallestabs');

    [~,arg_min] = min(diag(DV));
    Vecteur_propre_H_periodic = EV(:,arg_min);

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre_H_periodic(fix(Nbpt_h/2))<0)
        Vecteur_propre_H_periodic = -Vecteur_propre_H_periodic;
    end

    Vecteur_propre_H_periodic = [Vecteur_propre_H_periodic;Vecteur_propre_H_periodic(1)];

    Vecteur_propre_H_periodic_mid = Vecteur_propre_H_periodic(Nbpt_h:2*Nbpt_h-1);

    % % Rescale pour une norme L^2 = 1
    Norme_L2 = sqrt(Vecteur_propre_H_periodic_mid'*MM_ref*Vecteur_propre_H_periodic_mid);
    vect_propre_eps = Vecteur_propre_H_periodic_mid/Norme_L2;
    % -------------------------------------------------------------------------

    X_0 = (noeud_i-1)*H;
    
    Sigma_epsXX = zeros((3*(Nbpt_h-1)+1),1);     %vecteur des coordonnées des noeuds
    for j=1:(3*(Nbpt_h-1)+1)
        Sigma_epsXX(j)=Sigma((noeud_i-1)*H+(j-0.5)*h,epsilon);
    end

    Moy_sigma_eps = mean(Sigma_epsXX);

    MM_sigma = sparse(Nbpt_per,Nbpt_per); % matrice de masse avec sigma
    KK = sparse(Nbpt_per,Nbpt_per); % matrice de rigidite

    for l=2:(Nbpt_per-1)   %La première et dernière ligne doivent rester nulle

        KK(l,l)=(1/h)*(A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon); %#ok<*SPRIX>
        KK(l+1,l)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        MM_sigma(l,l)=(h/3)*(Sigma_moy_nulle((noeud_i-1)*H+(l-1-0.5)*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation) ...
            +Sigma_moy_nulle((noeud_i-1)*H+(l-1+0.5)*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation));
        MM_sigma(l,l+1)=h*(1/6)*Sigma_moy_nulle((noeud_i-1)*H+(l-1+0.5)*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation);
        MM_sigma(l+1,l)=h*(1/6)*Sigma_moy_nulle((noeud_i-1)*H+(l-1+0.5)*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation);

    end % for l

    KK(1,1)=(1/h)*(A((noeud_i-1)*H+0.5*h,epsilon)+A((noeud_i+2)*H-0.5*h,epsilon));
    KK(1,2)=(-1/h)*A((noeud_i-1)*H+0.5*h,epsilon);
    KK(2,1)=(-1/h)*A((noeud_i-1)*H+0.5*h,epsilon);
    KK(1,Nbpt_per)=(-1/h)*A((noeud_i+2)*H-0.5*h,epsilon);
    KK(Nbpt_per,1)=(-1/h)*A((noeud_i+2)*H-0.5*h,epsilon);
    KK(Nbpt_per,Nbpt_per)=(1/h)*(A((noeud_i-1)*H+(Nbpt_per-1-0.5)*h,epsilon)+A((noeud_i-1)*H+(Nbpt_per-1+0.5)*h,epsilon));

    MM_sigma(1,1)=(h/3)*(Sigma_moy_nulle((noeud_i-1)*H+0.5*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation) ...
        +Sigma_moy_nulle((noeud_i+2)*H-0.5*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation));
    MM_sigma(1,2)=h*(1/6)*Sigma_moy_nulle((noeud_i-1)*H+0.5*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation);
    MM_sigma(2,1)=h*(1/6)*Sigma_moy_nulle((noeud_i-1)*H+0.5*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation);
    MM_sigma(1,Nbpt_per)=h*(1/6)*Sigma_moy_nulle((noeud_i+2)*H-0.5*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation);
    MM_sigma(Nbpt_per,1)=h*(1/6)*Sigma_moy_nulle((noeud_i+2)*H-0.5*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation);
    MM_sigma(Nbpt_per,Nbpt_per)=(h/3)*(Sigma_moy_nulle((noeud_i-1)*H+(Nbpt_per-1-0.5)*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation) ...
        +Sigma_moy_nulle((noeud_i-1)*H+(Nbpt_per-1+0.5)*h,epsilon,Moy_sigma_eps,Moy_sigma_y,Largeur_dirac,X_0, taille_domaine,bool_moyennisation));


    AA = MM_sigma + KK*(epsilon^2) ;

    [EV,DV] = eigs(AA,MM,1,'smallestabs');

    [~,arg_min] = min(diag(DV));
    Vecteur_propre_H_periodic = EV(:,arg_min);

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre_H_periodic(fix(Nbpt_h/2))<0)
        Vecteur_propre_H_periodic = -Vecteur_propre_H_periodic;
    end

    Vecteur_propre_H_periodic = [Vecteur_propre_H_periodic;Vecteur_propre_H_periodic(1)];

    Vecteur_propre_H_periodic_mid = Vecteur_propre_H_periodic(Nbpt_h:2*Nbpt_h-1);

    % % Rescale pour une norme L^2 = 1
    Norme_L2 = sqrt(Vecteur_propre_H_periodic_mid'*MM_ref*Vecteur_propre_H_periodic_mid);
    vect_propre_droite = Vecteur_propre_H_periodic_mid/Norme_L2;
    % -------------------------------------------------------------------------

    % -------------------------------------------------------------------------

    KK_g = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_g = zeros(Nbpt_h,1);     % vecteur second membre
    KK_d = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_d = zeros(Nbpt_h,1);     % vecteur second membre

    % boucle pour les matrices EF
    % ------------------------
    for l=2:(Nbpt_h-2)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));

        KK_g(l,l)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_g(l,l+1)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK_g(l+1,l)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        KK_d(l,l)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_d(l,l+1)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);
        KK_d(l+1,l)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);

    end % for l
    phi_avant_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1-1));
    phi_apres_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1+1));
    phi_avant_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1-1));
    phi_apres_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1+1));

    KK_g(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
    KK_d(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)

    % Calcul du second membre F
    % -------------------------
    for l=2:(Nbpt_h-1)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));
        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        LL_g(l)=(-1/H)*(-1*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)+phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon));

        LL_d(l)=(1/H)*(-phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)+phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon));

    end % for l

    AA_g = KK_g ;
    AA_g(1,1)=1; AA_g(Nbpt_h,Nbpt_h)=1;

    AA_d = KK_d ;
    AA_d(1,1)=1; AA_d(Nbpt_h,Nbpt_h)=1;
    % inversion
    % ----------

    UU_g = AA_g\LL_g;
    UU_d = AA_d\LL_d;

    % réhaussement
    % ----------
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)+(l-1)*(h/H);
        UU_d(l)=UU_d(l)-(l-1)*(h/H)+1;

    end % for l

    %Multiplication par Psi_eps
    for l=1:Nbpt_h

%         UU_g(l)=UU_g(l)*vect_propre_eps(l)/vect_propre_eps(Nbpt_h);
%         UU_d(l)=UU_d(l)*vect_propre_droite(l)/vect_propre_droite(1);
                UU_g(l)=UU_g(l)*vect_propre_eps(l);
                UU_d(l)=UU_d(l)*vect_propre_droite(l);

    end % for l
    %val = [ UU_g ; UU_d(2:Nbpt_h)];
    %     coeff = UU_g(Nbpt_h);
    %     val = val./coeff ;
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
elseif num_idee==23

    MM_ref = matM_ref(Nbpt_h);
    Nbpt_OS = 3*(Nbpt_h-1)+1; % Nombre de points dans le patch d'OS
    %Matrice EF
    MM = sparse(Nbpt_OS,Nbpt_OS); % matrice de masse
    MM_sigma = sparse(Nbpt_OS,Nbpt_OS); % matrice de masse avec sigma
    KK = sparse(Nbpt_OS,Nbpt_OS); % matrice de rigidite
    R = 3*H;
    %On construit les matrice EF avec le filtre
    for l=1:(Nbpt_OS-1)   %La première et dernière ligne doivent rester nulle

        MM(l,l)=h*(1/3)*((1/R)*Phi((1/R)*(l-1-0.5)*h) + (1/R)*Phi((1/R)*(l-1+0.5)*h));
        MM(l,l+1)=h*(1/6)*(1/R)*Phi((1/R)*(l-1+0.5)*h);
        MM(l+1,l)=h*(1/6)*(1/R)*Phi((1/R)*(l-1+0.5)*h);

        KK(l,l)=(1/h)*(A((noeud_i-2)*H+(l-1-0.5)*h,epsilon)*(1/R)*Phi((1/R)*(l-1-0.5)*h)+A((noeud_i-2)*H+(l-1+0.5)*h,epsilon)*(1/R)*Phi((1/R)*(l-1+0.5)*h)  ); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-2)*H+(l-1+0.5)*h,epsilon)*(1/R)*Phi((1/R)*(l-1+0.5)*h); %#ok<*SPRIX>
        KK(l+1,l)=(-1/h)*A((noeud_i-2)*H+(l-1+0.5)*h,epsilon)*(1/R)*Phi((1/R)*(l-1+0.5)*h);

        MM_sigma(l,l)=(h/3)*(Sigma((noeud_i-2)*H+(l-1-0.5)*h,epsilon)*(1/R)*Phi((1/R)*(l-1-0.5)*h)+Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon)*(1/R)*Phi((1/R)*(l-1+0.5)*h) );
        MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon)*(1/R)*Phi((1/R)*(l-1+0.5)*h);
        MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i-2)*H+(l-1+0.5)*h,epsilon)*(1/R)*Phi((1/R)*(l-1+0.5)*h);

    end % for l

    MM(1,1)=h*(1/3)*(1/R)*Phi((1/R)*0.5*h);
    MM(Nbpt_OS,Nbpt_OS)=h*(1/3)*(1/R)*Phi((1/R)*(Nbpt_OS-1-0.5)*h);
    
    KK(1,1)=(1/h)*A((noeud_i-2)*H+0.5*h,epsilon)*(1/R)*Phi((1/R)*0.5*h);
    KK(Nbpt_OS,Nbpt_OS)=(1/h)*A((noeud_i+1)*H-0.5*h,epsilon)*(1/R)*Phi((1/R)*(Nbpt_OS-1-0.5)*h);

    MM_sigma(1,1)=(h/3)*Sigma((noeud_i-2)*H+0.5*h,epsilon)*(1/R)*Phi((1/R)*0.5*h);
    MM_sigma(Nbpt_OS,Nbpt_OS)=(h/3)*Sigma((noeud_i+1)*H-0.5*h,epsilon)*(1/R)*Phi((1/R)*(Nbpt_OS-1-0.5)*h);

    AA = MM_sigma + KK*(epsilon^2) ;

    BB = zeros(Nbpt_OS,1); %Matrice contrainte périodicité
    
    for l=2:(Nbpt_OS-1)
        BB(l)= (1/R)*Phi((1/R)*(l-1-0.5)*h)-(1/R)*Phi((1/R)*(l-1+0.5)*h);
    end
    BB(1)=-1*(1/R)*Phi((1/R)*(0.5)*h);
    BB(Nbpt_OS)=1*(1/R)*Phi((1/R)*(Nbpt_OS-1-0.5)*h);

    AA_f = [AA , BB ; BB' , 0];
    
    MM_f = [MM , 0*BB ; 0*BB' , 0] ;

    [EV,DV] = eigs(AA_f,MM_f,1,'smallestabs');

    [~,arg_min] = min(diag(DV));
    Vecteur_propre_H_periodic = EV(:,arg_min);

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre_H_periodic(fix(Nbpt_h/2))<0)
        Vecteur_propre_H_periodic = -Vecteur_propre_H_periodic;
    end

    Vecteur_propre_H_periodic = Vecteur_propre_H_periodic(1:Nbpt_OS);
    Vecteur_propre_H_periodic_mid = Vecteur_propre_H_periodic(Nbpt_h:2*Nbpt_h-1);

    % % Rescale pour une norme L^2 = 1
    Norme_L2 = sqrt(Vecteur_propre_H_periodic_mid'*MM_ref*Vecteur_propre_H_periodic_mid);
    vect_propre_eps = Vecteur_propre_H_periodic_mid/Norme_L2;
    % -------------------------------------------------------------------------

    MM_sigma = sparse(Nbpt_OS,Nbpt_OS); % matrice de masse avec sigma
    KK = sparse(Nbpt_OS,Nbpt_OS); % matrice de rigidite

    for l=1:(Nbpt_OS-1)   %La première et dernière ligne doivent rester nulle

        KK(l,l)=(1/h)*(A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)*(1/R)*Phi((1/R)*(l-1-0.5)*h)+A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)*(1/R)*Phi((1/R)*(l-1+0.5)*h)  ); %approximation de l'intégrale de a(.)
        KK(l,l+1)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)*(1/R)*Phi((1/R)*(l-1+0.5)*h); %#ok<*SPRIX>
        KK(l+1,l)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)*(1/R)*Phi((1/R)*(l-1+0.5)*h);

        MM_sigma(l,l)=(h/3)*(Sigma((noeud_i-1)*H+(l-1-0.5)*h,epsilon)*(1/R)*Phi((1/R)*(l-1-0.5)*h)+Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon)*(1/R)*Phi((1/R)*(l-1+0.5)*h) );
        MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon)*(1/R)*Phi((1/R)*(l-1+0.5)*h);
        MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon)*(1/R)*Phi((1/R)*(l-1+0.5)*h);

    end % for l

    KK(1,1)=(1/h)*A((noeud_i-1)*H+0.5*h,epsilon)*(1/R)*Phi((1/R)*0.5*h);
    KK(Nbpt_OS,Nbpt_OS)=(1/h)*A((noeud_i+2)*H-0.5*h,epsilon)*(1/R)*Phi((1/R)*(Nbpt_OS-1-0.5)*h);

    MM_sigma(1,1)=(h/3)*Sigma((noeud_i-1)*H+0.5*h,epsilon)*(1/R)*Phi((1/R)*0.5*h);
    MM_sigma(Nbpt_OS,Nbpt_OS)=(h/3)*Sigma((noeud_i+2)*H-0.5*h,epsilon)*(1/R)*Phi((1/R)*(Nbpt_OS-1-0.5)*h);

    AA = MM_sigma + KK*(epsilon^2) ;

    AA_f = [AA , BB ; BB' , 0];
    
    MM_f = [MM , 0*BB ; 0*BB' , 0] ;

    [EV,DV] = eigs(AA_f,MM_f,1,'smallestabs');

    [~,arg_min] = min(diag(DV));
    Vecteur_propre_H_periodic = EV(:,arg_min);

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre_H_periodic(fix(Nbpt_h/2))<0)
        Vecteur_propre_H_periodic = -Vecteur_propre_H_periodic;
    end

    Vecteur_propre_H_periodic = Vecteur_propre_H_periodic(1:Nbpt_OS);
    Vecteur_propre_H_periodic_mid = Vecteur_propre_H_periodic(Nbpt_h:2*Nbpt_h-1);

    % % Rescale pour une norme L^2 = 1
    Norme_L2 = sqrt(Vecteur_propre_H_periodic_mid'*MM_ref*Vecteur_propre_H_periodic_mid);
    vect_propre_droite = Vecteur_propre_H_periodic_mid/Norme_L2;
    % -------------------------------------------------------------------------

    % -------------------------------------------------------------------------

    KK_g = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_g = zeros(Nbpt_h,1);     % vecteur second membre
    KK_d = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
    LL_d = zeros(Nbpt_h,1);     % vecteur second membre

    % boucle pour les matrices EF
    % ------------------------
    for l=2:(Nbpt_h-2)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));

        KK_g(l,l)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_g(l,l+1)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
        KK_g(l+1,l)=(-1/h)*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        KK_d(l,l)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
        KK_d(l,l+1)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);
        KK_d(l+1,l)=(-1/h)*phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon);

    end % for l
    phi_avant_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1-1));
    phi_apres_d = 0.5*(vect_propre_droite(Nbpt_h-1)+vect_propre_droite(Nbpt_h-1+1));
    phi_avant_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1-1));
    phi_apres_g = 0.5*(vect_propre_eps(Nbpt_h-1)+vect_propre_eps(Nbpt_h-1+1));

    KK_g(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
    KK_d(Nbpt_h-1,Nbpt_h-1)=(1/h)*(phi_avant_d*phi_avant_d*A((noeud_i)*H+(Nbpt_h-1-1-0.5)*h,epsilon)  +  phi_apres_d*phi_apres_d*A((noeud_i)*H+(Nbpt_h-1-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)

    % Calcul du second membre F
    % -------------------------
    for l=2:(Nbpt_h-1)   %La première et dernière ligne doivent rester nulle

        phi_avant_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l-1));
        phi_apres_g = 0.5*(vect_propre_eps(l)+vect_propre_eps(l+1));
        phi_avant_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l-1));
        phi_apres_d = 0.5*(vect_propre_droite(l)+vect_propre_droite(l+1));

        LL_g(l)=(-1/H)*(-1*phi_apres_g*phi_apres_g*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)+phi_avant_g*phi_avant_g*A((noeud_i-1)*H+(l-1-0.5)*h,epsilon));

        LL_d(l)=(1/H)*(-phi_apres_d*phi_apres_d*A((noeud_i)*H+(l-1+0.5)*h,epsilon)+phi_avant_d*phi_avant_d*A((noeud_i)*H+(l-1-0.5)*h,epsilon));

    end % for l

    AA_g = KK_g ;
    AA_g(1,1)=1; AA_g(Nbpt_h,Nbpt_h)=1;

    AA_d = KK_d ;
    AA_d(1,1)=1; AA_d(Nbpt_h,Nbpt_h)=1;
    % inversion
    % ----------

    UU_g = AA_g\LL_g;
    UU_d = AA_d\LL_d;

    % réhaussement
    % ----------
    for l=1:Nbpt_h

        UU_g(l)=UU_g(l)+(l-1)*(h/H);
        UU_d(l)=UU_d(l)-(l-1)*(h/H)+1;

    end % for l

    %Multiplication par Psi_eps
    for l=1:Nbpt_h

%         UU_g(l)=UU_g(l)*vect_propre_eps(l)/vect_propre_eps(Nbpt_h);
%         UU_d(l)=UU_d(l)*vect_propre_droite(l)/vect_propre_droite(1);
                UU_g(l)=UU_g(l)*vect_propre_eps(l);
                UU_d(l)=UU_d(l)*vect_propre_droite(l);

    end % for l
    %val = [ UU_g ; UU_d(2:Nbpt_h)];
    %     coeff = UU_g(Nbpt_h);
    %     val = val./coeff ;
    % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
