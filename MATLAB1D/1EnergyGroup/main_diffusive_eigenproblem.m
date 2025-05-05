% =====================================================
%
% principal_reaction_diffusion;
%
% une routine qui calcul pour un pas H donné l'erreur entre la solution approchée par la méthode de résolution de notre choix
% et la solution exacte, calculée avec les éléments P1 classiques de pas de
% maillage très petit


% ---------------------------
%Paramètres
% ---------------------------
epsilon=0.05166;
Nbpt_H=11;   %Maillage grossier
Nbpt_ref = 10001;  %Maillage fin

num_idee = 23 ; %Voir la description dans la fonction Khi_i() pour les détails des numéros
% Num 23 est la methode qui marche in fine (methode du filtre)

affichage_erreurs = 1;
%Affichages des erreurs entre solution approchee et solution de référence :
% 0 Si pas d'affichage
% 1 pour affichage

trace_solution_Tf = 1;
%Tracé des solutions exactes et approchées :
% 0 Si pas de tracé
% 1 pour le tracé solution exacte et solution MsFEM

% ---------------------------


h_ref = 1/(Nbpt_ref -1);
H = 1/(Nbpt_H -1);

% Maillage
% ----------------------
XX_ref = zeros(Nbpt_ref,1);     %vecteur des coordonnées des noeuds du maillage de pas H
for i=1:Nbpt_ref
    XX_ref(i)=(i-1)*h_ref;
end
XX_H = zeros(Nbpt_H,1);     %vecteur des coordonnées des noeuds du maillage de pas H
for i=1:Nbpt_H
    XX_H(i)=(i-1)*H;
end
tic;

[val_propre_ref,Vecteur_propre_ref] = fonction_propre_eps_Dirichlet(2,Nbpt_ref,epsilon,1);

tempsEcoule = toc;
disp(['Le temps de calcul de la solution de reference est de ', num2str(tempsEcoule), ' secondes.']);
% calcul des matrices EF pour solution approchee
% ----------------------------------------------
tic;

Nbpt_h = (Nbpt_ref - 1)/(Nbpt_H-1) +1;
Khi_global = zeros(Nbpt_h,Nbpt_H-2,2);
for l=1:(Nbpt_H-2)
    [fct_forme_gauche , fct_forme_droite] = Khi_i(Nbpt_H,Nbpt_h,epsilon,l,num_idee);

    Khi_global(:,l,1)= fct_forme_gauche ;
    Khi_global(:,l,2)= fct_forme_droite ;
end
MM=matM(Khi_global,Nbpt_H,Nbpt_h);
MM_sigma=matM_sigma(Khi_global,Nbpt_H,Nbpt_h,epsilon);
KK=matK(Khi_global,Nbpt_H,Nbpt_h,epsilon);


AA = MM_sigma + KK*(epsilon^2) ;
AA(1,1)=1; AA(Nbpt_H,Nbpt_H)=1;

%[EV,DV] = eigs(AA,MM,5,'smallestabs');
[EV,DV] = eigs(AA,MM,1,'smallestabs');


[Valeur_propre_min,arg_min] = min(diag(DV));
Vecteur_propre_H = EV(:,arg_min);

%On a un VP de norme L^2 unitaire. On choisis le VP positif.
if (Vecteur_propre_H(fix(Nbpt_H/2))<0)
    Vecteur_propre_H = -Vecteur_propre_H;
end

% % Rescale pour une norme L^2 = 1
Norme_L2 = sqrt(Vecteur_propre_H'*MM*Vecteur_propre_H);
Vecteur_propre_H = Vecteur_propre_H/Norme_L2;

tempsEcoule = toc;
disp(['Le temps de calcul de la solution approchee est de ', num2str(tempsEcoule), ' secondes.']);

if trace_solution_Tf ==1
    figure;
    plot(XX_ref,Vecteur_propre_ref,'o')
    hold on ;

    for i=2:(Nbpt_H-2)
        % On trace la solution entre (i-1)H et iH
        solution_intervalle_i = zeros(Nbpt_ref,1);
        solution_intervalle_i = solution_intervalle_i*NaN;
        for j=1:Nbpt_h
            solution_intervalle_i((i-1)*(Nbpt_h-1)+j) = Vecteur_propre_H(i)*Khi_global(j,i-1,2) + Vecteur_propre_H(i+1)*Khi_global(j,i,1);
        end
        plot(XX_ref,solution_intervalle_i,LineWidth=2,color='r')
        hold on;
    end

    %Premier intervalle
    solution_intervalle_i = zeros(Nbpt_ref,1);
    solution_intervalle_i = solution_intervalle_i*NaN;
    i = 1; % On trace la solution entre 0 et H
    for j=1:Nbpt_h
        solution_intervalle_i((i-1)*(Nbpt_h-1)+j) = Vecteur_propre_H(i+1)*Khi_global(j,i,1);
    end
    plot(XX_ref,solution_intervalle_i,LineWidth=2,color='r')

    %Dernier intervalle
    solution_intervalle_i = zeros(Nbpt_ref,1);
    solution_intervalle_i = solution_intervalle_i*NaN;
    i = Nbpt_H-1; % On trace la solution entre 1-H et 1
    for j=1:Nbpt_h
        solution_intervalle_i((i-1)*(Nbpt_h-1)+j) = Vecteur_propre_H(i)*Khi_global(j,i-1,2);
    end
    plot(XX_ref,solution_intervalle_i,LineWidth=2,color='r')

    axis([0,1, 0 , 1.2*max(Vecteur_propre_ref)]);
    legend('solution exacte','solution approchée')
    title(sprintf("eps=%g / H=%g",epsilon,H))
    set(gca,'FontSize',30)
    hold off;
end


% -------------------------------------------------------------------
% Calcul de toutes les erreurs
if affichage_erreurs==1
    Erreur_H1 = 0;
    KK_ref_A_constant = matK_unitaire(Nbpt_h,h_ref);
    Norme_H1_vect_propre_ref = 0;

    for i=2:(Nbpt_H-2)
        Vecteur_propre_ref_local = Vecteur_propre_ref((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h);
        solution_intervalle_i = zeros(Nbpt_h,1);
        for j=1:Nbpt_h
            solution_intervalle_i(j) = Vecteur_propre_H(i)*Khi_global(j,i-1,2) + Vecteur_propre_H(i+1)*Khi_global(j,i,1);
        end
        Norme_H1_vect_propre_ref = Norme_H1_vect_propre_ref + (Vecteur_propre_ref_local'*KK_ref_A_constant*Vecteur_propre_ref_local);
        Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant*(solution_intervalle_i-Vecteur_propre_ref_local));
        Erreur_H1 = Erreur_H1 + Erreur_H1_locale;
    end

    %Premier intervalle
    i = 1; % On trace la solution entre 0 et H
    Vecteur_propre_ref_local = Vecteur_propre_ref((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h);
    solution_intervalle_i = zeros(Nbpt_h,1);
    for j=1:Nbpt_h
        solution_intervalle_i(j) = Vecteur_propre_H(i+1)*Khi_global(j,i,1);
    end
    Norme_H1_vect_propre_ref = Norme_H1_vect_propre_ref + (Vecteur_propre_ref_local'*KK_ref_A_constant*Vecteur_propre_ref_local);
    Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant*(solution_intervalle_i-Vecteur_propre_ref_local));
    Erreur_H1 = Erreur_H1 + Erreur_H1_locale;

    %Dernier intervalle
    i = Nbpt_H-1; % On trace la solution entre 1-H et 1
    Vecteur_propre_ref_local = Vecteur_propre_ref((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h);
    solution_intervalle_i = zeros(Nbpt_h,1);
    for j=1:Nbpt_h
        solution_intervalle_i(j) = Vecteur_propre_H(i)*Khi_global(j,i-1,2);
    end
    Norme_H1_vect_propre_ref = Norme_H1_vect_propre_ref + (Vecteur_propre_ref_local'*KK_ref_A_constant*Vecteur_propre_ref_local);
    Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant*(solution_intervalle_i-Vecteur_propre_ref_local));
    Erreur_H1 = Erreur_H1 + Erreur_H1_locale;

    Erreur_H1_relative = Erreur_H1/Norme_H1_vect_propre_ref;
    Erreur_H1_relative = sqrt(Erreur_H1_relative) %#ok<NOPTS>




    [~,Psi] = Solution_pb_spectral(Nbpt_ref);
    Psi_eps = zeros(Nbpt_ref,1);

    for i=1:Nbpt_ref
        xi_g = (i-1)*h_ref ;
        yi_g = rem(xi_g/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk_g = fix(yi_g*(Nbpt_ref-1)) +1; %partie entière de y*(Nbpt-1) (indice pour les matrices)

        t_g = yi_g*(Nbpt_h-1) - fix(yi_g*(Nbpt_h-1)) ;
        Psi_eps(i) = Psi(yk_g)*(1-t_g) +Psi(yk_g+1)*t_g;
    end

    Vecteur_propre_ref_sur_Psi = Vecteur_propre_ref./Psi_eps;

    Erreur_H1 = 0;
    Norme_H1_vect_propre_ref_sur_Psi = 0;

    for i=2:(Nbpt_H-2)
        Vecteur_propre_ref_local = Vecteur_propre_ref_sur_Psi((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h);
        solution_intervalle_i = zeros(Nbpt_h,1);
        for j=1:Nbpt_h

            xi_g = (j-1)*h_ref+(i-1)*H ;
            yi_g = rem(xi_g/epsilon,1) ; % y vaut x/epsilon modulo 1
            yk_g = fix(yi_g*(Nbpt_ref-1)) +1; %partie entière de y*(Nbpt-1) (indice pour les matrices)

            t_g = yi_g*(Nbpt_h-1) - fix(yi_g*(Nbpt_h-1)) ;
            Psi_eps_j = Psi(yk_g)*(1-t_g) +Psi(yk_g+1)*t_g;
            solution_intervalle_i(j) = (Vecteur_propre_H(i)*Khi_global(j,i-1,2) + Vecteur_propre_H(i+1)*Khi_global(j,i,1))/Psi_eps_j;
        end
        Norme_H1_vect_propre_ref_sur_Psi = Norme_H1_vect_propre_ref_sur_Psi + (Vecteur_propre_ref_local'*KK_ref_A_constant*Vecteur_propre_ref_local);
        Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant*(solution_intervalle_i-Vecteur_propre_ref_local));
        Erreur_H1 = Erreur_H1 + Erreur_H1_locale;
    end

    %Premier intervalle
    i = 1; % On trace la solution entre 0 et H
    Vecteur_propre_ref_local = Vecteur_propre_ref_sur_Psi((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h);
    solution_intervalle_i = zeros(Nbpt_h,1);
    for j=1:Nbpt_h

        xi_g = (j-1)*h_ref ;
        yi_g = rem(xi_g/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk_g = fix(yi_g*(Nbpt_ref-1)) +1; %partie entière de y*(Nbpt-1) (indice pour les matrices)

        t_g = yi_g*(Nbpt_h-1) - fix(yi_g*(Nbpt_h-1)) ;
        Psi_eps_j = Psi(yk_g)*(1-t_g) +Psi(yk_g+1)*t_g;

        solution_intervalle_i(j) = Vecteur_propre_H(i+1)*Khi_global(j,i,1)/Psi_eps_j;
    end
    Norme_H1_vect_propre_ref_sur_Psi = Norme_H1_vect_propre_ref_sur_Psi + (Vecteur_propre_ref_local'*KK_ref_A_constant*Vecteur_propre_ref_local);
    Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant*(solution_intervalle_i-Vecteur_propre_ref_local));
    Erreur_H1 = Erreur_H1 + Erreur_H1_locale;

    %Dernier intervalle
    i = Nbpt_H-1; % On trace la solution entre 1-H et 1
    Vecteur_propre_ref_local = Vecteur_propre_ref_sur_Psi((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h);
    solution_intervalle_i = zeros(Nbpt_h,1);
    for j=1:Nbpt_h
        xi_g = (j-1)*h_ref+(i-1)*H ;
        yi_g = rem(xi_g/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk_g = fix(yi_g*(Nbpt_ref-1)) +1; %partie entière de y*(Nbpt-1) (indice pour les matrices)

        t_g = yi_g*(Nbpt_h-1) - fix(yi_g*(Nbpt_h-1)) ;
        Psi_eps_j = Psi(yk_g)*(1-t_g) +Psi(yk_g+1)*t_g;

        solution_intervalle_i(j) = Vecteur_propre_H(i)*Khi_global(j,i-1,2)/Psi_eps_j;
    end
    Norme_H1_vect_propre_ref_sur_Psi = Norme_H1_vect_propre_ref_sur_Psi + (Vecteur_propre_ref_local'*KK_ref_A_constant*Vecteur_propre_ref_local);
    Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant*(solution_intervalle_i-Vecteur_propre_ref_local));
    Erreur_H1 = Erreur_H1 + Erreur_H1_locale;

    Erreur_H1_relative_sur_Psi = Erreur_H1/Norme_H1_vect_propre_ref_sur_Psi;
    Erreur_H1_relative_sur_Psi = sqrt(Erreur_H1_relative_sur_Psi) %#ok<NOPTS>

end
%%%