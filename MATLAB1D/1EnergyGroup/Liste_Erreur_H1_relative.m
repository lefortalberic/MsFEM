function Erreur_H1_list = Liste_Erreur_H1_relative(num_idee,Nbpt_H_list,epsilon,Nbpt_ref,Vecteur_propre_ref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul de la solution MsFEM par la méthode 'num_idee'
% Renvoie l'erreur H1 entre cette solution et la solution de reference vect_propre_ref


fprintf('\n Calcul de l evolution des erreurs H1 en fonction de H pour la methode %d \n', num_idee);

h_ref = 1/(Nbpt_ref -1);

Erreur_H1_list=[];

for Nbpt_H=Nbpt_H_list

    tic;
    fprintf('\n Calcul en cours pour %d nombre de points \n', Nbpt_H);

    H=1/(Nbpt_H-1);
    Nbpt_h = (Nbpt_ref - 1)/(Nbpt_H-1) +1;

    % Matrice des fonctions de forme Khi_i
    % ----------------------
    Khi_global = zeros(Nbpt_h,Nbpt_H-2,2);
    for l=1:(Nbpt_H-2)
        [fct_forme_gauche , fct_forme_droite] = Khi_i(Nbpt_H,Nbpt_h,epsilon,l,num_idee);

        Khi_global(:,l,1)= fct_forme_gauche ;
        Khi_global(:,l,2)= fct_forme_droite ;
    end

    % calcul des matrices EF
    % ----------------------
    MM_H=matM(Khi_global,Nbpt_H,Nbpt_h);
    MM_sigma_H=matM_sigma(Khi_global,Nbpt_H,Nbpt_h,epsilon);
    KK_H=matK(Khi_global,Nbpt_H,Nbpt_h,epsilon);

    % Matrice à inverser
    % --------------

    %Solution MsFEM
    AA_H = MM_sigma_H + KK_H*(epsilon^2) ;
    AA_H(1,1)=1; AA_H(Nbpt_H,Nbpt_H)=1;

    [EV_H,DV_H] = eigs(AA_H,MM_H,1,'smallestabs');

    [~,arg_min_H] = min(diag(DV_H));
    Vecteur_propre_H = EV_H(:,arg_min_H);

    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre_H(fix(Nbpt_H/2))<0)
        Vecteur_propre_H = -Vecteur_propre_H;
    end

    % % Rescale pour une norme L^2 = 1
    Norme_L2 = sqrt(Vecteur_propre_H'*MM_H*Vecteur_propre_H);
    Vecteur_propre_H = Vecteur_propre_H/Norme_L2;

    %Calculs de toutes les erreurs

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
    Erreur_H1_relative = sqrt(Erreur_H1_relative) ;

    Erreur_H1_list = [Erreur_H1_list , Erreur_H1_relative]; %#ok<AGROW> 

    tempsEcoule = toc;
    disp(['Le temps de calcul de la solution approchee pour ', num2str(Nbpt_H), ' points macroscopiques est de ', num2str(tempsEcoule), ' secondes.']);

end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
