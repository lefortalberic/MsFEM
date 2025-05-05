function Erreur_H1_list = Liste_Erreur_H1_relative_pb_veps(num_idee,Nbpt_H_list,epsilon,Nbpt_ref,Vecteur_propre_ref_1,Vecteur_propre_ref_2)
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
    h = 1/(Nbpt_h -1);
    Khi_global_1 = zeros(Nbpt_h,Nbpt_H-2,2);
    Khi_global_2 = zeros(Nbpt_h,Nbpt_H-2,2);
    for l=1:(Nbpt_H-2)
        [ UU_g_1, UU_g_2 , UU_d_1, UU_d_2] = Khi_i(Nbpt_H,Nbpt_h,epsilon,l,num_idee);
    
        Khi_global_1(:,l,1)= UU_g_1 ;
        Khi_global_1(:,l,2)= UU_d_1 ;
    
        Khi_global_2(:,l,1)= UU_g_2 ;
        Khi_global_2(:,l,2)= UU_d_2 ;
    end
    
    [~,Psi_1_adjoint,Psi_2_adjoint] = Solution_pb_spectral_adjoint(Nbpt_ref);
    [~,Psi_1,Psi_2] = Solution_pb_spectral(Nbpt_ref);
    
    Q_tilde_11 = Q_tilde(Psi_1,Psi_2,Psi_1_adjoint,Psi_2_adjoint,1);
    Q_tilde_22 = Q_tilde(Psi_1,Psi_2,Psi_1_adjoint,Psi_2_adjoint,4);
    Q_tilde_21 = - Q_tilde_22;
    Q_tilde_12 = - Q_tilde_11;
    
    vecteur_J1 = J(Psi_1,Psi_1_adjoint,1);
    vecteur_J2 = J(Psi_2,Psi_2_adjoint,2);
    
    MM_1=matM(Khi_global_1,Psi_1 , Psi_1_adjoint,Nbpt_H,Nbpt_h, epsilon);
    MM_2=matM(Khi_global_2,Psi_2 , Psi_2_adjoint,Nbpt_H,Nbpt_h, epsilon);
    
    KK_1=matK(Khi_global_1,Psi_1 , Psi_1_adjoint,Nbpt_H,Nbpt_h,epsilon, 1);
    KK_2=matK(Khi_global_2,Psi_2 , Psi_2_adjoint,Nbpt_H,Nbpt_h,epsilon, 2);
    
    MM_Q_tilde_11=matM_Q_tilde(Khi_global_1,Khi_global_1, Q_tilde_11, Nbpt_H,Nbpt_h,epsilon);
    MM_Q_tilde_12=matM_Q_tilde(Khi_global_1,Khi_global_2,Q_tilde_12, Nbpt_H,Nbpt_h,epsilon);
    MM_Q_tilde_21=matM_Q_tilde(Khi_global_2,Khi_global_1, Q_tilde_21, Nbpt_H,Nbpt_h,epsilon);
    MM_Q_tilde_22=matM_Q_tilde(Khi_global_2,Khi_global_2,Q_tilde_22, Nbpt_H,Nbpt_h,epsilon);
    
    CC_1 = matC(Khi_global_1,vecteur_J1,Nbpt_H,Nbpt_h, epsilon);
    CC_2 = matC(Khi_global_2,vecteur_J2,Nbpt_H,Nbpt_h, epsilon);
    
    %Assemblage de toutes les matrices
    
    Zeros = sparse(Nbpt_H,Nbpt_H); % matrice de masse
    
    KK = [KK_1 , Zeros ; Zeros ,KK_2 ];
    MM_Q_tilde = [MM_Q_tilde_11 , MM_Q_tilde_12 ; MM_Q_tilde_21 , MM_Q_tilde_22];
    MM = [MM_1 , Zeros ; Zeros ,MM_2 ];
    CC = [CC_1 , Zeros ; Zeros ,CC_2];
    
    AA = KK + MM_Q_tilde/(epsilon^2) + CC/epsilon;
    AA(1,1)=1; AA(2*Nbpt_H,2*Nbpt_H)=1;
    AA(Nbpt_H,Nbpt_H)=1 ; AA(Nbpt_H+1,Nbpt_H+1)=1;
    
    %[EV,DV] = eigs(AA,MM,5,'smallestabs');
    [EV,DV] = eigs(AA,MM,1,'smallestabs');
    
    
    [~,arg_min] = min(diag(DV));
    Vecteur_propre_H = EV(:,arg_min);
    
    Vecteur_propre_H_1 = Vecteur_propre_H(1:Nbpt_H);
    Vecteur_propre_H_2 = Vecteur_propre_H((Nbpt_H+1):2*Nbpt_H);
    
    %On a un VP de norme L^2 unitaire. On choisis le VP positif.
    if (Vecteur_propre_H_1(fix(Nbpt_H/2))<0)
        Vecteur_propre_H_1 = -Vecteur_propre_H_1;
    end
    if (Vecteur_propre_H_2(fix(Nbpt_H/2))<0)
        Vecteur_propre_H_2 = -Vecteur_propre_H_2;
    end
    % % Rescale pour une norme L^2 = 1
    MM_1_constant=matM(Khi_global_1, ones(Nbpt_ref,1) , ones(Nbpt_ref,1),Nbpt_H,Nbpt_h, epsilon);
    MM_2_constant=matM(Khi_global_2, ones(Nbpt_ref,1) , ones(Nbpt_ref,1),Nbpt_H,Nbpt_h, epsilon);
    
    Norme_L2_1 = sqrt(Vecteur_propre_H_1'*MM_1_constant*Vecteur_propre_H_1);
    Vecteur_propre_H_1 = Vecteur_propre_H_1/Norme_L2_1;
    
    Norme_L2_2 = sqrt(Vecteur_propre_H_2'*MM_2_constant*Vecteur_propre_H_2);
    Vecteur_propre_H_2 = Vecteur_propre_H_2/Norme_L2_2;
    
    Erreur_H1_1 = 0;
    KK_ref_A_constant = matK_unitaire(Nbpt_h,h_ref);
    Norme_H1_vect_propre_ref_1 = 0;

    for i=2:(Nbpt_H-2)
        Vecteur_propre_ref_local = Vecteur_propre_ref_1((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h);
        solution_intervalle_i = zeros(Nbpt_h,1);
        for j=1:Nbpt_h
            solution_intervalle_i(j) = Vecteur_propre_H_1(i)*Khi_global_1(j,i-1,2) + Vecteur_propre_H_1(i+1)*Khi_global_1(j,i,1);
        end
        Norme_H1_vect_propre_ref_1 = Norme_H1_vect_propre_ref_1 + (Vecteur_propre_ref_local'*KK_ref_A_constant*Vecteur_propre_ref_local);
        Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant*(solution_intervalle_i-Vecteur_propre_ref_local));
        Erreur_H1_1 = Erreur_H1_1 + Erreur_H1_locale;
    end

    %Premier intervalle
    i = 1; % On trace la solution entre 0 et H
    Vecteur_propre_ref_local = Vecteur_propre_ref_1((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h);
    solution_intervalle_i = zeros(Nbpt_h,1);
    for j=1:Nbpt_h
        solution_intervalle_i(j) = Vecteur_propre_H_1(i+1)*Khi_global_1(j,i,1);
    end
    Norme_H1_vect_propre_ref_1 = Norme_H1_vect_propre_ref_1 + (Vecteur_propre_ref_local'*KK_ref_A_constant*Vecteur_propre_ref_local);
    Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant*(solution_intervalle_i-Vecteur_propre_ref_local));
    Erreur_H1_1 = Erreur_H1_1 + Erreur_H1_locale;

    %Dernier intervalle
    i = Nbpt_H-1; % On trace la solution entre 1-H et 1
    Vecteur_propre_ref_local = Vecteur_propre_ref_1((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h);
    solution_intervalle_i = zeros(Nbpt_h,1);
    for j=1:Nbpt_h
        solution_intervalle_i(j) = Vecteur_propre_H_1(i)*Khi_global_1(j,i-1,2);
    end
    Norme_H1_vect_propre_ref_1 = Norme_H1_vect_propre_ref_1 + (Vecteur_propre_ref_local'*KK_ref_A_constant*Vecteur_propre_ref_local);
    Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant*(solution_intervalle_i-Vecteur_propre_ref_local));
    Erreur_H1_1 = Erreur_H1_1 + Erreur_H1_locale;

    Erreur_H1_relative_1 = Erreur_H1_1/Norme_H1_vect_propre_ref_1;

    Erreur_H1_2 = 0;
    KK_ref_A_constant = matK_unitaire(Nbpt_h,h_ref);
    Norme_H1_vect_propre_ref_2 = 0;

    for i=2:(Nbpt_H-2)
        Vecteur_propre_ref_local = Vecteur_propre_ref_2((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h);
        solution_intervalle_i = zeros(Nbpt_h,1);
        for j=1:Nbpt_h
            solution_intervalle_i(j) = Vecteur_propre_H_2(i)*Khi_global_2(j,i-1,2) + Vecteur_propre_H_2(i+1)*Khi_global_2(j,i,1);
        end
        Norme_H1_vect_propre_ref_2 = Norme_H1_vect_propre_ref_2 + (Vecteur_propre_ref_local'*KK_ref_A_constant*Vecteur_propre_ref_local);
        Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant*(solution_intervalle_i-Vecteur_propre_ref_local));
        Erreur_H1_2 = Erreur_H1_2 + Erreur_H1_locale;
    end

    %Premier intervalle
    i = 1; % On trace la solution entre 0 et H
    Vecteur_propre_ref_local = Vecteur_propre_ref_2((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h);
    solution_intervalle_i = zeros(Nbpt_h,1);
    for j=1:Nbpt_h
        solution_intervalle_i(j) = Vecteur_propre_H_2(i+1)*Khi_global_2(j,i,1);
    end
    Norme_H1_vect_propre_ref_2 = Norme_H1_vect_propre_ref_2 + (Vecteur_propre_ref_local'*KK_ref_A_constant*Vecteur_propre_ref_local);
    Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant*(solution_intervalle_i-Vecteur_propre_ref_local));
    Erreur_H1_2 = Erreur_H1_2 + Erreur_H1_locale;

    %Dernier intervalle
    i = Nbpt_H-1; % On trace la solution entre 1-H et 1
    Vecteur_propre_ref_local = Vecteur_propre_ref_2((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h);
    solution_intervalle_i = zeros(Nbpt_h,1);
    for j=1:Nbpt_h
        solution_intervalle_i(j) = Vecteur_propre_H_2(i)*Khi_global_2(j,i-1,2);
    end
    Norme_H1_vect_propre_ref_2 = Norme_H1_vect_propre_ref_2 + (Vecteur_propre_ref_local'*KK_ref_A_constant*Vecteur_propre_ref_local);
    Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant*(solution_intervalle_i-Vecteur_propre_ref_local));
    Erreur_H1_2 = Erreur_H1_2 + Erreur_H1_locale;

    Erreur_H1_relative_2 = Erreur_H1_2/Norme_H1_vect_propre_ref_2;

    Erreur_H1_relative = (1/sqrt(2))*sqrt(Erreur_H1_relative_1 + Erreur_H1_relative_2); 
    
    Erreur_H1_list = [Erreur_H1_list , Erreur_H1_relative]; %#ok<AGROW> 

    tempsEcoule = toc;
    disp(['Le temps de calcul de la solution approchee pour ', num2str(Nbpt_H), ' points macroscopiques est de ', num2str(tempsEcoule), ' secondes.']);

end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
