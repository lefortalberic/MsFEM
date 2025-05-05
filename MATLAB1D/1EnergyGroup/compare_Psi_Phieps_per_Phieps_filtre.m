% ------------% ------------% ------------% ------------% ------------
% Code pour comparer Psi_eps et Phi_eps avec methode du filtre.
% Il permet de comparer/tracer différentes estimées d'erreurs sur les
% vecteurs propres et sur les valeurs propres
% 
% ------------% ------------% ------------% ------------% ------------

Nbpt_h = 6001;
N_Liste = (3.4:0.5:30);
epsilon_liste =1./N_Liste;
Erreur_sur_omega =0;

h = 1/(Nbpt_h -1);
Nbpt_ref = (Nbpt_h-1)+1;
h_ref = 1/(Nbpt_ref-1);
%Maillage
XX = zeros(Nbpt_h,1);     %vecteur des coordonnées des noeuds
for j=1:Nbpt_ref
    XX(j)=(j-1)*h;
end
MM_ref = matM_ref(Nbpt_h);

%Maillage
Vect_Un = zeros(Nbpt_ref,1);
PhiXX = zeros(Nbpt_ref,1);     %vecteur des coordonnées des noeuds
for j=1:Nbpt_ref
    PhiXX(j)=Phi((j-1)*h);
    Vect_Un(j) = 1;
end
Norme_Phi = sqrt(PhiXX'*MM_ref*PhiXX);
%tracé du coefficient de diffusion
AXX = zeros(Nbpt_h,1);     %vecteur des coordonnées des noeuds
Sigma_epsXX = zeros(Nbpt_h,1);     %vecteur des coordonnées des noeuds
Sigma_yXX = zeros(Nbpt_h,1);     %vecteur des coordonnées des noeuds
for j=1:Nbpt_h
    epsilontmp = epsilon_liste(1);
    AXX(j)=A((j-1)*h,epsilontmp);
    Sigma_epsXX(j)=Sigma((j-1)*h,epsilontmp);
    Sigma_yXX(j)=Sigma((j-1)*h,1);
end
Sigma_epsXXmean = mean(Sigma_epsXX);
Sigma_yXXmean = mean(Sigma_yXX);


KK_ref_A_constant = matK_unitaire((1*(Nbpt_h-1)+1),h_ref);

KK_ref_A_constant_filtre = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
for l=1:(Nbpt_h-1)   %La première et dernière ligne doivent rester nulle
    KK_ref_A_constant_filtre(l,l)=(1/h)*(Phi((l-1-0.5)*h) + Phi((l-1+0.5)*h)); %approximation de l'intégrale de a(.)
    KK_ref_A_constant_filtre(l,l+1)=(-1/h)*Phi((l-1+0.5)*h); %#ok<*SPRIX>
    KK_ref_A_constant_filtre(l+1,l)=(-1/h)*Phi((l-1+0.5)*h);
end % for l

KK_ref_A_constant_filtre(1,1)=(1/h)*Phi((0.5)*h);
KK_ref_A_constant_filtre(Nbpt_h,Nbpt_h)=(1/h)*Phi(1-(0.5)*h);

KK_interieur_filtre = sparse((Nbpt_h-1)/3+1,(Nbpt_h-1)/3+1); % matrice de rigidite
for l=(Nbpt_h-1)/3:2*(Nbpt_h-1)/3-1  %La première et dernière ligne doivent rester nulle
    j=l-(Nbpt_h-1)/3 + 1;
    KK_interieur_filtre(j,j)=(1/h)*(Phi((l-1-0.5)*h) + Phi((l-1+0.5)*h)); %approximation de l'intégrale de a(.)
    KK_interieur_filtre(j,j+1)=(-1/h)*Phi((l-1+0.5)*h); %#ok<*SPRIX>
    KK_interieur_filtre(j+1,j)=(-1/h)*Phi((l-1+0.5)*h);
end % for l

KK_interieur_filtre(1,1)=(1/h)*Phi(((Nbpt_h-1)/3 -1+ 0.5)*h);
KK_interieur_filtre((Nbpt_h-1)/3 +1,(Nbpt_h-1)/3 +1)=(1/h)*Phi((2*(Nbpt_h-1)/3 -1 - 0.5)*h);

KK_interieur = KK_ref_A_constant((Nbpt_h-1)/3:2*(Nbpt_h-1)/3,(Nbpt_h-1)/3:2*(Nbpt_h-1)/3);
KK_interieur(1,1)=KK_interieur(1,1)/2;
KK_interieur((Nbpt_h-1)/3 +1,(Nbpt_h-1)/3 +1) = KK_interieur((Nbpt_h-1)/3 +1,(Nbpt_h-1)/3 +1)/2;

Liste_Erreur_H1_periodique=[];
Liste_Erreur_H1_filtre=[];
Liste_multiplicateur_Lagrange_filtre = [];
Liste_Erreur_H1_filtre_norme_filtre=[];
Liste_Erreur_VP_periodique=[];
Liste_Erreur_max_periodique=[];
Liste_Erreur_max_filtre=[];
Liste_Erreur_max_derive_filtre=[];

Liste_Erreur_max_filtre_norme_filtre=[];
Liste_Erreur_max_filtre_norme_racine_filtre=[];

Liste_Erreur_VP_filtre=[];
Liste_Erreur_H1_periodique_2cellules=[];
Liste_Erreur_H1_filtre_2cellules=[];

Liste_w_prime_mid = [];
Liste_derive_vect_propre_eps_filtre_interieur_mid = [];
for epsilon=epsilon_liste


        % Phieps_periodique
        % ------------% ------------% ------------% ------------% ------------
        
        Nbpt_per = 1*(Nbpt_h-1)+1-1;
        %Matrice EF
        MM = matM_per(Nbpt_per);
        MM_sigma = sparse(Nbpt_per,Nbpt_per); % matrice de masse avec sigma
        KK = sparse(Nbpt_per,Nbpt_per); % matrice de rigidite
        
        for l=2:(Nbpt_per-1)   %La première et dernière ligne doivent rester nulle
        
            KK(l,l)=(1/h)*(A((l-1-0.5)*h,epsilon)+A((l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
            KK(l,l+1)=(-1/h)*A((l-1+0.5)*h,epsilon); %#ok<*SPRIX>
            KK(l+1,l)=(-1/h)*A((l-1+0.5)*h,epsilon);
        
            MM_sigma(l,l)=(h/3)*(Sigma((l-1-0.5)*h,epsilon) ...
                +Sigma((l-1+0.5)*h,epsilon));

            MM_sigma(l,l+1)=h*(1/6)*Sigma((l-1+0.5)*h,epsilon);
            MM_sigma(l+1,l)=h*(1/6)*Sigma((l-1+0.5)*h,epsilon);
        
        end % for l
        
        KK(1,1)=(1/h)*(A(0.5*h,epsilon)+A(1-0.5*h,epsilon));
        KK(1,2)=(-1/h)*A(0.5*h,epsilon);
        KK(2,1)=(-1/h)*A(0.5*h,epsilon);
        KK(1,Nbpt_per)=(-1/h)*A(1-0.5*h,epsilon);
        KK(Nbpt_per,1)=(-1/h)*A(1-0.5*h,epsilon);
        KK(Nbpt_per,Nbpt_per)=(1/h)*(A((Nbpt_per-1-0.5)*h,epsilon)+A((Nbpt_per-1+0.5)*h,epsilon));
        
        MM_sigma(1,1)=(h/3)*(Sigma(0.5*h,epsilon) ...
            +Sigma(1-0.5*h,epsilon));
        MM_sigma(1,2)=h*(1/6)*Sigma(0.5*h,epsilon);
        MM_sigma(2,1)=h*(1/6)*Sigma(0.5*h,epsilon);
        MM_sigma(1,Nbpt_per)=h*(1/6)*Sigma(1-0.5*h,epsilon);
        MM_sigma(Nbpt_per,1)=h*(1/6)*Sigma(1-0.5*h,epsilon);
        MM_sigma(Nbpt_per,Nbpt_per)=(h/3)*(Sigma((Nbpt_per-1-0.5)*h,epsilon) ...
            +Sigma((Nbpt_per-1+0.5)*h,epsilon));
        
        AA = MM_sigma + KK*(epsilon*epsilon) ;
        
        [EV,DV_periodique] = eigs(AA,MM,1,'smallestabs');
        
        [~,arg_min] = min(diag(DV_periodique));
        Vecteur_propre_H_periodic = EV(:,arg_min);
        
        %On a un VP de norme L^2 unitaire. On choisis le VP positif.
        if (Vecteur_propre_H_periodic(fix(Nbpt_h/2))<0)
            Vecteur_propre_H_periodic = -Vecteur_propre_H_periodic;
        end
        
        Vecteur_propre_H_periodic = [Vecteur_propre_H_periodic;Vecteur_propre_H_periodic(1)]; %#ok<AGROW> 
        
        % % Rescale pour une norme L^2 = 1
        Norme_L2 = sqrt(Vecteur_propre_H_periodic'*MM_ref*Vecteur_propre_H_periodic);
        vect_propre_eps_periodique = Vecteur_propre_H_periodic/Norme_L2;
        
        % Phieps_avec_filtrage
        % ------------% ------------% ------------% ------------% ------------
        MM = sparse(Nbpt_h,Nbpt_h); % matrice de masse
        KK = sparse(Nbpt_h,Nbpt_h); % matrice de rigidité
        MM_sigma = sparse(Nbpt_h,Nbpt_h); % matrice de masse
        % boucle pour les matrices EF
        % ------------------------
        for l=1:(Nbpt_h-1)   %La première et dernière ligne doivent rester nulle
        
            MM(l,l)=h*(1/3)*(Phi((l-1-0.5)*h) + Phi((l-1+0.5)*h));
            MM(l,l+1)=h*(1/6)*Phi((l-1+0.5)*h);
            MM(l+1,l)=h*(1/6)*Phi((l-1+0.5)*h);

            KK(l,l)=(1/h)*(A((l-1-0.5)*h,epsilon)*Phi((l-1-0.5)*h)+A((l-1+0.5)*h,epsilon)*Phi((l-1+0.5)*h)); %approximation de l'intégrale de a(.)
            KK(l,l+1)=(-1/h)*A((l-1+0.5)*h,epsilon)*Phi((l-1+0.5)*h);
            KK(l+1,l)=(-1/h)*A((l-1+0.5)*h,epsilon)*Phi((l-1+0.5)*h);

            MM_sigma(l,l)=(h/3)*(Sigma((l-1-0.5)*h,epsilon)*Phi((l-1-0.5)*h)+Sigma((l-1+0.5)*h,epsilon)*Phi((l-1+0.5)*h));
            MM_sigma(l,l+1)=h*(1/6)*Sigma((l-1+0.5)*h,epsilon)*Phi((l-1+0.5)*h);
            MM_sigma(l+1,l)=h*(1/6)*Sigma((l-1+0.5)*h,epsilon)*Phi((l-1+0.5)*h);
            
        end % for l
        MM(1,1)=h*(1/3)*Phi((0.5)*h);
        MM(Nbpt_h,Nbpt_h)=h*(1/3)*Phi(1-(0.5)*h);
        
        KK(1,1)=(1/h)*(A(0.5*h,epsilon))*Phi((0.5)*h);
        KK(Nbpt_h,Nbpt_h)=(1/h)*(A(1-0.5*h,epsilon))*Phi(1-(0.5)*h);

        MM_sigma(1,1)= (h/3)*(Sigma((0.5)*h,epsilon))*Phi((0.5)*h);
        MM_sigma(Nbpt_h,Nbpt_h)=(h/3)*(Sigma((Nbpt_h-1-0.5)*h,epsilon))*Phi(1-(0.5)*h);

        BB = zeros(Nbpt_h,1); %Matrice contrainte périodicité
        
        for l=2:(Nbpt_h-1)
            BB(l)= Phi((l-1-0.5)*h)-Phi((l-1+0.5)*h);
        end
        BB(1)=-1*Phi((0.5)*h);
        BB(Nbpt_h)=1*Phi(1-(0.5)*h);

        AA = MM_sigma + KK*(epsilon*epsilon) ;
        
        AA_f = [AA , BB ; BB' , 0];
        
        MM_f = [MM , 0*BB ; 0*BB' , 0] ;

        [EV,DV_filtre] = eigs(AA_f,MM_f,1,'smallestabs');
        
        %[EV,DV_filtre] = eigs(AA,MM,1,'smallestabs');
        
        [~,arg_min] = min(diag(DV_filtre));
        Vecteur_propre_H_periodic_filtre = EV(:,arg_min);
        %TMP = AA_f*Vecteur_propre_H_periodic_filtre - valeur_propre_eps_min*MM_f*Vecteur_propre_H_periodic_filtre;
        
        Multiplicateur_Lagrange_filtre = Vecteur_propre_H_periodic_filtre(Nbpt_h+1);
        Vecteur_propre_H_periodic_filtre = Vecteur_propre_H_periodic_filtre(1:Nbpt_h);
        %TMP = Vecteur_propre_H_periodic_filtre'*MM*Vecteur_propre_H_periodic_filtre
        %On a un VP de norme L^2 unitaire. On choisis le VP positif.
        if (Vecteur_propre_H_periodic_filtre(fix(Nbpt_h/2))<0)
            Vecteur_propre_H_periodic_filtre = -Vecteur_propre_H_periodic_filtre;
        end

%         % % Rescale pour une norme L^2 = 1
%         Norme_L2 = sqrt(Vecteur_propre_H_periodic_filtre'*MM_ref*Vecteur_propre_H_periodic_filtre);
%         vect_propre_eps_filtre = Vecteur_propre_H_periodic_filtre/Norme_L2;

        test = zeros(Nbpt_h,1);
        for i=1:(Nbpt_h-1)
            grad_WRfiltre = (Vecteur_propre_H_periodic_filtre(i+1) - Vecteur_propre_H_periodic_filtre(i))/h;
            test(i) = Phi((i-1-0.5)*h)*grad_WRfiltre*A((i-1+0.5)*h,epsilon);
        end

        vect_propre_eps_filtre = Vecteur_propre_H_periodic_filtre;
        vect_propre_eps_filtre_tmp =zeros(Nbpt_h,1);
        for l=1:(Nbpt_h)
            vect_propre_eps_filtre_tmp(l)= vect_propre_eps_filtre(l)*sqrt(Phi((l-1)*h));
        end
        % Psi_eps
        % ------------% ------------% ------------% ------------% ------------
        
        MM = matM_per(Nbpt_h);
        MM_sigma = matM_sigma_per(Nbpt_h); % matrice de masse avec sigma
        KK = matK_per(Nbpt_h); % matrice de rigidite
        
        AA = MM_sigma + KK ;
        
        [EV,DV_Psi] = eigs(AA,MM,1,'smallestabs');
        
        [Valeur_propre_min,arg_min] = min(diag(DV_Psi));
        Vecteur_propre = EV(:,arg_min);
        %On a un VP de norme L^2 unitaire. On choisis le VP positif.
        if (Vecteur_propre(1)<0)
            Vecteur_propre = -Vecteur_propre;
        end
        VV_eps = zeros(Nbpt_h,1);
        for i=1:Nbpt_per+1
            xi_g = (i-1)*h ;
            yi_g = rem(xi_g/epsilon,1) ; % y vaut x/epsilon modulo 1
            yk_g = fix(yi_g*(Nbpt_h-1)) +1; %partie entière de y*(Nbpt-1) (indice pour les matrices)
    
            t_g = yi_g*(Nbpt_h-1) - fix(yi_g*(Nbpt_h-1)) ;
            VV_eps(i) = Vecteur_propre(yk_g)*(1-t_g) +Vecteur_propre(yk_g+1)*t_g;

        end
        
        % % Rescale pour une norme L^2 = 1
        Norme_L2 = sqrt(VV_eps'*MM_ref*VV_eps);
        VV_eps = VV_eps/Norme_L2;
        % % Rescale pour une norme L^2 = 1
        Norme_L2 = sqrt(vect_propre_eps_filtre'*MM_ref*vect_propre_eps_filtre);
        vect_propre_eps_filtre = vect_propre_eps_filtre/Norme_L2;
% 
%         Developpement_approchee_psi_prime=zeros(Nbpt_h,1);
%         for j=2:Nbpt_h
%             Developpement_approchee_psi_prime(j) = Developpement_approchee_psi_prime(j-1)+ h*(Sigma_epsXX(j)-Sigma_yXXmean);
%         end
%         derivepsi = (VV_eps(2:length(VV_eps))-VV_eps(1:(length(VV_eps)-1)))/h;
%         
        if (Erreur_sur_omega)
            Erreur_H1_Psi_Phi_periodique = sqrt(((VV_eps-vect_propre_eps_periodique)'*KK_ref_A_constant*(VV_eps-vect_propre_eps_periodique))/(VV_eps'*KK_ref_A_constant*VV_eps)); %#ok<UNRCH> 
    
            Erreur_H1_Psi_Phi_filtre = sqrt(((VV_eps-vect_propre_eps_filtre)'*KK_ref_A_constant*(VV_eps-vect_propre_eps_filtre))/(VV_eps'*KK_ref_A_constant*VV_eps));
    
            Erreur_H1_Psi_Phi_filtre_norme_filtre = sqrt(((VV_eps-vect_propre_eps_filtre)'*KK_ref_A_constant_filtre*(VV_eps-vect_propre_eps_filtre))/(VV_eps'*KK_ref_A_constant_filtre*VV_eps));

            Erreur_max_Psi_Phi_periodique = max( abs(VV_eps-vect_propre_eps_periodique) ) / max( abs(VV_eps) );

            Erreur_max_Psi_Phi_filtre = max( abs(VV_eps-vect_propre_eps_filtre) ) / max( abs(VV_eps) );

            diff_filtre =zeros(Nbpt_h,1);
            diff_racine_filtre =zeros(Nbpt_h,1);
            psi_filtre = zeros(Nbpt_h,1);
            psi_filtre_racine = zeros(Nbpt_h,1);
            for l=1:(Nbpt_h)
                diff_filtre(l)= (VV_eps(l)-vect_propre_eps_filtre(l))*Phi((l-1)*h);
                diff_racine_filtre(l) = (VV_eps(l)-vect_propre_eps_filtre(l))*sqrt(Phi((l-1)*h));
                psi_filtre(l) = (VV_eps(l))*Phi((l-1)*h);
                psi_filtre_racine(l) = (VV_eps(l))*sqrt(Phi((l-1)*h));
            end

            %Erreur_max_Psi_Phi_filtre_norme_filtre = max( abs(diff_filtre) ) / max( abs(psi_filtre) );
            %Erreur_max_Psi_Phi_filtre_norme_filtre_racine = max( abs(diff_racine_filtre) ) / max( abs(psi_filtre_racine) );

        else
            MM_ref_interieur = matM_ref(2*(Nbpt_h-1)/3-(Nbpt_h-1)/3 +1); %#ok<UNRCH> 
            VV_eps_interieur = VV_eps((Nbpt_h-1)/3:2*(Nbpt_h-1)/3);
            vect_propre_eps_periodique_interieur = vect_propre_eps_periodique((Nbpt_h-1)/3:2*(Nbpt_h-1)/3);
            vect_propre_eps_filtre_interieur = vect_propre_eps_filtre((Nbpt_h-1)/3:2*(Nbpt_h-1)/3);

            % % Rescale pour une norme L^2 = 1
            Norme_L2_VV_eps = sqrt(VV_eps_interieur'*MM_ref_interieur*VV_eps_interieur);
            VV_eps_interieur = VV_eps_interieur/Norme_L2_VV_eps;

            Norme_L2_vect_propre_eps_periodique = sqrt(vect_propre_eps_periodique_interieur'*MM_ref_interieur*vect_propre_eps_periodique_interieur);
            vect_propre_eps_periodique_interieur = vect_propre_eps_periodique_interieur/Norme_L2_vect_propre_eps_periodique;

            Norme_L2_vect_propre_eps_filtre = sqrt(vect_propre_eps_filtre_interieur'*MM_ref_interieur*vect_propre_eps_filtre_interieur);
            vect_propre_eps_filtre_interieur = vect_propre_eps_filtre_interieur/Norme_L2_vect_propre_eps_filtre;

            Erreur_H1_Psi_Phi_periodique = sqrt(((VV_eps_interieur-vect_propre_eps_periodique_interieur)'*KK_interieur*(VV_eps_interieur-vect_propre_eps_periodique_interieur))/(VV_eps_interieur'*KK_interieur*VV_eps_interieur)); %#ok<UNRCH> 
    
            Erreur_H1_Psi_Phi_filtre = sqrt(((VV_eps_interieur-vect_propre_eps_filtre_interieur)'*KK_interieur*(VV_eps_interieur-vect_propre_eps_filtre_interieur))/(VV_eps_interieur'*KK_interieur*VV_eps_interieur));
    
            Erreur_H1_Psi_Phi_filtre_norme_filtre = sqrt(((VV_eps_interieur-vect_propre_eps_filtre_interieur)'*KK_interieur_filtre*(VV_eps_interieur-vect_propre_eps_filtre_interieur))/(VV_eps_interieur'*KK_interieur_filtre*VV_eps_interieur));
        
            Erreur_max_Psi_Phi_periodique = max( abs(VV_eps_interieur-vect_propre_eps_periodique_interieur) ) / max( abs(VV_eps_interieur) );

            Erreur_max_Psi_Phi_filtre = max( abs(VV_eps_interieur-vect_propre_eps_filtre_interieur) ) / max( abs(VV_eps_interieur) );
        
            a = (VV_eps_interieur-vect_propre_eps_filtre_interieur);
            derive = (a(2:length(a))-a(1:(length(a)-1)))/h;
            derivepsi = (VV_eps_interieur(2:length(VV_eps_interieur))-VV_eps_interieur(1:(length(VV_eps_interieur)-1)))/h;
            derivephi_filtre = (vect_propre_eps_filtre_interieur(2:length(vect_propre_eps_filtre_interieur))-vect_propre_eps_filtre_interieur(1:(length(vect_propre_eps_filtre_interieur)-1)))/h;
            Erreur_max_derive_Psi_Phi_filtre = max( abs(derive) ) / max( abs(derivepsi) );

            yl_min = (Nbpt_ref-1)/2 - (fix(epsilon*((Nbpt_ref-1)))+1);
            yl_max = (Nbpt_ref-1)/2 + (fix(epsilon*((Nbpt_ref-1))));
            
            VV_eps_2_cellules = VV_eps(yl_min:yl_max);
            vect_propre_eps_periodique_2_celulles = vect_propre_eps_periodique(yl_min:yl_max);
            vect_propre_eps_filtre_2_cellules = vect_propre_eps_filtre(yl_min:yl_max);
            
            MM_ref_2_cellules = matM_ref(yl_max - yl_min + 1);

            % % Rescale pour une norme L^2 = 1
            Norme_L2_VV_eps = sqrt(VV_eps_2_cellules'*MM_ref_2_cellules*VV_eps_2_cellules);
            VV_eps_2_cellules = VV_eps_2_cellules/Norme_L2_VV_eps;
            Norme_L2_vect_propre_eps_periodique = sqrt(vect_propre_eps_periodique_2_celulles'*MM_ref_2_cellules*vect_propre_eps_periodique_2_celulles);
            vect_propre_eps_periodique_2_celulles = vect_propre_eps_periodique_2_celulles/Norme_L2_vect_propre_eps_periodique;
            Norme_L2_vect_propre_eps_filtre = sqrt(vect_propre_eps_filtre_2_cellules'*MM_ref_2_cellules*vect_propre_eps_filtre_2_cellules);
            vect_propre_eps_filtre_2_cellules = vect_propre_eps_filtre_2_cellules/Norme_L2_vect_propre_eps_filtre;


            KK_interieur_tmp = KK_ref_A_constant(yl_min:yl_max,yl_min:yl_max);
            KK_interieur_tmp(1,1)=KK_interieur_tmp(1,1)/2;
            KK_interieur_tmp(yl_max - yl_min + 1,yl_max - yl_min + 1) = KK_interieur_tmp(yl_max - yl_min + 1,yl_max - yl_min + 1)/2;

            Erreur_H1_Psi_Phi_periodique_2_cellules = sqrt(((VV_eps_2_cellules-vect_propre_eps_periodique_2_celulles)'*KK_interieur_tmp*(VV_eps_2_cellules-vect_propre_eps_periodique_2_celulles))/(VV_eps_2_cellules'*KK_interieur_tmp*VV_eps_2_cellules)); %#ok<UNRCH>
            
            Erreur_H1_Psi_Phi_filtre_2_cellules = sqrt(((VV_eps_2_cellules-vect_propre_eps_filtre_2_cellules)'*KK_interieur_tmp*(VV_eps_2_cellules-vect_propre_eps_filtre_2_cellules))/(VV_eps_2_cellules'*KK_interieur_tmp*VV_eps_2_cellules)); %#ok<UNRCH>


        end
        Erreur_VP_periodique = abs(DV_periodique-DV_Psi)/abs(DV_Psi);
        Erreur_VP_filtre = abs(DV_filtre-DV_Psi)/abs(DV_Psi);

        Liste_Erreur_H1_periodique=[Liste_Erreur_H1_periodique,Erreur_H1_Psi_Phi_periodique]; %#ok<AGROW>
        Liste_Erreur_H1_filtre=[Liste_Erreur_H1_filtre,Erreur_H1_Psi_Phi_filtre]; %#ok<AGROW> 
        Liste_Erreur_H1_filtre_norme_filtre=[Liste_Erreur_H1_filtre_norme_filtre,Erreur_H1_Psi_Phi_filtre_norme_filtre]; %#ok<AGROW> 

        Liste_Erreur_max_periodique=[Liste_Erreur_max_periodique,Erreur_max_Psi_Phi_periodique]; %#ok<AGROW>
        Liste_Erreur_max_filtre=[Liste_Erreur_max_filtre,Erreur_max_Psi_Phi_filtre]; %#ok<AGROW>
        Liste_Erreur_max_derive_filtre = [Liste_Erreur_max_derive_filtre,Erreur_max_derive_Psi_Phi_filtre];

        %Liste_Erreur_max_filtre_norme_filtre=[Liste_Erreur_max_filtre_norme_filtre,Erreur_max_Psi_Phi_filtre_norme_filtre];
        %Liste_Erreur_max_filtre_norme_racine_filtre=[Liste_Erreur_max_filtre_norme_racine_filtre, Erreur_max_Psi_Phi_filtre_norme_filtre_racine ];

        Liste_w_prime_mid = [Liste_w_prime_mid, abs(derive(floor(length(derive)/2)))];
        Liste_derive_vect_propre_eps_filtre_interieur_mid = [Liste_derive_vect_propre_eps_filtre_interieur_mid, abs(derivephi_filtre(floor(length(derivephi_filtre)/2)))];
        
        Liste_Erreur_VP_periodique=[Liste_Erreur_VP_periodique,Erreur_VP_periodique]; %#ok<AGROW>
        Liste_Erreur_VP_filtre=[Liste_Erreur_VP_filtre,Erreur_VP_filtre]; %#ok<AGROW>

        Liste_Erreur_H1_periodique_2cellules=[Liste_Erreur_H1_periodique_2cellules,Erreur_H1_Psi_Phi_periodique_2_cellules]; %#ok<AGROW>
        Liste_Erreur_H1_filtre_2cellules=[Liste_Erreur_H1_filtre_2cellules,Erreur_H1_Psi_Phi_filtre_2_cellules]; %#ok<AGROW>

        Liste_multiplicateur_Lagrange_filtre = [Liste_multiplicateur_Lagrange_filtre,Multiplicateur_Lagrange_filtre];
end

Sigma_moyXX = zeros(Nbpt_h,1);     %vecteur des coordonnées des noeuds
for j=1:Nbpt_h
    epsilontmp = epsilon_liste(1);
    Sigma_moyXX(j)=Sigma((j-1)*h,epsilontmp);
end
Sigma_moyXXmean = mean(Sigma_moyXX);

% a = (Vecteur_propre);
% derive = (a(2:length(a))-a(1:length(a)-1))/h;

plot(XX,Sigma_moyXX,'LineWidth',4)
hold on;
xlabel("X",'FontSize',25)
title(sprintf("tilde Sigma_{eps}"),'FontSize',25)
set(gca,'FontSize',30)
hold off;

%%%Trace de l'erreur H1%%%%%%%%%%%%%%%%%%%%%
plot(log10(1./(epsilon_liste)),log10(Liste_Erreur_H1_periodique),'LineWidth',4,'marker','x','markersize',18)
hold on;
plot(log10(1./(epsilon_liste)),log10(Liste_Erreur_H1_filtre),'LineWidth',4,'marker','x','markersize',18)
hold on;
% plot(log10(1./(epsilon_liste)),log10(Liste_Erreur_H1_filtre_norme_filtre),'LineWidth',4,'marker','x','markersize',18)
% hold on;
% plot(log10(1./(epsilon_liste)),+0*log10(1./(epsilon_liste))-0.7,'LineWidth',4,'marker','x','markersize',10)
% hold on;
plot(log10(1./(epsilon_liste)),-1*log10(1./(epsilon_liste))-0.3,'LineWidth',4,'marker','x','markersize',10)
hold on;
% Tracé des lignes verticales
ylims = get(gca,'YLim'); %#ok<NASGU> 
hold on;
%legend('CL periodique','avec filtre','avec filtre + norme filtre','pente +0','pente -1','FontSize',25)
legend('CL periodique','avec filtre','pente -1','FontSize',25)
xlabel("-log10(eps)",'FontSize',25)
ylabel("log10(Err_{relative} H1)",'FontSize',25)
title(sprintf("Evolution de l'erreur relative H1(1/3,2/3) en fonction de eps"),'FontSize',25)
set(gca,'FontSize',30)
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Norme L infini%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(log10(1./(epsilon_liste)),log10((Liste_Erreur_max_periodique)),'LineWidth',4,'marker','x','markersize',18)
hold on;
plot(log10(1./(epsilon_liste)),log10((Liste_Erreur_max_derive_filtre)),'LineWidth',4,'marker','x','markersize',18)
hold on;
plot(log10(1./(epsilon_liste)),-0.5*log10(1./(epsilon_liste))-1.7,'LineWidth',4,'marker','x','markersize',10)
hold on;
plot(log10(1./(epsilon_liste)),-1*log10(1./(epsilon_liste))-1.47,'LineWidth',4,'marker','x','markersize',10)
hold on;
% Tracé des lignes verticales
ylims = get(gca,'YLim'); %#ok<NASGU> 
hold on;
legend('CL periodique','avec filtre','pente +1','pente 0','FontSize',25)
xlabel("-log10(eps)",'FontSize',25)
ylabel("log10(Err_{relative} L-inf)",'FontSize',25)
title(sprintf("Evolution de l'erreur relative L-infini(1/3,2/3) en fonction de eps"),'FontSize',25)
set(gca,'FontSize',30)
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%Trace de l'erreur VP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(log10(1./(epsilon_liste)),log10(Liste_Erreur_VP_periodique),'LineWidth',4,'marker','x','markersize',18)
hold on;
plot(log10(1./(epsilon_liste)),log10(Liste_Erreur_VP_filtre),'LineWidth',4,'marker','x','markersize',18)
hold on;
plot(log10(1./(epsilon_liste)),-1*log10(1./(epsilon_liste))-0.9,'LineWidth',4,'marker','x','markersize',18)
hold on;
% plot(log10(1./(epsilon_liste)),-3*log10(1./(epsilon_liste))-0.5,'LineWidth',4,'marker','x','markersize',18)
% hold on;
plot(log10(1./(epsilon_liste)),-2*log10(1./(epsilon_liste))-0.5,'LineWidth',4,'marker','x','markersize',18)
hold on;
% Tracé des lignes verticales
ylims = get(gca,'YLim');
hold on;
legend('CL periodique','avec filtre','pente -1','pente -3','FontSize',25)
xlabel("-log10(eps)",'FontSize',25)
ylabel("log10(Err_{relative} VP)",'FontSize',25)
title(sprintf("Evolution de l'erreur relative VP en fonction de eps"),'FontSize',25)
set(gca,'FontSize',30)
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Trace du multiplicateur Lagrange %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(log10(1./(epsilon_liste)),log10(abs(Liste_multiplicateur_Lagrange_filtre)),'LineWidth',4,'marker','x','markersize',18)
hold on;
plot(log10(1./(epsilon_liste)),-1*log10(1./(epsilon_liste))-0.5,'LineWidth',4,'marker','x','markersize',18)
hold on;
% Tracé des lignes verticales
ylims = get(gca,'YLim');
hold on;
legend('Multiplicateur de Lagrange','pente -1','FontSize',25)
xlabel("-log10(eps)",'FontSize',25)
ylabel("log10(Err_{relative} VP)",'FontSize',25)
title(sprintf("Evolution de l'erreur relative VP en fonction de eps"),'FontSize',25)
set(gca,'FontSize',30)
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(XX,vect_propre_eps_periodique,'LineWidth',4)
hold on;
plot(XX,VV_eps,'LineWidth',4)
hold on;
plot(XX,vect_propre_eps_filtre,'LineWidth',4)
hold on;
legend('Phi_{eps,per}','Psi_{eps}','Phi_{eps,filtre}','FontSize',25)
xlabel("X",'FontSize',25)
ylabel("Vecteurs propres",'FontSize',25)
title(sprintf("1/eps=%g",1/epsilontmp),'FontSize',25)
set(gca,'FontSize',30)
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(vect_propre_eps_periodique_2_celulles,'LineWidth',4)
hold on;
plot(VV_eps_2_cellules,'LineWidth',4)
hold on;
plot(vect_propre_eps_filtre_2_cellules,'LineWidth',4)
hold on;
legend('Phi_{eps,per}','Psi_{eps}','Phi_{eps,filtre}','FontSize',25)
xlabel("X",'FontSize',25)
ylabel("Vecteurs propres",'FontSize',25)
title(sprintf("1/eps=%g",1/epsilontmp),'FontSize',25)
set(gca,'FontSize',30)
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Trace de l'erreur H1 sur 2 cellules %%%%%%%%%%%%%%%%%%%%%
plot(log10(1./(epsilon_liste)),log10(Liste_Erreur_H1_periodique_2cellules),'LineWidth',4,'marker','x','markersize',18)
hold on;
plot(log10(1./(epsilon_liste)),log10(Liste_Erreur_H1_filtre_2cellules),'LineWidth',4,'marker','x','markersize',18)
hold on;
plot(log10(1./(epsilon_liste)),-1*log10(1./(epsilon_liste))+0.11,'LineWidth',4,'marker','x','markersize',18)
hold on;
plot(log10(1./(epsilon_liste)),-2*log10(1./(epsilon_liste))+0.2,'LineWidth',4,'marker','x','markersize',18)
hold on;
% Tracé des lignes verticales
ylims = get(gca,'YLim'); %#ok<NASGU> 
hold on;
legend('CL periodique','avec filtre','pente -1','pente -2','FontSize',25)
xlabel("-log10(eps)",'FontSize',25)
ylabel("log10(Err_{relative} H1)",'FontSize',25)
title(sprintf("Evolution de l'erreur relative H1(0.5-eps;0.5+eps) en fonction de eps"),'FontSize',25)
set(gca,'FontSize',30)
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
