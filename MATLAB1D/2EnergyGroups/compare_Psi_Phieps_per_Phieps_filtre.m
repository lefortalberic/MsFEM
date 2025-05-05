% ------------% ------------% ------------% ------------% ------------
% Code pour comparer Psi_eps et Phi_eps avec methode du filtre.
% Il permet de comparer/tracer différentes estimées d'erreurs sur les
% vecteurs propres et sur les valeurs propres
% ------------% ------------% ------------% ------------% ------------

Nbpt_h = 3001;
%N_Liste = (55:0.01:56);
N_Liste = (3.4:0.4:30);
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
Zeros = sparse(Nbpt_h,Nbpt_h); % matrice de masse

%Maillage
Vect_Un = zeros(Nbpt_ref,1);
PhiXX = zeros(Nbpt_ref,1);     %vecteur des coordonnées des noeuds
for j=1:Nbpt_ref
    PhiXX(j)=Phi((j-1)*h);
    Vect_Un(j) = 1;
end
Norme_Phi = sqrt(PhiXX'*MM_ref*PhiXX);


Psi_1_eps = zeros(Nbpt_h,1); Psi_2_eps = zeros(Nbpt_h,1);
Psi_1_eps_adjoint = zeros(Nbpt_h,1); Psi_2_eps_adjoint = zeros(Nbpt_h,1);

[~,Psi_1_adjoint,Psi_2_adjoint] = Solution_pb_spectral_adjoint(Nbpt_h);
[~,Psi_1,Psi_2] = Solution_pb_spectral(Nbpt_h);

%tracé du coefficient de diffusion
KK_interieur = matK_unitaire((Nbpt_h-1)/3+1,h); % matrice de rigidite
KK_interieur_filtre = sparse((Nbpt_h-1)/3+1,(Nbpt_h-1)/3+1); % matrice de rigidite

for l=(Nbpt_h-1)/3:2*(Nbpt_h-1)/3-1  %La première et dernière ligne doivent rester nulle
    j=l-(Nbpt_h-1)/3 + 1;
    KK_interieur_filtre(j,j)=(1/h)*(Phi((l-1-0.5)*h) + Phi((l-1+0.5)*h)); %approximation de l'intégrale de a(.)
    KK_interieur_filtre(j,j+1)=(-1/h)*Phi((l-1+0.5)*h); %#ok<*SPRIX>
    KK_interieur_filtre(j+1,j)=(-1/h)*Phi((l-1+0.5)*h);
end % for l

KK_interieur_filtre(1,1)=(1/h)*Phi(((Nbpt_h-1)/3 -1+ 0.5)*h);
KK_interieur_filtre((Nbpt_h-1)/3 +1,(Nbpt_h-1)/3 +1)=(1/h)*Phi((2*(Nbpt_h-1)/3 -1 - 0.5)*h);


Liste_Erreur_H1_filtre_1=[];
Liste_Erreur_H1_filtre_2=[];
Liste_Erreur_H1_filtre_1_adjoint=[];
Liste_Erreur_H1_filtre_2_adjoint=[];
Liste_multiplicateur_Lagrange_filtre = [];
Liste_Erreur_H1_filtre_norme_filtre=[];

Liste_Erreur_max_filtre=[];
Liste_Erreur_max_derive_filtre=[];

Liste_Erreur_max_filtre_norme_filtre=[];
Liste_Erreur_max_filtre_norme_racine_filtre=[];

Liste_Erreur_VP_filtre=[];
Liste_Erreur_H1_filtre_2cellules=[];

for epsilon=epsilon_liste
        
        % Phieps_avec_filtrage
        % ------------% ------------% ------------% ------------% ------------
        MM_1 = matM_filtre(Nbpt_h,h,1);
        MM_2 = MM_1;

        KK_1 = matK_filtre(Nbpt_h,epsilon, 0, h, 1, 1);
        KK_2 = matK_filtre(Nbpt_h,epsilon, 0, h, 1, 2);

        MM_sigma_11 = matM_sigma_filtre(Nbpt_h,epsilon, 0, h, 1, 1);
        MM_sigma_12 = matM_sigma_filtre(Nbpt_h,epsilon, 0, h, 1, 2);
        MM_sigma_21 = matM_sigma_filtre(Nbpt_h,epsilon, 0, h, 1, 3);
        MM_sigma_22 = matM_sigma_filtre(Nbpt_h,epsilon, 0, h, 1, 4);

        %Assemblage de toutes les matrices
        KK = [KK_1 , Zeros ; Zeros ,KK_2 ];
        MM_sigma = [MM_sigma_11 , MM_sigma_12 ; MM_sigma_21 , MM_sigma_22];
        MM = [MM_1 , Zeros ; Zeros ,MM_2 ];

        BB = zeros(Nbpt_h,1); %Matrice contrainte périodicité
        for l=2:(Nbpt_h-1)
            BB(l)= Phi((l-1-0.5)*h)-Phi((l-1+0.5)*h);
        end
        BB(1)=-1*Phi((0.5)*h);
        BB(Nbpt_h)=1*Phi(1-(0.5)*h);

        BB_1 = [BB; 0*BB];
        BB_2 = [0*BB; BB];

        AA = MM_sigma + KK*(epsilon*epsilon) ;
        
        AA_f = [AA , BB_1 , BB_2 ; BB_1' , 0 ,0 ; BB_2' , 0 ,0];
        
        MM_f = [MM , 0*BB_1 , 0*BB_2 ; 0*BB_1' , 0 ,0 ; 0*BB_2' , 0 ,0] ;

        [EV,DV_filtre] = eigs(AA_f,MM_f,1,'smallestabs');
        
        [~,arg_min] = min(diag(DV_filtre));
        Vecteur_propre_H_periodic_filtre = EV(:,arg_min);

        Vecteur_propre_filtre_1 = Vecteur_propre_H_periodic_filtre(1:Nbpt_h);
        Vecteur_propre_filtre_2 = Vecteur_propre_H_periodic_filtre((Nbpt_h+1):2*Nbpt_h);
        
        %On a un VP de norme L^2 unitaire. On choisis le VP positif.
        if (Vecteur_propre_filtre_1(fix(Nbpt_h/2))<0)
            Vecteur_propre_filtre_1 = -Vecteur_propre_filtre_1;
        end
        if (Vecteur_propre_filtre_2(fix(Nbpt_h/2))<0)
            Vecteur_propre_filtre_2 = -Vecteur_propre_filtre_2;
        end
        
        a1 = norm(Vecteur_propre_filtre_1,2)/sqrt(Nbpt_h-1); % Rescale pour une norme L^2 = 1
        Vecteur_propre_filtre_1 = Vecteur_propre_filtre_1./a1;
        a2 = norm(Vecteur_propre_filtre_2,2)/sqrt(Nbpt_h-1); % Rescale pour une norme L^2 = 1
        Vecteur_propre_filtre_2 = Vecteur_propre_filtre_2./a2;
        
        Multiplicateur_Lagrange_filtre_1 = Vecteur_propre_H_periodic_filtre(2*Nbpt_h+1);
        Multiplicateur_Lagrange_filtre_2 = Vecteur_propre_H_periodic_filtre(2*Nbpt_h+2);

        %Phieps adjoint _avec_filtrage 
        MM_sigma_transpose = [MM_sigma_11 , MM_sigma_21 ; MM_sigma_12 , MM_sigma_22];
        AA = MM_sigma_transpose + KK*(epsilon*epsilon) ;
        
        AA_f = [AA , BB_1 , BB_2 ; BB_1' , 0 ,0 ; BB_2' , 0 ,0];
        
        MM_f = [MM , 0*BB_1 , 0*BB_2 ; 0*BB_1' , 0 ,0 ; 0*BB_2' , 0 ,0] ;

        [EV,DV_filtre] = eigs(AA_f,MM_f,1,'smallestabs');
        
        [~,arg_min] = min(diag(DV_filtre));
        Vecteur_propre_H_periodic_filtre_adjoint = EV(:,arg_min);

        Vecteur_propre_filtre_1_adjoint = Vecteur_propre_H_periodic_filtre_adjoint(1:Nbpt_h);
        Vecteur_propre_filtre_2_adjoint = Vecteur_propre_H_periodic_filtre_adjoint((Nbpt_h+1):2*Nbpt_h);
        
        %On a un VP de norme L^2 unitaire. On choisis le VP positif.
        if (Vecteur_propre_filtre_1_adjoint(fix(Nbpt_h/2))<0)
            Vecteur_propre_filtre_1_adjoint = -Vecteur_propre_filtre_1_adjoint;
        end
        if (Vecteur_propre_filtre_2_adjoint(fix(Nbpt_h/2))<0)
            Vecteur_propre_filtre_2_adjoint = -Vecteur_propre_filtre_2_adjoint;
        end
        
        a1 = norm(Vecteur_propre_filtre_1_adjoint,2)/sqrt(Nbpt_h-1); % Rescale pour une norme L^2 = 1
        Vecteur_propre_filtre_1_adjoint = Vecteur_propre_filtre_1_adjoint./a1;
        a2 = norm(Vecteur_propre_filtre_2_adjoint,2)/sqrt(Nbpt_h-1); % Rescale pour une norme L^2 = 1
        Vecteur_propre_filtre_2_adjoint = Vecteur_propre_filtre_2_adjoint./a2;
        
        Multiplicateur_Lagrange_filtre_1_adjoint = Vecteur_propre_H_periodic_filtre_adjoint(2*Nbpt_h+1);
        Multiplicateur_Lagrange_filtre_2_adjoint = Vecteur_propre_H_periodic_filtre_adjoint(2*Nbpt_h+2);
        % Psi_eps
        % ------------% ------------% ------------% ------------% ------------
        for i=1:Nbpt_h
            xi = (i-1)*h ;
            yi = rem(xi/epsilon,1) ; % y vaut x/epsilon modulo 1
            yk = fix(yi*(Nbpt_h-1)) +1; %partie entière de y*(Nbpt-1) (indice pour les matrices)
            t_g = yi*(Nbpt_h-1) - fix(yi*(Nbpt_h-1)) ;
            Psi_1_eps(i) = Psi_1(yk)*(1-t_g) +Psi_1(yk+1)*t_g;
            Psi_2_eps(i) = Psi_2(yk)*(1-t_g) +Psi_2(yk+1)*t_g;
            Psi_1_eps_adjoint(i) = Psi_1_adjoint(yk)*(1-t_g) +Psi_1_adjoint(yk+1)*t_g;
            Psi_2_eps_adjoint(i) = Psi_2_adjoint(yk)*(1-t_g) +Psi_2_adjoint(yk+1)*t_g;
        end

        if (Erreur_sur_omega)

            Erreur_H1_Psi_Phi_filtre = sqrt(((VV_eps-vect_propre_eps_filtre)'*KK_ref_A_constant*(VV_eps-vect_propre_eps_filtre))/(VV_eps'*KK_ref_A_constant*VV_eps));
    
            Erreur_H1_Psi_Phi_filtre_norme_filtre = sqrt(((VV_eps-vect_propre_eps_filtre)'*KK_ref_A_constant_filtre*(VV_eps-vect_propre_eps_filtre))/(VV_eps'*KK_ref_A_constant_filtre*VV_eps));

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


        else
            MM_ref_interieur = matM_ref(2*(Nbpt_h-1)/3-(Nbpt_h-1)/3 +1);  
            Psi_1_eps_interieur = Psi_1_eps((Nbpt_h-1)/3:2*(Nbpt_h-1)/3);
            Psi_2_eps_interieur = Psi_2_eps((Nbpt_h-1)/3:2*(Nbpt_h-1)/3);
            Psi_1_eps_adjoint_interieur = Psi_1_eps_adjoint((Nbpt_h-1)/3:2*(Nbpt_h-1)/3);
            Psi_2_eps_adjoint_interieur = Psi_2_eps_adjoint((Nbpt_h-1)/3:2*(Nbpt_h-1)/3);

            Vecteur_propre_filtre_1_interieur = Vecteur_propre_filtre_1((Nbpt_h-1)/3:2*(Nbpt_h-1)/3);
            Vecteur_propre_filtre_2_interieur = Vecteur_propre_filtre_2((Nbpt_h-1)/3:2*(Nbpt_h-1)/3);
            Vecteur_propre_filtre_1_adjoint_interieur = Vecteur_propre_filtre_1_adjoint((Nbpt_h-1)/3:2*(Nbpt_h-1)/3);
            Vecteur_propre_filtre_2_adjoint_interieur = Vecteur_propre_filtre_2_adjoint((Nbpt_h-1)/3:2*(Nbpt_h-1)/3);

            % % Rescale pour une norme L^2 = 1
            Norme_L2_Psi_1_eps_interieur = sqrt(Psi_1_eps_interieur'*MM_ref_interieur*Psi_1_eps_interieur);
            Psi_1_eps_interieur = Psi_1_eps_interieur/Norme_L2_Psi_1_eps_interieur;
            Norme_L2_Psi_2_eps_interieur = sqrt(Psi_2_eps_interieur'*MM_ref_interieur*Psi_2_eps_interieur);
            Psi_2_eps_interieur = Psi_2_eps_interieur/Norme_L2_Psi_2_eps_interieur;
            Norme_L2_Psi_1_eps_adjoint_interieur = sqrt(Psi_1_eps_adjoint_interieur'*MM_ref_interieur*Psi_1_eps_adjoint_interieur);
            Psi_1_eps_adjoint_interieur = Psi_1_eps_adjoint_interieur/Norme_L2_Psi_1_eps_adjoint_interieur;
            Norme_L2_Psi_2_eps_adjoint_interieur = sqrt(Psi_2_eps_adjoint_interieur'*MM_ref_interieur*Psi_2_eps_adjoint_interieur);
            Psi_2_eps_adjoint_interieur = Psi_2_eps_adjoint_interieur/Norme_L2_Psi_2_eps_adjoint_interieur;

            Norme_L2_Vecteur_propre_filtre_1_interieur = sqrt(Vecteur_propre_filtre_1_interieur'*MM_ref_interieur*Vecteur_propre_filtre_1_interieur);
            Vecteur_propre_filtre_1_interieur = Vecteur_propre_filtre_1_interieur/Norme_L2_Vecteur_propre_filtre_1_interieur;
            Norme_L2_Vecteur_propre_filtre_2_interieur = sqrt(Vecteur_propre_filtre_2_interieur'*MM_ref_interieur*Vecteur_propre_filtre_2_interieur);
            Vecteur_propre_filtre_2_interieur = Vecteur_propre_filtre_2_interieur/Norme_L2_Vecteur_propre_filtre_2_interieur;
            Norme_L2_Vecteur_propre_filltre_1_adjoint_interieur = sqrt(Vecteur_propre_filtre_1_adjoint_interieur'*MM_ref_interieur*Vecteur_propre_filtre_1_adjoint_interieur);
            Vecteur_propre_filtre_1_adjoint_interieur = Vecteur_propre_filtre_1_adjoint_interieur/Norme_L2_Vecteur_propre_filltre_1_adjoint_interieur;
            Norme_L2_Vecteur_propre_filtre_2_adjoint_interieur = sqrt(Vecteur_propre_filtre_2_adjoint_interieur'*MM_ref_interieur*Vecteur_propre_filtre_2_adjoint_interieur);
            Vecteur_propre_filtre_2_adjoint_interieur = Vecteur_propre_filtre_2_adjoint_interieur/Norme_L2_Vecteur_propre_filtre_2_adjoint_interieur;


            Erreur_H1_Psi_Phi_filtre_interieur_1 = sqrt(((Psi_1_eps_interieur-Vecteur_propre_filtre_1_interieur)'*KK_interieur*(Psi_1_eps_interieur-Vecteur_propre_filtre_1_interieur))/(Psi_1_eps_interieur'*KK_interieur*Psi_1_eps_interieur));
            Erreur_H1_Psi_Phi_filtre_interieur_2 = sqrt(((Psi_2_eps_interieur-Vecteur_propre_filtre_2_interieur)'*KK_interieur*(Psi_2_eps_interieur-Vecteur_propre_filtre_2_interieur))/(Psi_2_eps_interieur'*KK_interieur*Psi_2_eps_interieur));
            Erreur_H1_Psi_Phi_filtre_interieur_1_adjoint = sqrt(((Psi_1_eps_adjoint_interieur-Vecteur_propre_filtre_1_adjoint_interieur)'*KK_interieur*(Psi_1_eps_adjoint_interieur-Vecteur_propre_filtre_1_adjoint_interieur))/(Psi_1_eps_adjoint_interieur'*KK_interieur*Psi_1_eps_adjoint_interieur));
            Erreur_H1_Psi_Phi_filtre_interieur_2_adjoint = sqrt(((Psi_2_eps_adjoint_interieur-Vecteur_propre_filtre_2_adjoint_interieur)'*KK_interieur*(Psi_2_eps_adjoint_interieur-Vecteur_propre_filtre_2_adjoint_interieur))/(Psi_2_eps_adjoint_interieur'*KK_interieur*Psi_2_eps_adjoint_interieur));

        end
%         Erreur_VP_filtre = abs(DV_filtre-DV_Psi)/abs(DV_Psi);

        Liste_Erreur_H1_filtre_1=[Liste_Erreur_H1_filtre_1,Erreur_H1_Psi_Phi_filtre_interieur_1]; %#ok<AGROW>
        Liste_Erreur_H1_filtre_2=[Liste_Erreur_H1_filtre_2,Erreur_H1_Psi_Phi_filtre_interieur_2]; %#ok<AGROW>
        Liste_Erreur_H1_filtre_1_adjoint=[Liste_Erreur_H1_filtre_1_adjoint,Erreur_H1_Psi_Phi_filtre_interieur_1_adjoint]; %#ok<AGROW>
        Liste_Erreur_H1_filtre_2_adjoint=[Liste_Erreur_H1_filtre_2_adjoint,Erreur_H1_Psi_Phi_filtre_interieur_2_adjoint]; %#ok<AGROW>

end

%%%Trace de l'erreur H1%%%%%%%%%%%%%%%%%%%%%
plot(log10(1./(epsilon_liste)),log10(Liste_Erreur_H1_filtre_1),'LineWidth',4,'marker','x','markersize',18)
hold on;
plot(log10(1./(epsilon_liste)),log10(Liste_Erreur_H1_filtre_2),'LineWidth',4,'marker','x','markersize',18)
hold on;
plot(log10(1./(epsilon_liste)),log10(Liste_Erreur_H1_filtre_1_adjoint),'LineWidth',4,'marker','x','markersize',18)
hold on;
plot(log10(1./(epsilon_liste)),log10(Liste_Erreur_H1_filtre_2_adjoint),'LineWidth',4,'marker','x','markersize',18)
hold on;
plot(log10(1./(epsilon_liste)),-1*log10(1./(epsilon_liste))-0.25,'LineWidth',4,'marker','x','markersize',10)
hold on;
% Tracé des lignes verticales
ylims = get(gca,'YLim');
hold on;
legend('avec filtre Energie1 ','avec filtre Energie2','avec filtre Energie1 adjoint','avec filtre Energie2 adjoint','pente -1','FontSize',25)
xlabel("-log10(eps)",'FontSize',25)
ylabel("log10(Err_{relative} H1)",'FontSize',25)
title(sprintf("Erreur relative H1(1/3,2/3) entre Psi et Phi en fonction de eps"),'FontSize',25)
set(gca,'FontSize',30)
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%Norme L infini%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(log10(1./(epsilon_liste)),log10((Liste_Erreur_max_derive_filtre)),'LineWidth',4,'marker','x','markersize',18)
% hold on;
% plot(log10(1./(epsilon_liste)),-0.5*log10(1./(epsilon_liste))-1.7,'LineWidth',4,'marker','x','markersize',10)
% hold on;
% plot(log10(1./(epsilon_liste)),-1*log10(1./(epsilon_liste))-1.47,'LineWidth',4,'marker','x','markersize',10)
% hold on;
% % Tracé des lignes verticales
% ylims = get(gca,'YLim'); %#ok<NASGU> 
% hold on;
% legend('avec filtre','pente +1','pente 0','FontSize',25)
% xlabel("-log10(eps)",'FontSize',25)
% ylabel("log10(Err_{relative} L-inf)",'FontSize',25)
% title(sprintf("Evolution de l'erreur relative L-infini(1/3,2/3) en fonction de eps"),'FontSize',25)
% set(gca,'FontSize',30)
% hold off;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %%%Trace de l'erreur VP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(log10(1./(epsilon_liste)),log10(Liste_Erreur_VP_filtre),'LineWidth',4,'marker','x','markersize',18)
% hold on;
% plot(log10(1./(epsilon_liste)),-1*log10(1./(epsilon_liste))-0.9,'LineWidth',4,'marker','x','markersize',18)
% hold on;
% plot(log10(1./(epsilon_liste)),-3*log10(1./(epsilon_liste))-0.5,'LineWidth',4,'marker','x','markersize',18)
% hold on;
% % Tracé des lignes verticales
% ylims = get(gca,'YLim');
% hold on;
% legend('avec filtre','pente -1','pente -3','FontSize',25)
% xlabel("-log10(eps)",'FontSize',25)
% ylabel("log10(Err_{relative} VP)",'FontSize',25)
% title(sprintf("Evolution de l'erreur relative VP en fonction de eps"),'FontSize',25)
% set(gca,'FontSize',30)
% hold off;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%Trace du multiplicateur Lagrange %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(log10(1./(epsilon_liste)),log10(abs(Liste_multiplicateur_Lagrange_filtre)),'LineWidth',4,'marker','x','markersize',18)
% hold on;
% plot(log10(1./(epsilon_liste)),-1*log10(1./(epsilon_liste))-0.5,'LineWidth',4,'marker','x','markersize',18)
% hold on;
% % Tracé des lignes verticales
% ylims = get(gca,'YLim');
% hold on;
% legend('Multiplicateur de Lagrange','pente -1','FontSize',25)
% xlabel("-log10(eps)",'FontSize',25)
% ylabel("log10(Err_{relative} VP)",'FontSize',25)
% title(sprintf("Evolution de l'erreur relative VP en fonction de eps"),'FontSize',25)
% set(gca,'FontSize',30)
% hold off;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(XX,VV_eps,'LineWidth',4)
% hold on;
% plot(XX,vect_propre_eps_filtre,'LineWidth',4)
% hold on;
% legend('Phi_{eps,per}','Psi_{eps}','Phi_{eps,filtre}','FontSize',25)
% xlabel("X",'FontSize',25)
% ylabel("Vecteurs propres",'FontSize',25)
% title(sprintf("1/eps=%g",1/epsilontmp),'FontSize',25)
% set(gca,'FontSize',30)
% hold off;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(VV_eps_2_cellules,'LineWidth',4)
% hold on;
% plot(vect_propre_eps_filtre_2_cellules,'LineWidth',4)
% hold on;
% legend('Phi_{eps,per}','Psi_{eps}','Phi_{eps,filtre}','FontSize',25)
% xlabel("X",'FontSize',25)
% ylabel("Vecteurs propres",'FontSize',25)
% title(sprintf("1/eps=%g",1/epsilontmp),'FontSize',25)
% set(gca,'FontSize',30)
% hold off;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%Trace de l'erreur H1 sur 2 cellules %%%%%%%%%%%%%%%%%%%%%
% plot(log10(1./(epsilon_liste)),log10(Liste_Erreur_H1_filtre_2cellules),'LineWidth',4,'marker','x','markersize',18)
% hold on;
% plot(log10(1./(epsilon_liste)),-1*log10(1./(epsilon_liste))+0.11,'LineWidth',4,'marker','x','markersize',18)
% hold on;
% plot(log10(1./(epsilon_liste)),-2*log10(1./(epsilon_liste))+0.2,'LineWidth',4,'marker','x','markersize',18)
% hold on;
% % Tracé des lignes verticales
% ylims = get(gca,'YLim'); %#ok<NASGU> 
% hold on;
% legend('avec filtre','pente -1','pente -2','FontSize',25)
% xlabel("-log10(eps)",'FontSize',25)
% ylabel("log10(Err_{relative} H1)",'FontSize',25)
% title(sprintf("Evolution de l'erreur relative H1(0.5-eps;0.5+eps) en fonction de eps"),'FontSize',25)
% set(gca,'FontSize',30)
% hold off;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
