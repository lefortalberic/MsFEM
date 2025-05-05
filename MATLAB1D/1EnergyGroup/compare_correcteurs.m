% ------------% ------------% ------------% ------------% ------------
% Code Qui va comparer le correcteur (du pb purement diffusif de l'article
% de CLB + Xavier Blanc)
% défini sur la cellule et évalué en 
% x/eps au correcteur défini sur 0,1 mais avec A_eps ou 1/eps n'est pas 
% entier
%
%
% ------------% ------------% ------------% ------------% ------------

Nbpt_h = 30001;
N_Liste = (1.5:1:20.5);
epsilon_liste =1./N_Liste;
%epsilon_liste = 1/9.5;
%epsilon_liste = [1/1.5 ,1/3, 1/5.5 ,1/10.5,1/15, 1/20.5,1/35,1/50.5,1/100.5];
taille_domaine = 1;

h = 1/(Nbpt_h -1);
Nbpt_ref = (Nbpt_h-1)+1;
Nbpt_per = Nbpt_h-1;
h_ref = 1/(Nbpt_ref-1);
%Maillage
XX = zeros(Nbpt_h,1);     %vecteur des coordonnées des noeuds
for j=1:Nbpt_ref
    XX(j)=(j-1)*h;
end
PhiXX = zeros(Nbpt_ref,1);     %vecteur des coordonnées des noeuds
for j=1:Nbpt_ref
    PhiXX(j)=Phi((j-1)*h);
end

%%% Correcteur sur la cellule %%%
% ----------------------
KK = matK_per(Nbpt_per); % matrice de rigidite avec A remplacé par Psi**2*A

H = zeros(Nbpt_per,1);
for i=2:(Nbpt_per)
    H(i)= A((i-2+0.5)*h,1) - A((i-1+0.5)*h,1);
end
H(1)= -A(0.5*h,1) + A((Nbpt_per-0.5)*h,1);

AA = KK ;

FF = -H;
warning('off', 'all');
W = AA\FF;
warning('on', 'all');
W = W - mean(W);
W = [W ; W(1)];

% ------------% ------------% ------------% ------------% ------------
%%% Correcteur approchée dépendant de eps %%%
KK_ref_A_constant = matK_unitaire((1*(Nbpt_h-1)+1),h_ref);

KK_ref_A_constant_filtre = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
for l=1:(Nbpt_h-1)   %La première et dernière ligne doivent rester nulle
    KK_ref_A_constant_filtre(l,l)=(1/h)*(2)*Phi((l-1)*h); %approximation de l'intégrale de a(.)
    KK_ref_A_constant_filtre(l,l+1)=(-1/h)*Phi((l-1+0.5)*h); %#ok<*SPRIX>
    KK_ref_A_constant_filtre(l+1,l)=(-1/h)*Phi((l-1+0.5)*h);
end % for l

KK_ref_A_constant_filtre(1,1)=(1/h)*Phi((0.5)*h);
KK_ref_A_constant_filtre(Nbpt_h,Nbpt_h)=(1/h)*Phi(1-(0.5)*h);

Liste_Erreur_H1_Wcell_WR=[];
Liste_Erreur_H1_Wcell_WRfiltre=[];
Liste_Erreur_max_Wcell_WR=[];
Liste_Erreur_max_Wcell_WRfiltre=[];
Liste_Erreur_H1_Wcell_WRfiltre_norme_sans_filtre=[];
Liste_Erreur_A_homog_Wcell_WR=[];
Liste_Erreur_A_homog_Wcell_WRfiltre=[];
for epsilon=epsilon_liste

        %%% Correcteur_eps^R (0,1) %%%
        % ----------------------
        KK = sparse(Nbpt_per,Nbpt_per); % matrice de rigidite
        for l=1:(Nbpt_per-1)   %La première et dernière ligne doivent rester nulle
            KK(l,l)=(1/h)*(A((l-1-0.5)*h,epsilon)+A((l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
            KK(l,l+1)=(-1/h)*A((l-1+0.5)*h,epsilon); %#ok<*SPRIX>
            KK(l+1,l)=(-1/h)*A((l-1+0.5)*h,epsilon);
        end % for l
        
        KK(1,1)=(1/h)*(A(0.5*h,epsilon)+A(1-0.5*h,epsilon));
        KK(1,Nbpt_per)=(-1/h)*A(1-0.5*h,epsilon);
        KK(Nbpt_per,1)=(-1/h)*A(1-0.5*h,epsilon);
        KK(Nbpt_per,Nbpt_per)=(1/h)*(A((Nbpt_per-1-0.5)*h,epsilon)+A((Nbpt_per-1+0.5)*h,epsilon));

        H = zeros(Nbpt_per,1);
        for i=2:(Nbpt_per)
            H(i)= A((i-2+0.5)*h,epsilon) - A((i-1+0.5)*h,epsilon);
        end
        H(1)= -A(0.5*h,epsilon) + A((Nbpt_per-0.5)*h,epsilon);
        
        AA = epsilon*KK ;
        
        FF = -H;
        warning('off', 'all');
        WR = AA\FF;
        warning('on', 'all');
        WR = WR - mean(WR);
        WR = [WR ; WR(1)]; %#ok<AGROW>


        %%% Correcteur_eps^R_Phi avec filtre (0,1) %%%
        % ----------------------
        KK = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
        for l=1:(Nbpt_h-1)   %La première et dernière ligne doivent rester nulle
            KK(l,l)=(1/h)*(A((l-1-0.5)*h,epsilon)*Phi((l-1-0.5)*h)+A((l-1+0.5)*h,epsilon)*Phi((l-1+0.5)*h)); %approximation de l'intégrale de a(.)
            KK(l,l+1)=(-1/h)*A((l-1+0.5)*h,epsilon)*Phi((l-1+0.5)*h); %#ok<*SPRIX>
            KK(l+1,l)=(-1/h)*A((l-1+0.5)*h,epsilon)*Phi((l-1+0.5)*h);
        end % for l
        
        KK(1,1)=(1/h)*(A(0.5*h,epsilon))*Phi((0.5)*h);
        KK(Nbpt_h,Nbpt_h)=(1/h)*(A(1-(0.5)*h,epsilon))*Phi(1-(0.5)*h);

        H = zeros(Nbpt_h,1);
        for i=2:(Nbpt_h-1)
            H(i)= A((i-2+0.5)*h,epsilon)*Phi((i-2+0.5)*h) - A((i-1+0.5)*h,epsilon)*Phi((i-1+0.5)*h);
        end
        H(1)= -A(0.5*h,epsilon)*Phi((0.5)*h);
        H(Nbpt_h)= A((Nbpt_h-0.5)*h,epsilon)*Phi((Nbpt_h-0.5)*h);

        BB = zeros(Nbpt_h,1); %Matrice contrainte périodicité
        
        for l=2:(Nbpt_h-1)
            BB(l)= Phi((l-1-0.5)*h)-Phi((l-1+0.5)*h);
        end
        BB(1)=-1*Phi((0.5)*h);
        BB(Nbpt_h)=1*Phi(1-(0.5)*h);

        AA = epsilon*[KK , BB ; BB' , 0];
        
        FF = [-H ; 0];

        warning('off', 'all');
        WRfiltre = AA\FF;
        warning('on', 'all');

        WRfiltre(Nbpt_h+1)

        WRfiltre = WRfiltre(1:Nbpt_h);
        WRfiltre = WRfiltre - mean(WRfiltre);
        %%%%%%%%%%%%%%%%%%%%%
        W_eps = zeros(Nbpt_h,1);
        for i=1:Nbpt_per+1
            xi_g = (i-1)*h ;
            yi_g = rem(xi_g/epsilon,1) ; % y vaut x/epsilon modulo 1
            yk_g = fix(yi_g*(Nbpt_h-1)) +1; %partie entière de y*(Nbpt-1) (indice pour les matrices)
        
            t_g = yi_g*(Nbpt_h-1) - fix(yi_g*(Nbpt_h-1)) ;
            W_eps(i) = W(yk_g)*(1-t_g) +W(yk_g+1)*t_g;
        end
        
        grad_W_eps_vect = zeros(Nbpt_h-1,1);
        grad_W_filtre_vect = zeros(Nbpt_h-1,1);
        grad_W_R_vect = zeros(Nbpt_h-1,1);
        %Calcul de la matrice homogénéisée
        A_homog_cell = 0; A_homog_R = 0; A_homog_R_filtre = 0; A_homog_cell_eps = 0; test = zeros(Nbpt_h,1);
        for i=1:(Nbpt_h-1)
            grad_W_cell = (W(i+1) - W(i))/h;
            A_homog_cell=A_homog_cell + h*A((i-1+0.5)*h,1)*(1+grad_W_cell); %formulation avec les \tilde w_i

            grad_W_eps = (W_eps(i+1) - W_eps(i))/h;
            A_homog_cell_eps=A_homog_cell_eps + h*A((i-1+0.5)*h,epsilon)*(1+epsilon*grad_W_eps); %formulation avec les \tilde w_i

            grad_WR = (WR(i+1) - WR(i))/h;
            A_homog_R=A_homog_R + h*A((i-1+0.5)*h,epsilon)*(1+epsilon*grad_WR); %formulation avec les \tilde w_i

            grad_WRfiltre = (WRfiltre(i+1) - WRfiltre(i))/h;
            A_homog_R_filtre=A_homog_R_filtre + h*A((i-1+0.5)*h,epsilon)*(1+epsilon*grad_WRfiltre)*Phi((i-1+0.5)*h); %formulation avec les \tilde w_i


            grad_W_eps_vect(i)=(W_eps(i+1) - W_eps(i))/h;
            grad_W_R_vect(i) = (WR(i+1) - WR(i))/h;
            grad_W_filtre_vect(i) = (WRfiltre(i+1) - WRfiltre(i))/h;
            test(i) = A((i-1+0.5)*h,epsilon)*(1+epsilon*grad_WRfiltre);
        end

        A_homog_R_filtre = A_homog_R_filtre/(sum(PhiXX)/Nbpt_h);

        KK_interieur = KK_ref_A_constant((Nbpt_h-1)/3:2*(Nbpt_h-1)/3,(Nbpt_h-1)/3:2*(Nbpt_h-1)/3);
        KK_interieur(1,1)=KK_interieur(1,1)/2;
        KK_interieur((Nbpt_h-1)/3 +1,(Nbpt_h-1)/3 +1) = KK_interieur((Nbpt_h-1)/3 +1,(Nbpt_h-1)/3 +1)/2;

        KK_interieur_filtre = KK_ref_A_constant_filtre((Nbpt_h-1)/3:2*(Nbpt_h-1)/3,(Nbpt_h-1)/3:2*(Nbpt_h-1)/3);
        KK_interieur_filtre(1,1)=KK_interieur_filtre(1,1)/2;
        KK_interieur_filtre((Nbpt_h-1)/3 +1,(Nbpt_h-1)/3 +1) = KK_interieur_filtre((Nbpt_h-1)/3 +1,(Nbpt_h-1)/3 +1)/2;


        Erreur_H1_Wcell_WR = sqrt(((W_eps((Nbpt_h-1)/3:2*(Nbpt_h-1)/3)-WR((Nbpt_h-1)/3:2*(Nbpt_h-1)/3))'*KK_interieur*(W_eps((Nbpt_h-1)/3:2*(Nbpt_h-1)/3)-WR((Nbpt_h-1)/3:2*(Nbpt_h-1)/3)))/(W_eps((Nbpt_h-1)/3:2*(Nbpt_h-1)/3)'*KK_interieur*W_eps((Nbpt_h-1)/3:2*(Nbpt_h-1)/3)));
        Erreur_H1_Wcell_WRfiltre = sqrt(((W_eps((Nbpt_h-1)/3:2*(Nbpt_h-1)/3)-WRfiltre((Nbpt_h-1)/3:2*(Nbpt_h-1)/3))'*KK_interieur_filtre*(W_eps((Nbpt_h-1)/3:2*(Nbpt_h-1)/3)-WRfiltre((Nbpt_h-1)/3:2*(Nbpt_h-1)/3)))/(W_eps((Nbpt_h-1)/3:2*(Nbpt_h-1)/3)'*KK_interieur_filtre*W_eps((Nbpt_h-1)/3:2*(Nbpt_h-1)/3)));
        Erreur_H1_Wcell_WRfiltre_norme_sans_filtre = sqrt(((W_eps((Nbpt_h-1)/3:2*(Nbpt_h-1)/3)-WRfiltre((Nbpt_h-1)/3:2*(Nbpt_h-1)/3))'*KK_interieur*(W_eps((Nbpt_h-1)/3:2*(Nbpt_h-1)/3)-WRfiltre((Nbpt_h-1)/3:2*(Nbpt_h-1)/3)))/(W_eps((Nbpt_h-1)/3:2*(Nbpt_h-1)/3)'*KK_interieur*W_eps((Nbpt_h-1)/3:2*(Nbpt_h-1)/3)));

        Erreur_max_Wcell_WR = max(abs(W_eps((Nbpt_h-1)/3:2*(Nbpt_h-1)/3)-WR((Nbpt_h-1)/3:2*(Nbpt_h-1)/3)))/max(abs(W_eps((Nbpt_h-1)/3:2*(Nbpt_h-1)/3)));
        Erreur_max_Wcell_WR_filtre = max(abs(W_eps((Nbpt_h-1)/3:2*(Nbpt_h-1)/3)-WRfiltre((Nbpt_h-1)/3:2*(Nbpt_h-1)/3)))/max(abs(W_eps((Nbpt_h-1)/3:2*(Nbpt_h-1)/3)));
        Liste_Erreur_max_Wcell_WR=[Liste_Erreur_max_Wcell_WR,Erreur_max_Wcell_WR]; %#ok<AGROW> 
        Liste_Erreur_max_Wcell_WRfiltre=[Liste_Erreur_max_Wcell_WRfiltre,Erreur_max_Wcell_WR_filtre]; %#ok<AGROW> 

        Liste_Erreur_H1_Wcell_WR=[Liste_Erreur_H1_Wcell_WR,Erreur_H1_Wcell_WR]; %#ok<AGROW> 
        Liste_Erreur_H1_Wcell_WRfiltre=[Liste_Erreur_H1_Wcell_WRfiltre,Erreur_H1_Wcell_WRfiltre]; %#ok<AGROW> 
        Liste_Erreur_H1_Wcell_WRfiltre_norme_sans_filtre=[Liste_Erreur_H1_Wcell_WRfiltre_norme_sans_filtre,Erreur_H1_Wcell_WRfiltre_norme_sans_filtre]; %#ok<AGROW> 

        Liste_Erreur_A_homog_Wcell_WR=[Liste_Erreur_A_homog_Wcell_WR, abs(A_homog_cell-A_homog_R)]; %#ok<AGROW> 
        Liste_Erreur_A_homog_Wcell_WRfiltre=[Liste_Erreur_A_homog_Wcell_WRfiltre,abs(A_homog_cell-A_homog_R_filtre)]; %#ok<AGROW> 
end


plot(XX,W_eps,'LineWidth',4)
hold on;
plot(XX,WR,'LineWidth',4)
hold on;
plot(XX,WRfiltre,'LineWidth',4)
hold on;
legend('W_cell(x/\epsilon)','W_R','W_R filtre','FontSize',25)
xlabel("X",'FontSize',25)
ylabel("Vecteurs propres",'FontSize',25)
title(sprintf("1/eps=%g",1/epsilon),'FontSize',25)
set(gca,'FontSize',30)
hold off;


plot(grad_W_eps_vect,'LineWidth',4)
hold on;
plot(grad_W_R_vect,'LineWidth',4)
hold on;
plot(grad_W_filtre_vect,'LineWidth',4)
hold on;
legend('W_cell(x/\epsilon)','W_R sans filtre','W_R avec filtre','FontSize',25)
xlabel("X",'FontSize',25)
title(sprintf("1/eps=%g",1/epsilon),'FontSize',25)
set(gca,'FontSize',30)
hold off;


%%%Coeffs A_homog
plot(-log10(epsilon_liste),log10(Liste_Erreur_A_homog_Wcell_WR),'LineWidth',4,'marker','x','markersize',18)
hold on;
plot(-log10(epsilon_liste),log10(Liste_Erreur_A_homog_Wcell_WRfiltre),'LineWidth',4,'marker','x','markersize',18)
hold on;
% Tracé des lignes verticales
ylims = get(gca,'YLim');
hold on;
legend('Sans filtre','Avec filtre','FontSize',25)
xlabel("1/eps",'FontSize',25)
ylabel("Err_{relative} H1",'FontSize',25)
title(sprintf("Evolution de l'erreur relative H1 en fonction de eps"),'FontSize',25)
set(gca,'FontSize',30)
hold off;

%%%Norme L infini
plot(-log10(epsilon_liste),log10(Liste_Erreur_max_Wcell_WR),'LineWidth',4,'marker','x','markersize',18)
hold on;
plot(-log10(epsilon_liste),log10(Liste_Erreur_max_Wcell_WRfiltre),'LineWidth',4,'marker','x','markersize',18)
hold on;
% Tracé des lignes verticales
ylims = get(gca,'YLim');
hold on;
legend('Sans filtre','Avec filtre','FontSize',25)
xlabel("1/eps",'FontSize',25)
ylabel("Err_{relative} L infini",'FontSize',25)
title(sprintf("Evolution de l'erreur relative H1 en fonction de eps"),'FontSize',25)
set(gca,'FontSize',30)
hold off;

%%%Erreurs H1
plot(-log10(epsilon_liste),log10(Liste_Erreur_H1_Wcell_WR),'LineWidth',4,'marker','x','markersize',18)
hold on;
plot(-log10(epsilon_liste),log10(Liste_Erreur_H1_Wcell_WRfiltre),'LineWidth',4,'marker','x','markersize',18)
hold on;
plot(-log10(epsilon_liste),log10(Liste_Erreur_H1_Wcell_WRfiltre_norme_sans_filtre),'LineWidth',4,'marker','x','markersize',18)
hold on;
% Tracé des lignes verticales
ylims = get(gca,'YLim');
hold on;
legend('Sans filtre','Avec filtre','avecfiltre norme classique','FontSize',25)
xlabel("1/eps",'FontSize',25)
ylabel("Err_{relative} H1",'FontSize',25)
title(sprintf("Evolution de l'erreur relative H1 en fonction de eps"),'FontSize',25)
set(gca,'FontSize',30)
hold off;
