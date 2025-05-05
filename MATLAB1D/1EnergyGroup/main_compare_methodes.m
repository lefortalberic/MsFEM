
% ---------------------------
%Paramètres
% ---------------------------

%Nbpt_H_list=[3,5,6,11,21,26,51,101,126,201,251,501,1001];
%Nbpt_H_list=[9,11,21,26,51,101,126,201,251];
Nbpt_H_list=[5,9,11,26,51,101];
%Nbpt_H_list = 251;
epsilon = 0.02;

Nbpt_ref = 20001;  %Maillage fin

h_ref = 1/(Nbpt_ref -1);

%Solution de référence
fprintf('\n Calcul de la solution de reference \n');
tic;
[~,Vecteur_propre_ref] = fonction_propre_eps_Dirichlet(2,Nbpt_ref,epsilon,1);
tempsEcoule = toc;
disp(['Le temps d''exécution de la solution de reference est de ', num2str(tempsEcoule), ' secondes.']);

%Solutions MsFEM
List_erreur_H1_P1 = Liste_Erreur_H1_relative(6,Nbpt_H_list,epsilon,Nbpt_ref,Vecteur_propre_ref);
List_erreur_H1_triche = Liste_Erreur_H1_relative(8,Nbpt_H_list,epsilon,Nbpt_ref,Vecteur_propre_ref);
List_erreur_H1_0_3H = Liste_Erreur_H1_relative(11,Nbpt_H_list,epsilon,Nbpt_ref,Vecteur_propre_ref);
%List_erreur_H1_0_3H_moyenne = Liste_Erreur_H1_relative(17,Nbpt_H_list,epsilon,Nbpt_ref,Vecteur_propre_ref);
List_erreur_H1_0_3H_periodique = Liste_Erreur_H1_relative(19,Nbpt_H_list,epsilon,Nbpt_ref,Vecteur_propre_ref);
List_erreur_H1_0_9H_periodique = Liste_Erreur_H1_relative(21,Nbpt_H_list,eps,Nbpt_ref,Vecteur_propre_ref);


% visualisation
% ------------------------------------------------------------------------------

h_list = 1./(Nbpt_H_list -1);
plot(-log10(h_list),log10(List_erreur_H1_P1),'LineWidth',4,'marker','x','markersize',28)
hold on;
plot(-log10(h_list),log10(List_erreur_H1_triche),'LineWidth',4,'marker','x','markersize',28)
hold on;
plot(-log10(h_list),log10(List_erreur_H1_0_3H),'LineWidth',4,'marker','x','markersize',28)
hold on;
% plot(-log10(h_list),log10(List_erreur_H1_0_3H_moyenne),'LineWidth',4,'marker','x','markersize',28)
% hold on;
plot(-log10(h_list),log10(List_erreur_H1_0_3H_periodique),'LineWidth',4,'marker','x','markersize',28)
hold on;
plot(-log10(h_list),log10(List_erreur_H1_0_9H_periodique),'LineWidth',4,'marker','x','markersize',28)
hold on;
%plot(-log10(h_list),log10(h_list.^1)-0.15,'LineWidth',4)
hold on;

%legend('Methode P1','Methode triche','Methode MsFEM\_0\_3H','Methode MsFEM\_0\_3H\_moyenne','Methode MsFEM\_0\_3H\_sinus','FontSize',25)
legend('P1 method','Cheating method','MsFEM dirichlet method','MsFEM periodique 3H','MsFEM periodique 9H','FontSize',25)

xlabel("log(1/H)",'FontSize',25)
ylabel("log(Err_{relative})",'FontSize',25)
%title(sprintf("Evolution de l'erreur relative H1 en fonction de H / eps=%g",epsilon),'FontSize',25)
set(gca,'FontSize',30)
hold off;