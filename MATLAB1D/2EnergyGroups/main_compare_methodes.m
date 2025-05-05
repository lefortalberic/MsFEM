
% ---------------------------
%Paramètres
% ---------------------------

%Nbpt_H_list=[3,5,6,11,21,26,51,101,126,201,251,501,1001];
%Nbpt_H_list=[9,11,21,26,51,101,126,201,251];
Nbpt_H_list=[5,9,11,26,51,101,201];
%Nbpt_H_list = 201;
epsilon = 0.024;

Nbpt_ref = 5001;  %Maillage fin

h_ref = 1/(Nbpt_ref -1);

%Solution de référence
fprintf('\n Calcul de la solution de reference \n');
tic;
[~,Vecteur_propre_ref_1,Vecteur_propre_ref_2] = fonction_propre_ref_veps(Nbpt_ref,epsilon);
tempsEcoule = toc;
disp(['Le temps d''exécution de la solution de reference est de ', num2str(tempsEcoule), ' secondes.']);

%Solutions MsFEM
List_erreur_H1_P1 = Liste_Erreur_H1_relative_pb_veps(6,Nbpt_H_list,epsilon,Nbpt_ref,Vecteur_propre_ref_1,Vecteur_propre_ref_2);
List_erreur_H1_triche = Liste_Erreur_H1_relative_pb_veps(8,Nbpt_H_list,epsilon,Nbpt_ref,Vecteur_propre_ref_1,Vecteur_propre_ref_2);
List_erreur_H1_MsFEM = Liste_Erreur_H1_relative_pb_veps(10,Nbpt_H_list,epsilon,Nbpt_ref,Vecteur_propre_ref_1,Vecteur_propre_ref_2);

% visualisation
% ------------------------------------------------------------------------------

H_list = 1./(Nbpt_H_list -1);
plot(-log10(H_list),log10(List_erreur_H1_P1),'LineWidth',4,'marker','x','markersize',28)
hold on;
plot(-log10(H_list),log10(List_erreur_H1_triche),'Color','#006400','LineWidth',4,'marker','x','markersize',28)
hold on;
plot(-log10(H_list),log10(List_erreur_H1_MsFEM),'Color','red','LineWidth',4,'marker','x','markersize',28)
hold on;

% Tracé des lignes verticales
ylims = get(gca,'YLim');
hold on;
plot([-log10(epsilon) -log10(epsilon)],ylims,'--','Color','black','LineWidth',3);
text(-log10(epsilon),-0.5, 'H = \epsilon', 'HorizontalAlignment', 'center','FontSize', 38);
hold on;

legend('P1 method','Cheating method','MsFEM filtre','FontSize',25)

xlabel("log(1/H)",'FontSize',25)
ylabel("log(Err_{relative})",'FontSize',25)
title(sprintf("Erreur relative H1 pour le probleme multigroup en veps en fonction de H / eps=%g",epsilon),'FontSize',25)
set(gca,'FontSize',30)
hold off;