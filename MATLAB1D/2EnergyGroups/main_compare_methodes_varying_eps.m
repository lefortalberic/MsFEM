
% ---------------------------
%Paramètres
% ---------------------------

%Nbpt_H_list=[3,5,6,11,21,26,51,101,126,201,251,501,1001];
Nbpt_H_list=[11];
H=1/(Nbpt_H_list -1);
%Nbpt_H_list=[5,11,26,51,101,126,201,251];
%Nbpt_H_list = 251;
epsilon = [0.1250,0.10714, 0.0833333333,  0.075, 0.068181, 0.0576923, 0.05,0.04411764,0.03 ,0.0241935,0.01415094,0.0105633 ,0.00572519083];
%epsilon = [0.1250,0.10714, 0.0833333333,  0.075, 0.068181, 0.0576923, 0.05,0.04411764, 0.04166, 0.039473,0.03 ,0.0241935,0.01415094,0.0105633];
%epsilon = [0.1250,0.10714, 0.0833333333,  0.075, 0.068181, 0.0576923, 0.05,0.04411764, 0.04166, 0.039473,0.03 ,0.0241935];

Nbpt_ref = 50001;  %Maillage fin

h_ref = 1/(Nbpt_ref -1);

List_erreur_H1_P1 = [];
List_erreur_H1_triche = [];
List_erreur_H1_triche_bis = [];
List_erreur_H1_MsFEM = [];
List_erreur_VP_P1 = [];
List_erreur_VP_triche = [];
List_erreur_VP_triche_bis = [];
List_erreur_VP_MsFEM = [];

%Solutions MsFEM
for eps=epsilon

    %Solution de référence
    fprintf('\n Calcul de la solution de reference \n');
    tic;
    %[~,Vecteur_propre_ref_1,Vecteur_propre_ref_2] = fonction_propre_ref_veps(Nbpt_ref,eps);
    [val_propre_ref,Vecteur_propre_ref_1,Vecteur_propre_ref_2] = fonction_propre_ref_ueps(Nbpt_ref,eps);
    normL2carre = sqrt(norm(Vecteur_propre_ref_1,2)^2/(Nbpt_ref-1) + norm(Vecteur_propre_ref_2,2)^2/(Nbpt_ref-1)); % Rescale pour une norme L^2 = 1
    Vecteur_propre_ref_1 = Vecteur_propre_ref_1./normL2carre;
    Vecteur_propre_ref_2 = Vecteur_propre_ref_2./normL2carre;

    tempsEcoule = toc;
    disp(['Le temps d''exécution de la solution de reference est de ', num2str(tempsEcoule), ' secondes.']);
    
    [erreur_H1_P1,Erreur_VP_P1] = Liste_Erreur_H1_relative_pb_ueps(6,Nbpt_H_list,eps,Nbpt_ref,Vecteur_propre_ref_1,Vecteur_propre_ref_2,val_propre_ref);
    [erreur_H1_triche,Erreur_VP_triche] = Liste_Erreur_H1_relative_pb_ueps(9,Nbpt_H_list,eps,Nbpt_ref,Vecteur_propre_ref_1,Vecteur_propre_ref_2,val_propre_ref);
    [erreur_H1_triche_bis,Erreur_VP_triche_bis] = Liste_Erreur_H1_relative_vectoriel_pb_ueps(12,Nbpt_H_list,eps,Nbpt_ref,Vecteur_propre_ref_1,Vecteur_propre_ref_2,val_propre_ref);
    %erreur_H1_MsFEM = Liste_Erreur_H1_relative_pb_ueps(11,Nbpt_H_list,eps,Nbpt_ref,Vecteur_propre_ref_1,Vecteur_propre_ref_2);

    List_erreur_H1_P1 = [List_erreur_H1_P1,erreur_H1_P1]; %#ok<*AGROW> 
    List_erreur_H1_triche = [List_erreur_H1_triche,erreur_H1_triche];
    List_erreur_H1_triche_bis = [List_erreur_H1_triche_bis,erreur_H1_triche_bis];

    List_erreur_VP_P1 = [List_erreur_VP_P1,Erreur_VP_P1];
    List_erreur_VP_triche = [List_erreur_VP_triche,Erreur_VP_triche];
    List_erreur_VP_triche_bis = [List_erreur_VP_triche_bis,Erreur_VP_triche_bis];
    %List_erreur_H1_MsFEM = [List_erreur_H1_MsFEM,erreur_H1_MsFEM];

end

% visualisation
% ------------------------------------------------------------------------------

h_list = 1./epsilon;
plot(log10(h_list),log10(min(List_erreur_H1_P1,1)),'LineWidth',4,'marker','x','markersize',28)
hold on;
plot(log10(h_list),log10(List_erreur_H1_triche),'Color','#006400','LineWidth',4,'marker','x','markersize',28)
hold on;
plot(log10(h_list),log10(List_erreur_H1_triche_bis),'Color','red','LineWidth',4,'marker','x','markersize',28)
hold on;
% plot(log10(h_list),log10(List_erreur_H1_MsFEM),'Color','red','LineWidth',4,'marker','x','markersize',28)
% hold on;
% Tracé des lignes verticales
ylims = get(gca,'YLim');
hold on;
% plot([H/(0.03) H/(0.03)],ylims,'--','Color','black','LineWidth',3);
% text(H/(0.03),-0.03, ' \epsilon_{Min,2D}', 'HorizontalAlignment', 'center','FontSize', 28);
% hold on;

%legend('P1 method','Methode MsFEM filtre','FontSize',25)
legend('P1 method','Methode triche','Methode triche bis','FontSize',25)

xlabel("log(1/\epsilon)",'FontSize',25)
ylabel("log(Err_{relative})",'FontSize',25)
title(sprintf("Evolution de l'erreur relative H1 pour le probleme multigroup en fonction de eps / H=%g",H),'FontSize',25)
set(gca,'FontSize',30)
hold off;

% visualisation
% ------------------------------------------------------------------------------

plot(log10(h_list),log10(min(List_erreur_VP_P1,1)),'LineWidth',4,'marker','x','markersize',28)
hold on;
plot(log10(h_list),log10(List_erreur_VP_triche),'Color','#006400','LineWidth',4,'marker','x','markersize',28)
hold on;
plot(log10(h_list),log10(List_erreur_VP_triche_bis),'Color','red','LineWidth',4,'marker','x','markersize',28)
hold on;
% plot(log10(h_list),log10(List_erreur_VP_MsFEM),'Color','red','LineWidth',4,'marker','x','markersize',28)
% hold on;
% Tracé des lignes verticales
ylims = get(gca,'YLim');
hold on;
% plot([H/(0.03) H/(0.03)],ylims,'--','Color','black','LineWidth',3);
% text(H/(0.03),-0.03, ' \epsilon_{Min,2D}', 'HorizontalAlignment', 'center','FontSize', 28);
% hold on;

%legend('P1 method','Methode MsFEM filtre','FontSize',25)
legend('P1 method','Methode triche','Methode triche bis','FontSize',25)

xlabel("log(1/\epsilon)",'FontSize',25)
ylabel("log(Err_{relative})",'FontSize',25)
title(sprintf("Evolution de l'erreur relative VP pour le probleme multigroup en fonction de eps / H=%g",H),'FontSize',25)
set(gca,'FontSize',30)
hold off;