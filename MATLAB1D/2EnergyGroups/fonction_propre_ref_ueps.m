function [val_propre_ref,Vecteur_propre_ref_1,Vecteur_propre_ref_2] = fonction_propre_ref_ueps(Nbpt_h,epsilon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pour appeler cette fct il faut faire : [a,b]=fonction_propre_eps(Nbpt_H,Nbpt_h,epsilon,noeud_i)
% fonction de forme solution du
% problème aux valeurs propres

% INPUT - %noeud_i : qui va de 1 à Nbpt_H-1.
%noeud_i=1 correspond à la fct propre entre x=0 et x=H
%noeud_i=Nbpt_H-1 correspond à la fct propre entre x=(Nbpt_H-2)*H (=1-H) et x=1

% OUTPUT - val:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = 1/(Nbpt_h -1);

% ----------------------
% calcul des matrices EF
% ----------------------

Zeros = sparse(Nbpt_h,Nbpt_h); % matrice de masse

MM_1 = matM_ref(Nbpt_h);
MM_2 = MM_1;

KK_1 = matK_ref(Nbpt_h,epsilon, 1);
KK_2 = matK_ref(Nbpt_h,epsilon, 2);

MM_sigma_11 = matM_sigma_ref(Nbpt_h,epsilon, 1);
MM_sigma_22 = matM_sigma_ref(Nbpt_h,epsilon, 4);
MM_sigma_21 = matM_sigma_ref(Nbpt_h,epsilon, 3);
MM_sigma_12 = matM_sigma_ref(Nbpt_h,epsilon, 2);


%Assemblage de toutes les matrices
KK = [KK_1 , Zeros ; Zeros ,KK_2 ];
MM_sigma = [MM_sigma_11 , MM_sigma_12 ; MM_sigma_21 , MM_sigma_22];
MM = [MM_1 , Zeros ; Zeros ,MM_2 ];

%Maillage
XX = zeros(Nbpt_h,1);     %vecteur des coordonnées des noeuds
for j=1:Nbpt_h
    XX(j)=(j-1)*h;
end


% boucle pour les matrices EF
% ------------------------

AA = KK*(epsilon^2) + MM_sigma ;
AA(1,1)=1; AA(2*Nbpt_h,2*Nbpt_h)=1;
AA(Nbpt_h,Nbpt_h)=1 ; AA(Nbpt_h+1,Nbpt_h+1)=1;

%[EV,DV] = eigs(AA,MM,5,'smallestabs');
[EV,DV] = eigs(AA,MM,1,'smallestabs');


[Valeur_propre_min,arg_min] = min(diag(DV));
Vecteur_propre = EV(:,arg_min);

Vecteur_propre_ref_1 = Vecteur_propre(1:Nbpt_h);
Vecteur_propre_ref_2 = Vecteur_propre((Nbpt_h+1):2*Nbpt_h);

%On a un VP de norme L^2 unitaire. On choisis le VP positif.
if (Vecteur_propre_ref_1(fix(Nbpt_h/2))<0)
    Vecteur_propre_ref_1 = -Vecteur_propre_ref_1;
end
if (Vecteur_propre_ref_2(fix(Nbpt_h/2))<0)
    Vecteur_propre_ref_2 = -Vecteur_propre_ref_2;
end

% a1 = norm(Vecteur_propre_ref_1,2)/sqrt(Nbpt_h-1); % Rescale pour une norme L^2 = 1
% Vecteur_propre_ref_1 = Vecteur_propre_ref_1./a1;
% a2 = norm(Vecteur_propre_ref_2,2)/sqrt(Nbpt_h-1); % Rescale pour une norme L^2 = 1
% Vecteur_propre_ref_2 = Vecteur_propre_ref_2./a2;

val_propre_ref = Valeur_propre_min;

end
