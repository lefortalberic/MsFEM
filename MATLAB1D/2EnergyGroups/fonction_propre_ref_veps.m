function [val_propre_ref,Vecteur_propre_ref_1,Vecteur_propre_ref_2] = fonction_propre_ref_veps(Nbpt_h,epsilon)
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
Psi_1_eps = zeros(Nbpt_h,1); Psi_2_eps = zeros(Nbpt_h,1);
Psi_1_eps_adjoint = zeros(Nbpt_h,1); Psi_2_eps_adjoint = zeros(Nbpt_h,1);
Q_tilde_11_eps = zeros(Nbpt_h,1);
Q_tilde_22_eps = zeros(Nbpt_h,1);
Q_tilde_21_eps = zeros(Nbpt_h,1);
Q_tilde_12_eps = zeros(Nbpt_h,1);
vecteur_J1_eps = zeros(Nbpt_h,1);
vecteur_J2_eps = zeros(Nbpt_h,1);

[~,Psi_1_adjoint,Psi_2_adjoint] = Solution_pb_spectral_adjoint(Nbpt_h);
[~,Psi_1,Psi_2] = Solution_pb_spectral(Nbpt_h);

Q_tilde_11 = Q_tilde(Psi_1,Psi_2,Psi_1_adjoint,Psi_2_adjoint,1);
Q_tilde_22 = Q_tilde(Psi_1,Psi_2,Psi_1_adjoint,Psi_2_adjoint,4);
Q_tilde_21 = - Q_tilde_22;
Q_tilde_12 = - Q_tilde_11;

vecteur_J1 = J(Psi_1,Psi_1_adjoint,1);
vecteur_J2 = J(Psi_2,Psi_2_adjoint,2);

for i=1:Nbpt_h
    xi = (i-1)*h ;
    yi = rem(xi/epsilon,1) ; % y vaut x/epsilon modulo 1
    yk = fix(yi*(Nbpt_h-1)) +1; %partie entière de y*(Nbpt-1) (indice pour les matrices)
    t_g = yi*(Nbpt_h-1) - fix(yi*(Nbpt_h-1)) ;
    Psi_1_eps(i) = Psi_1(yk)*(1-t_g) +Psi_1(yk+1)*t_g;
    Psi_2_eps(i) = Psi_2(yk)*(1-t_g) +Psi_2(yk+1)*t_g;
    Psi_1_eps_adjoint(i) = Psi_1_adjoint(yk)*(1-t_g) +Psi_1_adjoint(yk+1)*t_g;
    Psi_2_eps_adjoint(i) = Psi_2_adjoint(yk)*(1-t_g) +Psi_2_adjoint(yk+1)*t_g;
    Q_tilde_11_eps(i) = Q_tilde_11(yk)*(1-t_g) +Q_tilde_11(yk+1)*t_g;
    Q_tilde_22_eps(i) = Q_tilde_22(yk)*(1-t_g) +Q_tilde_22(yk+1)*t_g;
    Q_tilde_21_eps(i) = Q_tilde_21(yk)*(1-t_g) +Q_tilde_21(yk+1)*t_g;
    Q_tilde_12_eps(i) = Q_tilde_12(yk)*(1-t_g) +Q_tilde_12(yk+1)*t_g;
    vecteur_J1_eps(i) = vecteur_J1(yk)*(1-t_g) +vecteur_J1(yk+1)*t_g;
    vecteur_J2_eps(i) = vecteur_J2(yk)*(1-t_g) +vecteur_J2(yk+1)*t_g;
end

Zeros = sparse(Nbpt_h,Nbpt_h); % matrice de masse

MM_1 = matM_ref_Psi(Nbpt_h, Psi_1_eps , Psi_1_eps_adjoint,h);
MM_2 = matM_ref_Psi(Nbpt_h, Psi_2_eps , Psi_2_eps_adjoint,h);

KK_1 = matK_ref_Psi(Nbpt_h,epsilon, Psi_1_eps , Psi_1_eps_adjoint, 1, 0, 1/(Nbpt_h-1));
KK_2 = matK_ref_Psi(Nbpt_h,epsilon, Psi_2_eps , Psi_2_eps_adjoint, 2, 0, 1/(Nbpt_h-1));

MM_Q_tilde_11 = matM_ref_Psi(Nbpt_h, Q_tilde_11_eps , ones(Nbpt_h,1),h);
MM_Q_tilde_22 = matM_ref_Psi(Nbpt_h, Q_tilde_22_eps , ones(Nbpt_h,1),h);
MM_Q_tilde_21 = matM_ref_Psi(Nbpt_h, Q_tilde_21_eps , ones(Nbpt_h,1),h);
MM_Q_tilde_12 = matM_ref_Psi(Nbpt_h, Q_tilde_12_eps , ones(Nbpt_h,1),h);

CC_1 = matC_ref(Nbpt_h, vecteur_J1_eps);
CC_2 = matC_ref(Nbpt_h, vecteur_J2_eps);

%Assemblage de toutes les matrices
KK = [KK_1 , Zeros ; Zeros ,KK_2 ];
MM_Q_tilde = [MM_Q_tilde_11 , MM_Q_tilde_12 ; MM_Q_tilde_21 , MM_Q_tilde_22];
MM = [MM_1 , Zeros ; Zeros ,MM_2 ];
CC = [CC_1 , Zeros ; Zeros ,CC_2];

%Maillage
XX = zeros(Nbpt_h,1);     %vecteur des coordonnées des noeuds
for j=1:Nbpt_h
    XX(j)=(j-1)*h;
end


% boucle pour les matrices EF
% ------------------------

AA = KK + MM_Q_tilde/(epsilon^2) + CC/epsilon;
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

a1 = norm(Vecteur_propre_ref_1,2)/sqrt(Nbpt_h-1); % Rescale pour une norme L^2 = 1
Vecteur_propre_ref_1 = Vecteur_propre_ref_1./a1;
a2 = norm(Vecteur_propre_ref_2,2)/sqrt(Nbpt_h-1); % Rescale pour une norme L^2 = 1
Vecteur_propre_ref_2 = Vecteur_propre_ref_2./a2;

val_propre_ref = Valeur_propre_min;

end
