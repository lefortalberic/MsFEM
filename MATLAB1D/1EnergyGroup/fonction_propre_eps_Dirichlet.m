function [val_propre,vect_propre] = fonction_propre_eps_Dirichlet(Nbpt_H,Nbpt_h,epsilon,noeud_i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pour appeler cette fct il faut faire : [a,b]=fonction_propre_eps(Nbpt_H,Nbpt_h,epsilon,noeud_i)
% fonction de forme solution du
% problème aux valeurs propres

% INPUT - %noeud_i : qui va de 1 à Nbpt_H-1.
%noeud_i=1 correspond à la fct propre entre x=0 et x=H
%noeud_i=Nbpt_H-1 correspond à la fct propre entre x=(Nbpt_H-2)*H (=1-H) et x=1

% OUTPUT - val:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


H=1/(Nbpt_H -1);
h = H/(Nbpt_h -1);

% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
MM = matM_ref(Nbpt_h); % matrice de masse
KK = sparse(Nbpt_h,Nbpt_h); % matrice de rigidite
MM_sigma = sparse(Nbpt_h,Nbpt_h); % matrice de masse

for l=2:(Nbpt_h-2)   %La première et dernière ligne doivent rester nulle

    KK(l,l)=(1/h)*(A((noeud_i-1)*H+(l-1-0.5)*h,epsilon)+A((noeud_i-1)*H+(l-1+0.5)*h,epsilon)); %approximation de l'intégrale de a(.)
    KK(l,l+1)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
    KK(l+1,l)=(-1/h)*A((noeud_i-1)*H+(l-1+0.5)*h,epsilon);

    MM_sigma(l,l)=(h/3)*(Sigma((noeud_i-1)*H+(l-1-0.5)*h,epsilon)+Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon));
    MM_sigma(l,l+1)=h*(1/6)*Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
    MM_sigma(l+1,l)=h*(1/6)*Sigma((noeud_i-1)*H+(l-1+0.5)*h,epsilon);
            
end % for l

MM_sigma(Nbpt_h-1,Nbpt_h-1)=(h/3)*(Sigma((noeud_i-1)*H+(Nbpt_h-2-0.5)*h,epsilon)+Sigma((noeud_i-1)*H+(Nbpt_h-2+0.5)*h,epsilon));
KK(Nbpt_h-1,Nbpt_h-1)=(1/h)*(A((noeud_i-1)*H+(Nbpt_h-2-0.5)*h,epsilon)+A((noeud_i-1)*H+(Nbpt_h-2+0.5)*h,epsilon));

%Maillage
XX = zeros(Nbpt_h,1);     %vecteur des coordonnées des noeuds
for j=1:Nbpt_h
    XX(j)=(j-1)*h;
end


% boucle pour les matrices EF
% ------------------------

AA = MM_sigma + KK*(epsilon^2) ;
AA(1,1)=1; AA(Nbpt_h,Nbpt_h)=1;

%[EV,DV] = eigs(AA,MM,5,'smallestabs');
[EV,DV] = eigs(AA,MM,1,'smallestabs');


[Valeur_propre_min,arg_min] = min(diag(DV));
Vecteur_propre = EV(:,arg_min);
Vecteur_propre = Vecteur_propre*(sqrt(H)); %On "renormalise" pour que le vecteur soit toujours aux alentours de 1

%On a un VP de norme L^2 unitaire. On choisis le VP positif.
if (Vecteur_propre(fix(Nbpt_h/2))<0)
    Vecteur_propre = -Vecteur_propre;
end

a = norm(Vecteur_propre,2)/sqrt(Nbpt_h-1); % Rescale pour une norme L^2 = 1
Vecteur_propre = Vecteur_propre./a;

val_propre = Valeur_propre_min;
vect_propre = Vecteur_propre;
%plot(vect_propre)

end
