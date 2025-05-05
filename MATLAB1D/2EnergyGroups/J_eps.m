function [vecteur_J_eps] = J_eps(h, epsilon, x_gauche, Phi,Phi_adjoint,Num_coeff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% INPUT - x : coordonnee du point ou on veut evaluer la fonction.
% OUTPUT - val: valeur de la matrice sur ce point.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nbpt = size(Phi,1);
Grad_Phi = gradient(Phi,h);
Grad_Phi_adjoint = gradient(Phi_adjoint,h);

vecteur_J_eps = zeros(Nbpt,1);

for i = 1:Nbpt
    vecteur_J_eps(i) = A(x_gauche + (i-1)*h,epsilon,Num_coeff)*(Phi(i)*Grad_Phi_adjoint(i) - Phi_adjoint(i)*Grad_Phi(i)) ;
end
vecteur_J_eps = epsilon *vecteur_J_eps;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
