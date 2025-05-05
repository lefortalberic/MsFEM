function [vecteur_J] = J(Psi,Psi_adjoint,Num_coeff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% INPUT - x : coordonnee du point ou on veut evaluer la fonction.
% OUTPUT - val: valeur de la matrice sur ce point.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nbpt = size(Psi,1);
Grad_Psi = gradient(Psi,1/(Nbpt-1));
Grad_Psi_adjoint = gradient(Psi_adjoint,1/(Nbpt-1));

vecteur_J = zeros(Nbpt,1);
dy = 1/(Nbpt-1);
for i = 1:Nbpt
    vecteur_J(i) = A((i-1)*dy,1,Num_coeff)*(Psi(i)*Grad_Psi_adjoint(i) - Psi_adjoint(i)*Grad_Psi(i)) ;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
