function [val] = check_symmetry(Psi_1,Psi_2,Psi_1_adjoint,Psi_2_adjoint)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% INPUT - x : coordonnee du point ou on veut evaluer la fonction.
% OUTPUT - val: valeur de la matrice sur ce point.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nbpt = size(Psi_1,1);
Grad_Psi_1 = gradient(Psi_1,1/(Nbpt-1));
Grad_Psi_2 = gradient(Psi_2,1/(Nbpt-1));
Grad_Psi_1_adjoint = gradient(Psi_1_adjoint,1/(Nbpt-1));
Grad_Psi_2_adjoint = gradient(Psi_2_adjoint,1/(Nbpt-1));

val_1 = 0;
val_2 = 0;
dy = 1/(Nbpt-1);
for i = 1:(Nbpt-1)
    val_1 = val_1 + A((i-1)*dy,1,1)*(Psi_1(i)*Grad_Psi_1_adjoint(i) - Psi_1_adjoint(i)*Grad_Psi_1(i))*dy;
    val_2 = val_2 + A((i-1)*dy,1,2)*(Psi_2(i)*Grad_Psi_2_adjoint(i) - Psi_2_adjoint(i)*Grad_Psi_2(i))*dy;
end

val = val_1 + val_2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
