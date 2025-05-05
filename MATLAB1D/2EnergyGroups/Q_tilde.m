function [vecteur_Q_tilde] = Q_tilde(Psi_1,Psi_2,Psi_1_adjoint,Psi_2_adjoint,Num_coeff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% INPUT - x : coordonnee du point ou on veut evaluer la fonction.
% OUTPUT - val: valeur de la matrice sur ce point.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nbpt = size(Psi_1,1);

vecteur_Q_tilde = zeros(Nbpt,1);
dy = 1/(Nbpt-1);
Q_12 = zeros(Nbpt,1);
Q_21 = zeros(Nbpt,1);
for i = 1:Nbpt
    Q_12(i) = Sigma((i-1)*dy,1,2)*(Psi_2(i)*Psi_1_adjoint(i)) ;
    Q_21(i) = Sigma((i-1)*dy,1,3)*(Psi_1(i)*Psi_2_adjoint(i)) ;
end
if Num_coeff==1 
    vecteur_Q_tilde = -Q_12;
elseif Num_coeff ==2
    vecteur_Q_tilde = Q_12;
elseif Num_coeff == 3
    vecteur_Q_tilde = Q_21;
elseif Num_coeff == 4
    vecteur_Q_tilde = -Q_21;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
