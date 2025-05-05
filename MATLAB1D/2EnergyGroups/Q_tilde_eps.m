function [vecteur_Q_tilde_eps] = Q_tilde_eps(h, epsilon, x_gauche, Phi_1,Phi_2,Phi_1_adjoint,Phi_2_adjoint,Num_coeff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% INPUT - x : coordonnee du point ou on veut evaluer la fonction.
% OUTPUT - val: valeur de la matrice sur ce point.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nbpt = size(Phi_1,1);

vecteur_Q_tilde_eps = zeros(Nbpt,1);

Q_12 = zeros(Nbpt,1);
Q_21 = zeros(Nbpt,1);
for i = 1:Nbpt
    Q_12(i) = Sigma(x_gauche+ (i-1)*h,epsilon,2)*(Phi_2(i)*Phi_1_adjoint(i)) ;
    Q_21(i) = Sigma(x_gauche+ (i-1)*h,epsilon,3)*(Phi_1(i)*Phi_2_adjoint(i)) ;
end
if Num_coeff==1 
    vecteur_Q_tilde_eps = -Q_12;
elseif Num_coeff ==2
    vecteur_Q_tilde_eps = Q_12;
elseif Num_coeff == 3
    vecteur_Q_tilde_eps = Q_21;
elseif Num_coeff == 4
    vecteur_Q_tilde_eps = -Q_21;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
