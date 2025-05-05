function val = Phi(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% INPUT - x : coordonnee du point ou on veut evaluer la fonction.
% OUTPUT - val: valeur de la matrice sur ce point.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if((x > 0) && (x < 1))
%     val = exp(-1/(1-(2*x-1)*(2*x-1)));
% else 
%     val = 0;
% end
% norme_Phi_carre = 0.412592943718215;
%val = val + 1;
%val = (val)/sqrt(norme_Phi_carre);
%val = 630*(1-x)*(1-x)*(1-x)*(1-x)*x*x*x*x; % Pour que la norme L^1 soit 1 
val = 30*(1-x)*(1-x)*x*x; % Pour que la norme L^1 soit 1 
%val = 6*(1-x)*x;
%val = 1.61805*sqrt(sqrt(x*(1-x)));
%val = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %Pour plot :
% Nbpt_h = 1001;
% h = 1/(Nbpt_h -1)
% 
% %Maillage
% PhiXX = zeros(Nbpt_h,1);     %vecteur des coordonnées des noeuds
% for j=1:Nbpt_h
%     PhiXX(j)=Phi((j-1)*h);
% end