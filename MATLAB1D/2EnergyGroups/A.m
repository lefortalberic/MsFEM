function val = A(x,epsilon, Num_coeff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% INPUT - x : coordonnee du point ou on veut evaluer la fonction.
%Num_coeff = 1 ou 2. Num_coeff = 1 Pour calculer la matrice de raideur
% correspondant au premier niveaux d'energie
% OUTPUT - val: valeur de la matrice sur ce point.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%val = [[1, 0];[0, 1]];

%val = 1;
%val = 1*(5*sin((pi*x)/epsilon)*sin(pi*(x/epsilon))+1);
%val = 5*(sin((pi*x)/epsilon) + sin(sqrt(2)*(pi*x)/epsilon))*(sin((pi*x)/epsilon) + sin(sqrt(2)*(pi*x)/epsilon))+2;
% val = 5*(sin((pi*x)/epsilon)*sin((pi*x)/epsilon) + sin(sqrt(2)*(pi*x)/epsilon)*sin(sqrt(2)*(pi*x)/epsilon) )+2; %Quasiperiodic
%val = 5*sin(2*(pi*x)/epsilon)+6;   %Cas pour W et Psi oscille
%val = 5*sin((pi*x)/epsilon+pi/4)*sin((pi*x)/epsilon+pi/4)+2;
%val = [5*sin((pi*x)/epsilon)*sin(pi*(x/epsilon))+1 , 5*cos((pi*x)/epsilon)*cos((pi*x)/epsilon)+6 ];


val = [5*sin((pi*x)/epsilon)*sin(pi*(x/epsilon))+3.5 , 5*cos((pi*x)/epsilon)*cos(pi*(x/epsilon))+3.5 ]; %Cas pour condition de symétrie

% val = [1*(sin((pi*x)/epsilon)*sin((pi*x)/epsilon) + sin(sqrt(3)*(pi*x)/epsilon)*sin(sqrt(3)*(pi*x)/epsilon) )+8 ...
%     , 1*(cos((pi*x)/epsilon)*cos((pi*x)/epsilon) + cos(sqrt(2)*(pi*x)/epsilon)*cos(sqrt(2)*(pi*x)/epsilon) )+8 ]; %Cas pour QP

%%%%%%%%%%%%%%%%%%
val = val(Num_coeff);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Pour plot A(X,epsilon) :
% Nbpt_h = 1001;
% h = 1/(Nbpt_h -1)
% epsilon = 0.1;
% %Maillage
% AXX = zeros(Nbpt_h,1);     %vecteur des coordonnées des noeuds
% for j=1:Nbpt_h
%     AXX(j)=A((j-1)*h,epsilon);
% end