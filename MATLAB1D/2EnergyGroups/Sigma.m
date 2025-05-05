function val = Sigma(x,epsilon, Num_coeff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% INPUT - x : coordonnee du point ou on veut evaluer la fonction.
% Num_coeff = 1 ou 2 ou 3 ou 4. Num_coeff = 1 Pour calculer la
% matrice de M_sigma
% correspondant a Sigma_11
% OUTPUT - val: valeur de la matrice sur ce point.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%val = [[1, 0];[0, 1]];

%val = 15*((sin((pi*x)/epsilon)*sin((pi*x)/epsilon) + sin(sqrt(2)*(pi*x)/epsilon)*sin(sqrt(2)*(pi*x)/epsilon) )+1); %Quasiperiodic


%val = 10*(sin((2*pi*x)/epsilon)+2);
%val = 20*(sin((2*pi*(x/epsilon)))+2);
% val = [10*cos((pi*x)/epsilon)*cos((pi*x)/epsilon)+1 , -10*sin((pi*x)/epsilon)*sin((pi*x)/epsilon)+1 , ...
%     -10*sin((pi*x)/epsilon)*sin((pi*x)/epsilon)+1 , 10*cos((pi*x)/epsilon)*cos((pi*x)/epsilon)+1]; %Cas pour W et Psi oscille


%Cas periodique pour condition de symétrie
val = [25*cos((pi*x)/epsilon)*cos((pi*x)/epsilon)+1 , -2*sin((pi*x)/epsilon)*sin((pi*x)/epsilon)-1 , ...
    -5*cos((pi*x)/epsilon)*cos((pi*x)/epsilon)-1 , 12*cos((pi*x)/epsilon)*cos((pi*x)/epsilon)+1]; 

% %Quasiperiodique symétrique sigma12 quelconque
% val = [15*((sin((pi*x)/epsilon)*sin((pi*x)/epsilon) + sin(sqrt(2)*(pi*x)/epsilon)*sin(sqrt(2)*(pi*x)/epsilon) )+10) , -15*sin((pi*x)/epsilon)*sin((pi*x)/epsilon)+10 , ...
%     -5*cos((pi*x)/epsilon)*cos((pi*x)/epsilon)+10 , 15*((sin((pi*x)/epsilon)*sin((pi*x)/epsilon) + sin(sqrt(3)*(pi*x)/epsilon)*sin(sqrt(3)*(pi*x)/epsilon) )+10)]; 

% %Quasiperiodique symétrique sigma12 negatif
% val = [20*((sin((pi*x)/epsilon)*sin((pi*x)/epsilon) + sin(sqrt(2)*(pi*x)/epsilon)*sin(sqrt(2)*(pi*x)/epsilon) )+10) , -15*sin((pi*x)/epsilon)*sin((pi*x)/epsilon)-1 , ...
%     -5*cos((pi*x)/epsilon)*cos((pi*x)/epsilon)-1 , 20*((sin((pi*x)/epsilon)*sin((pi*x)/epsilon) + sin(sqrt(3)*(pi*x)/epsilon)*sin(sqrt(3)*(pi*x)/epsilon) )+10)]; 

% %Cas periodique pour condition de non symétrie
% val = [10*cos((pi*x)/epsilon)*cos((pi*x)/epsilon)+1 , -10*sin((pi*x)/epsilon)*sin((pi*x)/epsilon)-1 , ...
%     -5*cos((pi*x)/epsilon + pi/4)*cos((pi*x)/epsilon+ pi/4)-1 , ...
%     10*cos((pi*x)/epsilon)*cos((pi*x)/epsilon)+1]; 

% %Quasiperiodique non symétrique sigma12 quelconque
% val = [15*((sin((pi*x)/epsilon)*sin((pi*x)/epsilon) + sin(sqrt(2)*(pi*x)/epsilon)*sin(sqrt(2)*(pi*x)/epsilon) )+10) , -15*sin((pi*x)/epsilon+ pi/4)*sin((pi*x)/epsilon+ pi/4)+10 , ...
%     -5*cos((pi*x)/epsilon)*cos((pi*x)/epsilon)+10 , 15*((sin((pi*x)/epsilon)*sin((pi*x)/epsilon) + sin(sqrt(3)*(pi*x)/epsilon)*sin(sqrt(3)*(pi*x)/epsilon) )+10)]; 

% %Quasiperiodique non symétrique sigma12 negatif
% val = [15*((sin((pi*x)/epsilon)*sin((pi*x)/epsilon) + sin(sqrt(2)*(pi*x)/epsilon)*sin(sqrt(2)*(pi*x)/epsilon) )+10) , -15*sin((pi*x)/epsilon+ pi/4)*sin((pi*x)/epsilon+ pi/4)-1 , ...
%     -5*cos((pi*x)/epsilon)*cos((pi*x)/epsilon)-1 , 15*((sin((pi*x)/epsilon)*sin((pi*x)/epsilon) + sin(sqrt(3)*(pi*x)/epsilon)*sin(sqrt(3)*(pi*x)/epsilon) )+10)]; 

val = val(Num_coeff);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% %Pour plot sigma(X,epsilon) :
% Nbpt_h = 1001;
% h = 1/(Nbpt_h -1)
% epsilon = 0.01;
% %Maillage
% AXX = zeros(Nbpt_h,1);     %vecteur des coordonnées des noeuds
% for j=1:Nbpt_h
%     AXX(j)=Sigma((j-1)*h,epsilon);
% end