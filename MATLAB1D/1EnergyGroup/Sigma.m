function val = Sigma(x,epsilon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% INPUT - x : coordonnee du point ou on veut evaluer la fonction.
% OUTPUT - val: valeur de la matrice sur ce point.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%val = [[1, 0];[0, 1]];

%val = 10*(sin((2*pi*(x/epsilon)))+2);
%val = 15*((sin((pi*x)/epsilon)*sin((pi*x)/epsilon) + sin(sqrt(2)*(pi*x)/epsilon)*sin(sqrt(2)*(pi*x)/epsilon) )+1); %Quasiperiodic
%val = 14*(sin((pi*x)/epsilon)*sin((pi*x)/epsilon)+1); 
%val = (1+0.5*(x-0.5)*(x-0.5))*0.8*(10*cos((pi*x)/epsilon)*cos((pi*x)/epsilon)+1); %Cas pour W et Psi oscille + comportement macro
val = 1*(10*cos((pi*x)/epsilon)*cos((pi*x)/epsilon)+1); %Cas pour W et Psi oscille
%val = 10*(sin((2*pi*x)/epsilon)+2);
%val= 0;
%val = 20*(sin((2*pi*(x/epsilon)))+2);
%val = (0.01*sin(2*(pi*(x/epsilon)))+2);

%%% Diphasique %%%

% y = mod(x/epsilon, 1);
% 
% if y >= 0 && y < 0.5
%     val = 4;
% else
%     val = 1; % Optionally handle cases where x is outside the range [0, 1]
% end
%%%%%%%%%%%%%%%%%%%

% if((x/epsilon)>30)
%     val = 5.5;
% end

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