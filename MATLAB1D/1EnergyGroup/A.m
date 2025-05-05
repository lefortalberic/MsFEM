function val = A(x,epsilon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% INPUT - x : coordonnee du point ou on veut evaluer la fonction.
% OUTPUT - val: valeur de la matrice sur ce point.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%val = [[1, 0];[0, 1]];

%val = 1;
val = 1*(5*sin((pi*x)/epsilon)*sin(pi*(x/epsilon))+1);
%val = 5*(sin((pi*x)/epsilon) + sin(sqrt(2)*(pi*x)/epsilon))*(sin((pi*x)/epsilon) + sin(sqrt(2)*(pi*x)/epsilon))+2;
%val = 5*(sin((pi*x)/epsilon)*sin((pi*x)/epsilon) + sin(sqrt(2)*(pi*x)/epsilon)*sin(sqrt(2)*(pi*x)/epsilon) )+2; %Quasiperiodic
%val = 5*sin(2*(pi*x)/epsilon)+6;   %Cas pour W et Psi oscille
%val = 5*sin((pi*x)/epsilon+pi/4)*sin((pi*x)/epsilon+pi/4)+2;

%%% Diphasique %%%

% y = mod(x/epsilon, 1+eps);
% 
% if y >= 0 && y < 0.5
%     val = 1;
% elseif y >= 0.5 && y <= 1
%     val = 3;
% else
%     val = 0; % Optionally handle cases where x is outside the range [0, 1]
% end

%%%%%%%%%%%%%%%%%%

%val = 1;
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