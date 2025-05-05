function UU_grand = decompose(UU_petit,Nbpt_grand)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% Attention, ne marche que lorsque (Nbpt_grand-1) divisible par (Nbpt_petit-1)
% Ex : Nbpt_grand = 101 , Nbpt_petit = 11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% declarations
% ------------
UU_grand = zeros(Nbpt_grand,1);
Nbpt_petit = length(UU_petit);

pas = (Nbpt_grand-1)/(Nbpt_petit-1);

for i=0:(Nbpt_petit-2)
    for j=1:pas
        UU_grand(i*pas+j) = UU_petit(i+1)*(1-((j-1)/pas)) + UU_petit(i+2)*((j-1)/pas);
    end
end
UU_grand(Nbpt_grand)=UU_petit(Nbpt_petit);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
