function KK = matK_per_Psi(Nbpt_per,Psi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matrice de raideur avec A remplacé par Psi**2*A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%h = 1/(Nbpt-1);
h = 1/(Nbpt_per);   %attention ici h est modifié car on précalcul une sous matrice
% declarations
% ------------

KK = sparse(Nbpt_per,Nbpt_per); % matrice de masse

% boucle pour les matrices EF
% ------------------------
for l=2:(Nbpt_per-1)   %La première et dernière ligne doivent rester nulle
    psi_avant = 0.5*(Psi(l)+Psi(l-1));
    psi_apres = 0.5*(Psi(l)+Psi(l+1));
    KK(l,l)=(1/h)*(psi_avant*psi_avant*A((l-1-0.5)*h,1) + psi_apres*psi_apres*A((l-1+0.5)*h,1)); %approximation de l'intégrale de a(.)
    KK(l,l+1)=(-1/h)*psi_apres*psi_apres*A((l-1+0.5)*h,1);
    KK(l+1,l)=(-1/h)*psi_apres*psi_apres*A((l-1+0.5)*h,1);
        
end % for l

% l=1
psi_avant = 0.5*(Psi(1)+Psi(Nbpt_per));
psi_apres = 0.5*(Psi(1)+Psi(2));
KK(1,1)=(1/h)*(psi_apres*psi_apres*A(0.5*h,1) + psi_avant*psi_avant*A(1-0.5*h,1));
KK(1,2)=(-1/h)*psi_apres*psi_apres*A(0.5*h,1);
KK(2,1)=(-1/h)*psi_apres*psi_apres*A(0.5*h,1);

% l=Nbpt_per
psi_avant = 0.5*(Psi(Nbpt_per)+Psi(Nbpt_per-1));
psi_apres = 0.5*(Psi(Nbpt_per)+Psi(1));
KK(Nbpt_per,Nbpt_per)=(1/h)*(psi_avant*psi_avant*A((Nbpt_per-1-0.5)*h,1) + psi_apres*psi_apres*A((Nbpt_per-1+0.5)*h,1));
KK(1,Nbpt_per)=(-1/h)*psi_apres*psi_apres*A(1-0.5*h,1);
KK(Nbpt_per,1)=(-1/h)*psi_apres*psi_apres*A(1-0.5*h,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
