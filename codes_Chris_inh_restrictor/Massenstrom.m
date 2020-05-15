% Berechnung von m0_punkt mittels p0 und ps 
% für die Formel von Tang zu verwenden muss M0 eingelesen werden!

function [m0_punkt] = Massenstrom(p0,ps,r0,h0,Rgas,T0,kappa)%,M0,ra)
pi =3.1416;
A0 = 2*pi*r0*h0;
p0_quer = p0/ps;
pc_quer = (2/(kappa+1))^(kappa/(kappa-1));

if p0_quer>=pc_quer
    phi_e = (p0_quer^(2/kappa)-p0_quer^((kappa+1)/kappa))^.5;
else
    phi_e = (pc_quer^(2/kappa)-pc_quer^((kappa+1)/kappa))^.5;
end
    
m0_punkt = A0*((2*kappa)/(kappa-1))^0.5*ps/((Rgas*T0)^0.5)*phi_e;   %Waumanns
%m0_punkt = M0*p0*2*3.1416*h0*r0*(kappa/(Rgas*T0))^0.5;             %Tang
%m0_punkt = pi*h0^3*(pt^2-pa^2)/(12*my0*Rgas*T0*log(ra/r0));        %Radial Flow
end