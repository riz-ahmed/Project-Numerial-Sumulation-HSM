% Lambda_e Berechnung für lumped Parameter Modell
function lae = lambda_e(a1,a2,a3,a4,a5,a6,a7)
lae = (a1.*a2./(a3.*a4.^2))*(a5.*a6./a7)^0.5;
end
