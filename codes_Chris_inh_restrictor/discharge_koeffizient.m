% Berechnung des Discharge Koeffizienten Cd
function [Cd] = discharge_koeffizient(p0,ps,pt,kappa)

p0_quer = p0/ps;
pt_quer = pt/ps;
pc_quer = (2/(kappa+1))^(kappa/(kappa-1)); %kritisches Druckverhältnis

    if p0_quer>=pc_quer
        phi_e = (p0_quer^(2/kappa)-p0_quer^((kappa+1)/kappa))^.5;
    else
        phi_e = (pc_quer^(2/kappa)-pc_quer^((kappa+1)/kappa))^.5;
    end

    if pt_quer>=pc_quer
        phi_e_pt = (pt_quer^(2/kappa)-pt_quer^((kappa+1)/kappa))^0.5;
    else
        phi_e_pt = (pc_quer^(2/kappa)-pc_quer^((kappa+1)/kappa))^0.5;
    end

Cd = phi_e/phi_e_pt;    
  
end