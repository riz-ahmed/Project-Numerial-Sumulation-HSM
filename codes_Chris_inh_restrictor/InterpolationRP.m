function [x,y,pt] = InterpolationRP(R,P,I,r0,ra,p0)

% Dieses Programm liest den normierten Druckverlauf ein, und führt vom
% Spaltausgang bis zur Stelle des Druckabfalls eine kubische Interpolation durch,
% um das extrapolierte p_t an der Stelle r0 zu ermitteln. 
% Hierzu wird die I-Funktion herangezogen. Sobald I=0.66 hat, handelt es
% sich um viskose Strömung. Wichtig dabei ist, dass dieses I erst nach
% dem lok Minimum der (R,I) Kurve gesucht wird.  

% Erzeuge einen Vektor, mit R und P, 
    PR = [R,P];

% Finde den viskosen Bereich
    % Stelle mit lok Minimum von I finden
        lok_min = find(I== min(I));
        
    % finde I mit Werten, die dem viskosen Bereich entsprechen
       viskos = find(I==0.66);
    % Bilde die Differenz der Stellen des lok Mins und den mit den Funktionswerten I=0.66   
       r_ast = viskos-lok_min(1);
    % Ist diese Differenz (wir interessieren uns für den Bereich hinter dem lok MIN)
    % dann lösche die Stelle aus r_ast
       li = r_ast<0;
       r_ast(li) =[];
    % Start des viskosen Bereiches ist dort, wo I=0.66 hat und hinter dem lok Min liegt   
       start_viskos = r_ast(1)+lok_min(1);

 % Entnormierung der Datei
    PR=[PR(start_viskos:end,1)*r0,PR(start_viskos:end,2)*p0];

% Kurven Fitting mit Polynom n.Grades
  p = polyfit(PR(:,1),PR(:,2),3);   % 3.Grad
  % p = polyfit(PR(:,1),PR(:,2),2); % 2.Grad
  % p = polyfit(PR(:,1),PR(:,2),1); % 1.Grad
  
  
% Erstellen des gesamten Funktionsverlaufs, um die Interpolierte Kurve zu plotten   
  i =1;
 
  for ki = r0:0.1e-3:ra
      y(i) = p(1)*ki^3+p(2)*ki^2+p(3)*ki+p(4);  % kubische Interpolation
      %y_exp(i) = p(1)*ki^2+p(2)*ki+p(3);       % quadratische Interpolation       
      %y_exp(i) = p(1)*ki+p(2);                 % lineare Interpolation
      x(i)= ki;
      i=i+1;
  end
 pt = y(1);
end