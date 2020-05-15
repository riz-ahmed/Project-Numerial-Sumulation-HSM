% Hauptprogramm - Berechnung P(R),Q(R),I(R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hauptporgramm zur automatischen Berechnung des Druckverlaufs 
% mit der Möglichkeit auch mehrere h zu bestimmen

clc;
clear all;
close all
tic
% h eingeben!
%h=[6.03e-6 5.222e-6 10.44e-6 9.342e-6 7.627e-6 6.963e-6];
 h=15e-6;

hold on
%grid on
for zaehler1 = 1:1
    h0 = h(zaehler1);
    %ps = startp(zaehler1);

% Eingabe der restlichen Größen
    ps      = 600000;     % Speißedruck in Pa
    p00     = ps-1000;    % gewählter Startwert am Spalteingang in Pa
    T0      = 293;        % Speißetemperatur in K
    Rgas    = 287.1;      % Gaskonstante trockene Luft
    kappa   = 1.4;        % Isentropenexponent
    my0     = 17.8e-6;    % Viskosität Luft in Pa*s
    %h0      = 20e-6;      % Spalthöhe am Spalteingang in µm
    h_r     = h0;      % Funktion der Spalthöhe, abhängig von R
    r0      = 187.5e-6;
    ra      = 3e-3;
    pa      = 101500;
    eps     = 1e-10;
bert=0;
% Programminterne Werte 

    direction = 1;
    schritt = 50000;%(ps-p00)/4;
    pa_stern_neu=ps;
    precision = 250; % Toleranz, als Abbrucckriterium am Ende wenn abs(pa_stern_neu-pa)<precision
    
% 1.Schritt:
%   Finde ein p0 mit einem resultierendem pa_stern<=pa, um den Suchbereich
%   einzuschränken.

A=zeros(2,2);
z=1;
j=0;


while ((pa - pa_stern_neu) < 0)
    disp(['------------ Durchlauf ',num2str(j),' ------------']); 
    pa_stern_alt = ps;
    br=false;

    % Berechne pa_stern von sinkenden Werten von p0, und verfeinere die Suche im Bereich
    % der p0 um das Minimum herum. Nur dort kann der Bereich sein, wo pa
    % liegen kann
    for i = 0:9
        p0      = p00-i*schritt; 
        
        % Um die Berechnung zu beschleunigen wird jede Kombination aus p0 und daraus 
        % berechneten pa_stern im Array A gespeichert
        S=find(A(:,1)==p0);
        if (isempty(S))
            %[R,P,Q,I,I_err] = entrance_flowvol2(ps,p0,T0,Rgas,kappa,my0,h0,h_r,r0,ra,eps);
            [R,P,Q,I,I_err,M0,re0stern] = RKV45_RPQI(ps,p0,T0,Rgas,kappa,my0,h0,h_r,r0,ra,eps);
            pa_stern_neu = P(end)*p0;
            A(z,1)=p0;
            A(z,2)=pa_stern_neu;
            z=z+1;
        else
            pa_stern_neu = (A(S(1),2));
        end

        % Anzeigen der Ergebnisse
        % disp(['m_D: ',num2str(m0_punkt)]); 
%         disp(['p0: ',num2str(p0)]); 
%         disp(['pa_stern: ',num2str(pa_stern_neu)]); 
         plot(R*r0,P*p0,'Color','blue','LineWidth',1.5)
        
        % Um unnötige Berechnungen zu sparen wird die Berechnungsschleife
        % verlassen wenn das neu berechnete pa_stern wieder größer wird.
        % Denn dann hat man den Bereich in dem pa liegt, ausreichend
        % eingeschränkt
        if ((pa_stern_neu > pa_stern_alt) && (pa_stern_neu>pa))
            bert=1;
            break;
        end
        
        if pa_stern_neu<=pa
            bert=2;
            br=true;
        end
        
        if (abs(pa_stern_neu-pa)<100)
            bert=3;
            break;
        end
        pa_stern_alt = pa_stern_neu;
    end

    if (br==true)
        break;
    end
    % Bestimme neuer Suchbereich: p0_min-1<p0_min<p0_min+1
    A
    B=sortrows(A)
    [~,im]=min(A(:,2));
    
    % Bei kleinen Spalten in p00=A(im-1,1) "-1" löschen;
    if h0<5e-6
        p00=A(im,1);
    else
        p00=A(im-1,1);
    end
    schritt=abs(p00-p0)/10;
    

    j=j+1;
end

% 2.Schritt:
%   Bestimme p0-Bereich mit: pa_stern_1<pa<pa_stern_2

% Sortiere A gemäß der ersten Spalte aufsteigend um ersten Wert ab p0_min zu finden
% der größer pa ist
B=sortrows(A);
[c,i] = min(B);

% Finde erstes p0 größer als das p0_min mit pa_stern ist größer gleich pa
K=find(B(:,2)>=pa & B(:,1)>=B(i(2),1));

% Suchgrenzen bzgl p0
% Untere Grenze B(i,1)
u_p0=B(K(1)-1,1);
u_pa_stern=B(K(1)-1,2);
% Obere Grenze B(K(1),1)
o_p0=B(K(1),1);
o_pa_stern=B(K(1),2);

schritt = (o_p0-u_p0)/4;

% 3.Schritt:
%   Bestimme p0 mit abs(pa_stern-pa)<=precision
% z=1;
j=0;

while (abs(pa_stern_neu-pa)>precision)
    
    % Bestimme das p0 das mit pa_stern näher an pa dran ist, denn dort wird
    % pa wahrscheinlich in weniger Berechnungsschritten erreichbar sein
    if (abs(u_pa_stern-pa) > abs(o_pa_stern-pa))
        p00=o_p0;
        direction=-1;
    else
        p00=u_p0;
        direction=1;
    end
    
    disp(['------------ Durchlauf ',num2str(j),' ------------']); 
    
    for i=0:9    
        p0      = p00+i*schritt*direction; 
     
        % Um die Berechnung zu beschleunigen wird jede Kombination aus p0 und daraus 
        % berechneten pa_stern im Array A gespeichert
        S=find(A(:,1)==p0);
        if (isempty(S))
            [R,P,Q,I,I_err,M0,re0stern] = RKV45_RPQI(ps,p0,T0,Rgas,kappa,my0,h0,h_r,r0,ra,eps);
            pa_stern_neu = P(end)*p0;
            A(z,1)=p0;
            A(z,2)=pa_stern_neu;
            z=z+1;
        else
            %disp(['S: ',num2str(S(1))]);
            pa_stern_neu = (A(S(1),2));
        end

        

        % Anzeigen der Ergebnisse
        % disp(['m_D: ',num2str(m0_punkt)]); 
        disp(['p0: ',num2str(p0)]); 
        disp(['pa_stern: ',num2str(pa_stern_neu)]); 
        plot(R*r0,P*p0,'Color','green')
        
        % Abbruchbedingung abhängig davon ob man von der oberen Grenze
        % runter zählt oder von der Unteren nach oben.
        %   "direction" bestimmt ob der Wert in Schritt addiert oder
        %       subtrahiert wird um neues p0 zu bestimmen
        if direction==-1
            if pa_stern_neu < pa
                break;
            end
        else
            if pa_stern_neu > pa
                break;
            end
        end
        
        % Abbruchbedingung: Gesuchter Wert gefunden
        if (abs(pa_stern_neu-pa)<precision)
            break;
        end
    end

    % Hier wird wie oben der Bereich um pa bestimmt, pa wird also immer
    % weiter eingegrenzt
    B=sortrows(A);
    [c,i] = min(B);

    % Finde erstes p0 mit pa_stern ist größer gleich pa
    K=find(B(:,2)>=pa & B(:,1)>=B(i(2),1));

    % Suchgrenzen bzgl p0
    % Untere Grenze B(i,1)
    u_p0=B(K(1)-1,1);
    u_pa_stern=B(K(1)-1,2);
    % Obere Grenze B(K(1),1)
    o_p0=B(K(1),1);
    o_pa_stern=B(K(1),2);
    
    schritt = (o_p0-u_p0)/4;
    j=j+1;
    
    if schritt<0.0001
        break;
    end
end       
  
%clc

%[m0_punkt] = m0punkt(p0,ps,kappa,r0,h0,Rgas,T0);
B
% disp('------------ Ergebnisse ------------'); 
% % Anzeigen der Ergebnisse
[K,l]=min(abs(B(:,2)-pa));
p0=B(l,1);
pa_stern=B(l,2);
%disp(['m_D: ',num2str(m0_punkt)]); 
disp(['p0: ',num2str(p0)]);
disp(['pa_stern: ',num2str(P(end)*p0)]); 

[R,P,Q,I,I_err,M0,re0stern] = RKV45_RPQI(ps,p0,T0,Rgas,kappa,my0,h0,h_r,r0,ra,eps);
plot(R*r0,P*p0,'Color','red')
grid on
ylabel('Pressure in Pa')
xlabel('r in m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Berechnung der Interpolation und der extrapolierten Druckkurve des
% viskosen Bereichs. Mit der Ausgabe von pt

[x,y,pt] = InterpolationRP(R,P,I,r0,ra,p0);
 pt_quer = pt/ps; 
 hold on
 plot(x,y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Berechnung des Massenstroms nach Waumanns, function enthält noch weitere 
% Formeln. Die Formel von Waumanns ist validiert, Tang68 passt
% einigermaßen, die klassische viskose Formel konnte nicht validiert werden
 
[m0_punkt] = Massenstrom(p0,ps,r0,h0,Rgas,T0,kappa);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Berechnung des Discharge Koeffizienten Cd

[Cd] = discharge_koeffizient(p0,ps,pt,kappa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Groß Lambda Berechnung für Pertubation-Method
 Grosslambda = lambda_e(my0,r0,ps,h0,Rgas,T0,kappa);
 
 toc
% optionales Schreiben einer Datei
werte = [h0,r0,ra,ps,p0,pt,Cd,Grosslambda];

dlmwrite('werte.txt',[werte],'delimiter',' ','newline','pc','-append');

end