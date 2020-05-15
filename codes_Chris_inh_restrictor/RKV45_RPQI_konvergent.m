function [R,P,Q,I,I_err,M0,re0stern] = RKV45_RPQI_konvergent(ps,p0,T0,Rgas,kappa,my0,h0,nue,r0,ra,eps)    

% Funktion der Spalthöhe, abhängig von R und einem Winkel "nue"

    h_r = @(h0,nue,R,r0) h0 - r0*(R-1)*tan(nue);
    
% Konstanten

    M0       = ((2/(kappa-1))*((ps/p0)^((kappa-1)/kappa)-1))^0.5; %0.568;%
    re0stern = (p0*M0*kappa^0.5*h0^2)/((Rgas*T0)^0.5*my0*r0); %57.6;%
    I_err = 0;

% Einlesen der n und m Dateien:
    n_I = dlmread('I_n_2010_0bis1.txt');
    m_I = dlmread('I_m_2010_0bis1.txt');
    m_n_I = [m_I(:,1),n_I];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runge Kutta Verfahren  

% Randbedingungen
    Q0 = 1;
    I0 = 1;
    P0 = 1;
    R0 = 1;
    
% Schrittweite  
    Ra  = ra/r0;
    intervalle = 1000000;      % Anzahl der Subintervalle
    j = 1;         % soll nur jeden j-ten Schritt später in Plot darstellen
    s = (Ra-R0)/(intervalle); % Schrittweite
    
% Setzen der Anfangswerte    
    R = R0;
    P = P0;
    Q = Q0;
    %I = I0;    
    
% finde Startwerte für m und n
    I_such = 1;
    Zeile = find(m_n_I(:,3) == I_such);
    m = m_n_I(Zeile,1);
    n = m_n_I(Zeile,2);
  
% Ergebnismatrix wird geschrieben
    nullmatrix = zeros(n,4);
    ergebnis = [R0,Q0,P0,I0;nullmatrix];
% DGLs 

    dQ = @(R,Q,P,h_quer,re0stern,m)(1/(P*(h_quer^2)*re0stern))*m;
    dP = @(R,Q,P,h_quer,re0stern,n,kappa,M0)(-kappa*M0^2*Q/(h_quer^2*re0stern))*n;
    % I = 1/(R*h_quer*P*Q);

for i = 1:intervalle
  
                R = R + s;   

        h_quer   = h_r(h0,nue,R,r0)/h0;
        if R>=Ra
            break;
        end;
        if R < Ra
            
            k1 = dQ(R,Q,P,h_quer,re0stern,m);
            l1 = dP(R,Q,P,h_quer,re0stern,n,kappa,M0);


            k2 = dQ(R+s/2,Q+s*k1/2,P+s*l1/2,h_quer,re0stern,m);
            l2 = dP(R+s/2,Q+s*k1/2,P+s*l1/2,h_quer,re0stern,n,kappa,M0); % 1. Mittelpunktssteigung


            k3 = dQ(R+s/2,Q+s*k2/2,P+s*l2/2,h_quer,re0stern,m);
            l3 = dP(R+s/2,Q+s*k2/2,P+s*l2/2,h_quer,re0stern,n,kappa,M0); % 2. Mittelpunktssteigung


            k4 = dQ(R+s/2,Q+s*k3/2,P+s*l3/2,h_quer,re0stern,m);
            l4 = dP(R+s/2,Q+s*k3/2,P+s*l3/2,h_quer,re0stern,n,kappa,M0); % rechte Steigung

            k5 = dQ(R+s/2,Q+s*k4/2,P+s*l4/2,h_quer,re0stern,m);
            l5 = dP(R+s/2,Q+s*k4/2,P+s*l4/2,h_quer,re0stern,n,kappa,M0);

            k = (k1 + 2*k2 + 2*k3 + k4)/6; % gemittelter k-Wert
            l = (l1 + 2*l2 + 2*l3 + l4)/6; % gemittelter l-Wert

            kf = (k1 + 2*k2 + 3*k3 + 2*k4 + k5)/9; % gemittelter k-Wert
            lf = (l1 + 2*l2 + 3*l3 + 2*l4 + l5)/9; % gemittelter l-Wert


            Q = Q + s*k;            % Eulerverfahrensschritt für neues Q   
            Qf= Q + s*kf;
            
            P = P + s*l;            % Eulerverfahrensschritt für neues P
            Pf= P + s*lf;
            
            I = 1/(R*P*Q);
            If = 1/(R*Pf*Qf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        % Schrittweitensteuerung
            
            % Bedingung Betragsdifferenz von Qf - Q, bzw.
            % Pf - P muss kleiner sein als neue
            % Schrittweite mal Epsilon
     
   
               if s < eps/(10*abs(lf-l)) && s < eps/(10*abs(kf-k))
                    s = 2*s;
               else
                if s > (Ra-R0)/(intervalle);
                    s = s/2;
                end
               end
                
            
             I = (round(1000*I))/1000;
%             
%               if R>r1/r0 && sprung == 1
%                      I=1;
% %                    M0       = ((2/(kappa-1))*((ps/(P*p0))^((kappa-1)/kappa)-1))^0.5; %0.568;%
% %                    re0stern = (P*p0*M0*kappa^0.5*h1^2)/((Rgas*T0)^0.5*my0*r1); 
%                      sprung = 0;
% %                end
%               end
%  
            
            
          % Abbruchbedingungen, wenn I > 1 wird oder <0,46 werden in den m und n Dateien keine Daten mehr gefunden!   
            
            if I > 1
               I_err = 1;
               break
            end
            if I < 0
               I_err = -1;
               break    
            end

            Zeile = find(m_n_I(:,3) == I);
            m = m_n_I(Zeile,1);
            n = m_n_I(Zeile,2);

    
          if floor(i/j) == i/j      % Schrittweite im Eregbnisvektor ganzzahlige Vielfache von j
            ergebnis(i+1,:) = [R,P,Q,I];
          end

       end

end
% Da der Vektor "ergebnis" aus Geschwindigkeitsgründen mit Nullen besetzt wird,
% wurde ein Filter eingebaut um Nullen herauszufiltern. Dies passiert, wenn
% j > 1 ist.
    R = ergebnis(:,1);    
    P = ergebnis(:,2);
    Q = ergebnis(:,3);
    I = ergebnis(:,4);    

    ri = R == 0;
    R(ri) = [];
    pi = P ==0;
    P(pi) = [];
    qi = Q ==0;
    Q(qi) = [];
    ii = I ==0;
    I(ii) = [];

end
