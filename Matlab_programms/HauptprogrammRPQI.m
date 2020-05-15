% Main programm for simulation of shallow bearing P(R), I(R),Q(R) characteristics
% includes function file RKM45_RPQI for solving IVP with Range-kutta 4th order method 
clc;
clear all;
close all
tic
hold on

%% Programm constants
    ps      = 500000;     % Static pressure in Pa
    p00     = ps-1000;    % Initial choosen pressure in Pa
    T0      = 293;        % Supply temperature in K
    Rgas    = 287.1;      % Universal Gas constant of air in KJ/Kg-k
    kappa   = 1.4;        % ratio of specific heats for air at T0
    my0     = 17.8e-6;    % Dynamic viscosity of air at T0
    dp      = 20e-6;      % Feed pocket depth in m
    h_r     = 11e-6;      % second restrictor height in m
    h0      = dp + h_r;   % entrance gap height in m
    r0      = 100e-6;     % entrance radius in m
    rp      = 2e-3;       % radius at recess in m
    ra      = 7.5e-3;     % radius at bearing exit in m
    pa      = 101500;     % atmospheric pressure in Pa
    eps     = 1e-10;      % precsion control for step size

%% Simulation of the pressure characteristics and displaying intial choosen results
% Simulation Parameters
steps = 50000;      % step's for choosing p0
pEnd_new=ps;        % pEnd is the bearing pressure at the end
    
% initializing solution matrix
A=zeros(2,2);       % initialized solution matrix for storing and sorting results
idx=1;              % row index
k=1;                % loop counter

while ((pa - pEnd_new) < 0)
    disp(['========== Trial: ',num2str(k),' ==========']); 
    pEnd_old = ps;

    % simlation with progressive step-size and p0 reduction for initial
    % solution determination
    for i = 0:9
        p0      = p00-i*steps;         
        % checking p0 from A, if not found execute RKM45
        S=find(A(:,1)==p0);
        if (isempty(S))            
            [R,P,Q,I,M0,re0stern] = RKV45_RPQI(ps,p0,T0,Rgas,kappa,my0,h0,h_r,r0,rp,ra,eps);
            pEnd_new = P(end)*p0;
            A(idx,1)=p0;
            A(idx,2)=pEnd_new;
            idx=idx+1;
        else
            % if p0 already exits, pEnd matches the corresponding value
            pEnd_new = (A(S(1),2));
        end

        % plotting the results w.r.t r
        plot(R*r0,P*p0,'Color','b','LineWidth',1.5)
        
        % break condition's for the loop
        if ((pEnd_new > pEnd_old) && (pEnd_new>pa))
            break
        end        
        if pEnd_new<=pa
            break
        end        
        if (abs(pEnd_new-pa)<100)
            bert=3;
            break;
        end
        pEnd_old = pEnd_new;
    end

    if pEnd_new<=pa
        break
    end
    % results and choosing step-size
    A
    B=sortrows(A)
    [~,im]=min(A(:,2));
    % choosing p00
    p00=A(im-1,1);
    steps=abs(p00-p0)/10;  

    k=k+1;
end

%% Second phase of the simulation - for simulating with smaller step-size of controlled results from above
B=sortrows(A);
[~,i] = min(B);

% narrowing results form first simulation based on pEnd > pa and p0
% corresponding to min p0
% K defines the logical array of the short-listed matrix indeces
K=find(B(:,2)>=pa & B(:,1)>=B(i(2),1));

% choosing the lower limit
l_p0=B(K(1)-1,1);           % limit is taken lower than the least pEnd but pEnd > pa
l_pEnd=B(K(1)-1,2);
% choosing the upper limit
u_p0=B(K(1),1);
u_pEnd=B(K(1),2);

steps = (u_p0-l_p0)/4;      % p0 is choosen in steps
k=1;                        % loop counter
direction = 1;              % simulation direction
precision = 250;            % set for limiting the choosing upper and lower limit
while (abs(pEnd_new-pa)>precision)
    
    if (abs(l_pEnd-pa) > abs(u_pEnd-pa))
        p00=u_p0;
        direction=-1;
    else
        p00=l_p0;
        direction=1;
    end
    
    disp(['========== Trial ',num2str(k),' ===========']); 
    
    for i=0:9    
        p0      = p00+i*steps*direction; 
     
        % checking p0 from A, if not found execute RKM45
        S=find(A(:,1)==p0);
        if (isempty(S))
            [R,P,Q,I,M0,re0stern] = RKV45_RPQI(ps,p0,T0,Rgas,kappa,my0,h0,h_r,r0,rp,ra,eps);
            pEnd_new = P(end)*p0;
            A(idx,1)=p0;
            A(idx,2)=pEnd_new;
            idx=idx+1;
        else
            pEnd_new = (A(S(1),2));
        end
        % displaying the results
        disp(['p0: ',num2str(p0)]); 
        disp(['pBearingEnd: ',num2str(pEnd_new)]); 
        plot(R*r0,P*p0,'Color','green')
        
        if direction==-1
            if pEnd_new < pa
                break;
            end
        else
            if pEnd_new > pa
                break;
            end
        end
        
        % loop break condition upon reaching the limit defined by precision
        if (abs(pEnd_new-pa)<precision)
            break;
        end
    end

    % choosing the newer limits for further refining the choosen limits of
    % solutions
    B=sortrows(A);
    [~,i] = min(B);

    % choosing the required simulation values based on pEnd < pa and p0
    % corresponding the smallest pEnd
    K=find(B(:,2)>=pa & B(:,1)>=B(i(2),1));

    % choosing lower limit
    l_p0=B(K(1)-1,1);
    l_pEnd=B(K(1)-1,2);
    % choosing upper limit
    u_p0=B(K(1),1);
    u_pEnd=B(K(1),2);
    % defining new step size
    steps = (u_p0-l_p0)/4;
    k=k+1;
    % defining the break condition
    if steps<0.0001
        break;
    end
end       
  
%% Final simulation phase for displaying final short-listed solution
B                               % display's the final sorted solution matrix
[K,l]=min(abs(B(:,2)-pa));
p0=B(l,1);
pa_stern=B(l,2);
% displaying the final p0 and pEnd at the end of the simulation 
disp(['p0: ',num2str(p0)]);
disp(['pBearingExit: ',num2str(P(end)*p0)]); 

[R,P,Q,I,M0,re0stern] = RKV45_RPQI(ps,p0,T0,Rgas,kappa,my0,h0,h_r,r0,rp,ra,eps);
plot(R*r0,P*p0,'Color','red')
grid on
ylabel('Pressure in Pa')
xlabel('r in m')

toc

% end of the simulation
