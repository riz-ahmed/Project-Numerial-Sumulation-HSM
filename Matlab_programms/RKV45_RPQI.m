function [R,P,Q,I,M0,re0stern] = RKV45_RPQI(ps,p0,T0,Rgas,kappa,my0,h0,h_r,r0,rp,ra,eps)    

% Mach's Number
M0       = ((2/(kappa-1))*((ps/p0)^((kappa-1)/kappa)-1))^0.5; %0.568;%
re0stern = (p0*M0*kappa^0.5*h0^2)/((Rgas*T0)^0.5*my0*r0); %57.6;%

% reading m and n values form a file
n_I = dlmread('I_n_2010_0bis1.txt');
m_I = dlmread('I_m_2010_0bis1.txt');
m_n_I = [m_I(:,1),n_I];

%% Range - Kutta Method

% Initial conditions
Q0 = 1;
I0 = 1;
P0 = 1;
R0 = 1;    
% simulation parameters 
Ra  = ra/r0;
N = 10000000;    % interval of R
s = (Ra-R0)/N; % step-size    
% applying the initial conditions     
R = R0;
P = P0;
Q = Q0;   
% finding the starting values of m and n
I_init = 1;
IDX = find(m_n_I(:,3) == I_init);
m = m_n_I(IDX,1);
n = m_n_I(IDX,2);
  
% Initializing the solution matrix
% nullmatrix = zeros(N+1,4);
% solution = [R0,Q0,P0,I0;nullmatrix];
solution = zeros(N+1,4);
solution(1,:) = [R0,P0,Q0,I0];
% ODE's 
dQ = @(R,Q,P,h_quer,re0stern,m)(1/(P*(h_quer^2)*re0stern))*m;
dP = @(R,Q,P,h_quer,re0stern,n,kappa,M0)(-kappa*M0^2*Q/(h_quer^2*re0stern))*n;
 
%% defining gradual step at the recess
Rp = rp /r0;
% with first and the last segment of r
r_seg = [r0 rp-0.3e-3 rp+0.3e-3 ra];
h_seg = [h0 h0 h_r h_r];
% creating cubic curve using pchip
r_chip = [rp-0.3e-3 :1e-6: ra];
h_chip = pchip(r_seg,h_seg,r_chip);
%% RK45 - with adaptive step-size
    for i = 1:N    
        R = R + s;
        % initializing gradual step w.r.t R
        if (R <= (rp-0.1e-3) / r0)
            H = 1;
        else
            h = pchip(r_chip, h_chip, R*r0);
            H = h/h0;
        end
        % exit condition
        if R>=Ra 
            break;
        end
        if R < Ra
            
            k1 = dQ(R,Q,P,H,re0stern,m);
            l1 = dP(R,Q,P,H,re0stern,n,kappa,M0);

            k2 = dQ(R+s/2,Q+s*k1/2,P+s*l1/2,H,re0stern,m);
            l2 = dP(R+s/2,Q+s*k1/2,P+s*l1/2,H,re0stern,n,kappa,M0);

            k3 = dQ(R+s/2,Q+s*k2/2,P+s*l2/2,H,re0stern,m);
            l3 = dP(R+s/2,Q+s*k2/2,P+s*l2/2,H,re0stern,n,kappa,M0);

            k4 = dQ(R+s/2,Q+s*k3/2,P+s*l3/2,H,re0stern,m);
            l4 = dP(R+s/2,Q+s*k3/2,P+s*l3/2,H,re0stern,n,kappa,M0);
            
            k5 = dQ(R+s/2,Q+s*k4/2,P+s*l4/2,H,re0stern,m);
            l5 = dP(R+s/2,Q+s*k4/2,P+s*l4/2,H,re0stern,n,kappa,M0);
            
            k = (k1 + 2*k2 + 2*k3 + k4)/6;
            l = (l1 + 2*l2 + 2*l3 + l4)/6;

            kf = (k1 + 2*k2 + 3*k3 + 2*k4 + k5)/9;
            lf = (l1 + 2*l2 + 3*l3 + 2*l4 + l5)/9;

            Q = Q + s*k;    
            Qf= Q + s*kf;
            
            P = P + s*l; 
            Pf= P + s*lf;
            
                      
            I = 1/(R*P*H*Q);
            If = 1/(R*Pf*H*Qf);
         
            %% adaptive step-size control
               if s < eps/(10*abs(lf-l)) && s < eps/(10*abs(kf-k))
                    s = 2*s;
               else
                if s > (Ra-R0)/N
                    s = s/2;
                end
               end                
            
            I = (round(1000*I))/1000;
            % break condition for I 
            if I > 1
               I_err = 1;
               break
            end
            if I < 0
               I_err = -1;
               break    
            end
            % finding m and n values for the next iteration
            IDX = find(m_n_I(:,3) == I);
            m = m_n_I(IDX,1);
            n = m_n_I(IDX,2);
            % updating solution matrix
            solution(i+1,:) = [R,P,Q,I];

       end

end
% sorting the resutls for the solution matrix
    R = solution(:,1);    
    P = solution(:,2);
    Q = solution(:,3);
    I = solution(:,4);    
    
% removing the zeros form the solution matrix, makes plotting more smooth  
    ri = R == 0;
    R(ri) = [];
    pi = P ==0;
    P(pi) = [];
    qi = Q ==0;
    Q(qi) = [];
    ii = I ==0;
    I(ii) = [];
    figure(1)

end