    % Main programm file for Shallow Pocket Feeding Static Bearing analysis
clear all; close all;clc; hold on; tic
%% Program constants

    ps      = 6e6;          % Static supply pressure in Pa
    p0      = 3.9e6;   % Starting choosen pressure in Pa
    T0      = 293;          % Standard room temprature in K
    Rgas    = 287.1;        % Universal gas constant for Air at STP
    kappa   = 1.4;          % Isentropic expansion co-efficient of air at STP
    my0     = 1.846e-5;     % Dynamic Viscosity of air at STP im Pa-s
    h       = 9e-6;        % orifice entrance height in �m 
    dp      = 20e-6;        % pocket height in �m
    h_r     = h;            % height as function of 'r'
    r0      = 100e-6;       % Orifice entrance radius in mm
    ra      = 6e-3;       % exit radius in mm
    rp      = 1e-3;       % Pocket radius in mm
    pa      = 101000;       % atmospheric pressure in Pa
    eps     = 1e-10;        % precision
    
%% Setting up the problem

% defining h0
h0 = h_r + dp;

[R,P,Q,I,M0,re0stern] = RKV45_RPQI(ps,p0,T0,Rgas,kappa,my0,h0,h_r,r0,rp,ra,eps);
pEnd = P(end)*p0;
plot(R*r0,P*p0,'b','LineWidth',1.5)

toc
