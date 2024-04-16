%MSE 855
%Homework 9
%4/16/2024
clear
clc
close all


T_pc = 911+273; %temp of phase change FCC to BCC (K)
T= 600:1:T_pc+25; %(K)
t =  [0.001:0.001:1]*10^0; %time in seconds

%data from Barin (converted molar G to volumetric G)
T_alpha = [298.15 300.00 400.00 500.00 600.00 700.00 800.00 900.00 1000.00 1042.00 1042.00 1100.00 1184.00]';
T_beta = [1184.00 1200.00 1300.00 1400.00 1500.00 1600.00 1665.00]';
G_alpha = ([-8.133 -8.184 -11.317 -15.140 -19.559 -24.513 -29.964 -35.891 -42.302 -45.153 -45.153 -49.252 -55.437])'; % kJ/mol
G_beta = ([-55.437 -56.656 -64.437 -72.486 -80.790 -89.338 -95.020])'; %kJ/mol

%create 2nd order polynomial fits
fit_alpha = polyfit(T_alpha, G_alpha, 3);
fit_beta = polyfit(T_beta, G_beta, 2);

%evaluate the fits over the temperature range desired
alpha_eval = polyval(fit_alpha, T-24); %shifted temp range slightly to set 0 at phase change of 911C
beta_eval = polyval(fit_beta, T-24); %shifted temp range slightly to set 0 at phase change of 911C

%finally calculate volumetric delta G
delta_g_v = (beta_eval-alpha_eval).*1000*1000000*7.8/55.835; %J/m^3


n_o = 10^20; %n zero (volumetric prefactor for n*) 1/m^3
gamma = 0.04; %J/m^2
k_b = 1.3806452*10^-23; %boltzman constant in J/K
a_fcc = 0.3571*10^-9; %meters


G_dot = 4.898*10^12 .* exp(-2.943./(T.*8.6173303*10^(-5))); %m/s

D = 0.49*10^-4.*exp(-284000./(8.314.*T)); %m^2/s
nu = 6.*D./(a_fcc^2); %jump frequency (1/s)

n_star = n_o.*exp(-16*pi*gamma^3./(3.*delta_g_v.^2.*T.*k_b)); %1/m^3
N = nu.*n_star; %volumetric nucleation rate (1/s*m^3)

Transf_rate = G_dot.*N;

figure
scatter(T_alpha-273, G_alpha,'blue')
hold on
scatter(T_beta-273, G_beta,'red')
plot(T-273, alpha_eval,'blue')
plot(T-273, beta_eval,'red')
xlabel ("Temperature (C)")
ylabel("Gibbs Free Energy (kJ/mol)")
legend("\alpha","\gamma",'fontsize',18)

figure
plot(T-273, delta_g_v)
hold on
yscale("log")
ylabel("delta G_v (J/m^3)")
xlabel ("Temp (C)")
legend("Curve Fit Estimation")



figure
plot(N./max(N), T-273, 'blue')
hold on
plot(G_dot./max(G_dot),T, 'red')
plot(Transf_rate./max(Transf_rate), T,'green')
xlabel("Normalized Rates")
ylabel("Temperature (C)")
legend("Nucleation Rate", "Growth Rate", "Normalized Transformation Rate")



for ii = 1:length(t)
    for jj = 1:length(T)
        f_beta(ii,jj) = 1-exp((-1/3).*pi.*G_dot(jj).^3.*N(jj).*t(ii).^4);
    end
end

figure

contour(t,T-273,f_beta', [0.01, .25, .5, .75, .99] ,'ShowText',true)
clim([0,1.1])
yline(911,'LineStyle','--','Color','r')
xlabel("time (s)",'FontSize',14)
ylabel("Temperature (C)",'FontSize',14)
xscale("log")
title("TTT Curve of \gamma-Fe to \alpha-Fe",'FontSize',16)



