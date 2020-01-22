%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  code  : Inversion_OSLsurface_exposure_age_resampling.m
% Version     : 20.01.2020
% Description : Inversion of unknown exposure age using sigmaphi and mu
% Specificity : division with standard deviation on the plateau of the experimental data
% Resampling  : Resampling of the likelihood by comparing to a random value
%               between 0 and 1
% Units       : SP  [a-1]
%             : mu  [mm-1]
%             : Age [a]
%             : Ddot_input [Gy ka-1] conversion to [Gy a-1] in line 53
%             : D0  [Gy]
%             : Magic number = Ddot/D0 [a-1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
clc;
close all;

TT  = 10000;

% Loading data file excel
SampleName     = 'OSL_MBTP7_cal';
[num]          = xlsread(SampleName);
n              = length(num(:,1));

% Input parameters from calibration
SP0 = 2.86e-06;                      % [s-1]
mu  = 0.565 ;                        % [mm-1]
SP  = SP0*365.25*24.*3600;           % [a-1]

ind_plt      = 40;                   % at which index the plateaus starts considering all cores as one single dataset

%% Parameterization
% Entrer de la variable temps

tmin = 0;
tmax = 20;

% Definitions de donnees experimentales

x(1:n)          = (num(:,1));
L(1:n)          = (num(:,2));
e(1:n)          = (num(:,3));
[x_s,ind]       = sort(x(:));
Ls_M              = L(ind);
a               = std(Ls_M(ind_plt:n));

Ddot_input      = num(1,5);
Ddot            = Ddot_input/(1e3);       % [Gy ka-1] ==> [Gy a-1]
D0              = 500;                    % [Gy]  
magicN          = Ddot/D0;                % [a-1]

%% Compute residuals (i.e. fit to data)

M         = nan(TT,1);
t_vec     = nan(TT,1);

h         = waitbar(0,'Put a good song and relax');  % creation of the loading bar (just for fun)
rand_vec  = rand(TT,1);

r_t1      = sort(tmin+(tmax-tmin)*rand_vec);

for i = 1:TT
            
            M(i)   = 0;           
            L_th   = exp(-SP*r_t1(i)*exp(-mu*x));        % Equation without taking the dose rate in account
%            L_th   = (SP.*exp(-mu.*x).*exp(-r_t1(i).*(SP.*exp(-mu.*x)+magicN))+magicN)./(SP.*exp(-mu.*x)+magicN);    % Sohbati et al. 2012a
            M(i)   = nansum((L-L_th).^2/a^2);             % L2 norm weighted with the noise a to calculate the misfit M         
            t_vec(i)  = r_t1(i);

        waitbar(i/TT,h)
end
    
close(h)

chi      = 1./exp(0.5*M);    % Likelihood non normalized
max_chi  = max(chi(:));      % Max value of the Likelihood
norm_chi = chi/max_chi;      % Likelihood normalized

%% Resampling the likelihood by comparing to a random value between 0 and 1

jt=0;
for it=1:TT
    R=rand;
    if (norm_chi(it)>R)
        jt        = jt+1;
        s_chi(jt) = norm_chi(it);
        s_t(jt)   = t_vec(it);
    end
end

%% extract 1d PDFs and confidence intervals

nbin             = 20;

[n,xout]         =   hist(s_t,nbin);
xwork            =   cumsum(n/sum(n));

ix               =   find(xwork>0.175,1);
t_1sd            =   xout(ix);
t_1sd_pd         =   n(ix)/sum(n);

ix               =   find(xwork>0.025,1);
t_2sd            =   xout(ix);
t_2sd_pd         =   n(ix)/sum(n);

ix               =   find(xwork>0.825,1);
t_1su            =   xout(ix);
t_su_pd          =   n(ix)/sum(n);

ix               =   find(xwork>0.925,1);
t_2su            =   xout(ix);
t_2su_pd         =   n(ix)/sum(n);

ix               =   find(xwork>0.50,1);
t_median         =   xout(ix);
t_median_pd      =   n(ix)/sum(n);

%% Print the results

fprintf('\nResult for the inversion\n')

disp(['t Median     = ' num2str(t_median,3) ' yr']);
disp(['t 1sigma sup = ' num2str(abs(t_median-t_1su),3) ' yr']);
disp(['t 1sigma inf = ' num2str(abs(t_median-t_1sd),3) ' yr']);
disp(['t 2sigma sup = ' num2str(abs(t_median-t_2su),3) ' yr']);
disp(['t 2sigma inf = ' num2str(abs(t_median-t_2sd),3) ' yr']);

%% Plotting

xs   = 0:0.5:40;
Ls_M = exp(-SP*t_median*exp(-mu*xs));
% Ls_M = (SP.*exp(-mu.*xs).*exp(-t_median.*(SP.*exp(-mu.*xs)+magicN))+magicN)./(SP.*exp(-mu.*xs)+magicN);    % Sohbati et al. 2012a

figure(102)
set(gcf,'units','points','position',[10,1200,1000,300])

subplot(1,2,1)
plot(x,L,'go','MarkerFaceColor','g')
hold on 
plot(xs,Ls_M,'r','LineWidth',1)
xlabel('Depth [mm]')
ylabel('Normalized IRSL Signal')
legend('Experimental values','Inversed solution Median','Inversed solution Bestfit','Location','Southeast')
axis([0 30 0 1.2])
title('Evolution of the luminescence signal with depth')

t_median_vec = t_median*ones(100,1);
likeH        = 0:1/(100-1):1;

subplot(1,2,2)
%set(gca,'XScale','log')   
plot(t_vec,norm_chi,'b','LineWidth',1)
hold on
plot(t_median_vec,likeH,'r','LineWidth',1)
% axis([0 exp(tmax) 0 1])
ylabel('Likelihood')
xlabel('Time [a]')
title('Probability distribution')
legend('Likelihood distribution','Median','Location','Northeast')

