%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  code  : Inversion_two_variables_t_unknown_resampling.m
% Version     : 20.01.2020
% Input       : Unknown exposure age of the samples [year]
%               Experimental luminescence signal (Lx/Tx normalized)
%               Depth related to Exp. Lumi. Signal [mm]
% Output      : Bleaching rate SP x exposure time 
%               Attenuation coeff. mu [mm-1]
% Inversion   : L1-norm waited over the experimental noise of the luminescence plateau
% Resampling  : Resampling of the likelihood by comparing to a random value
%               between 0 and 1
% Units       : SP  [a-1]
%             : mu  [mm-1]
%             : Age [a]
%             : Ddot_input [Gy ka-1] conversion to [Gy a-1]
%             : D0  [Gy]
%             : Magic number = Ddot/D0 [a-1]
% IMPORTANT   : Some variable have to be adapted for every different dataset
%               >> nC1 and nC2 are the last disc number for respectively core
%                  1 and core 2, these variables are just for plotting
%               >> ind_plt defines at which index the plateaus starts considering all cores as one single dataset
%                  can be determine by exploring LxTx_s
% Contact     : B.Lehmann (lehmann.benj@gmail)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
close all;
addpath(genpath('/Users/blehman1/Documents/0_Pro/5_Code'));  % Set the depository path to the actual working folder

TT      = 1000;                                              % Size of the inversion matrix(TT,TT)

%% Loading data file excel MBMV_.xlsx 

[num,txt,tab]    =  xlsread('OSL_MBTP7_cal.xlsx');
n                =  length(num(:,1));

%% Ordering of the imported data

nC1          = 23;                     % Disc number for core 1
nC2          = 47;                     % Disc number for core 2
ind_plt      = 40;                     % at which index the plateaus starts considering all cores as one single dataset

x(1:n)       = (num(1:n,1));           % 1st colomn of excell table = depth [mm]
LxTx(1:n)    = (num(1:n,2));           % 2nd colomn of excell table = Lx/Tx
LxTx_e(1:n)  = (num(1:n,3));           % 3rd colomn of excell table = Error Lx/Tx

[x_s,ind]    = sort(x(:));             % sorting increasingly the depth 
LxTx_s       = LxTx(ind);              % sorting Lx/Tx with increasing depth
LxTx_a       = std(LxTx_s(ind_plt:n)); % calculating the standard deviation on the plateau

%% Definition of parameter and domain of freedom 

x_model      = 0:0.25:25;           % create a synthetic depth version [mm]

log_SP_t_max = 4; 
log_SP_t_min = 1;     
SP_t_max     = 10^(log_SP_t_max);   % Dimensionless
SP_t_min     = 10^(log_SP_t_min);   % Dimensionless
         
mu_max       = 2.0;                 % mm-1
mu_min       = 0.1;                 % mm-1

%% Initialization of the variable and matrix for the inversion

M            = nan(TT,TT);    % creation of an empty misfit matrix

mu_matrix    = nan(TT,TT);    % creation of an empty mu matrix
SP_t_matrix  = nan(TT,TT);    % creation of an empty SP_t matrix

rand_vec1    = rand(TT,1);    % random vector of size TT to sample SP_t
rand_vec2    = rand(TT,1);    % random vector of size TT to sample mu

r_SP_t       = sort(10.^(log_SP_t_min+(log_SP_t_max-log_SP_t_min)*rand_vec1));    % r_SP_t is randomly sampled in log space
r_mu         = sort(mu_min+(mu_max-mu_min)*rand_vec2);                            % mu is randomly sampled in normal space

%% Inversion

h               = waitbar(0,'Put a good song and relax');  % creation of the loading bar (just for fun)

    for i = 1:TT
        for j = 1:TT                                       % "parfor" for paralellelization, for normal claculation put "for"
            
            M(i,j)   = 0;                                  % initialization of the misfit at i and j index
                                
            LxTx_th  = exp(-r_SP_t(j)*exp(-r_mu(i)*x));    % creation of synthetic luminescence signal
            M(i,j)   = sum((abs(LxTx-LxTx_th))./LxTx_a);   % calculation of the misfit between experimental data and synthetic signal                             
                                    
            mu_matrix(i,j)  = r_mu(i);                     % compile mu tested into matrix
            SP_t_matrix(i,j)  = r_SP_t(j);                 % compile SP_t tested into matrix
            
        end
        waitbar(i/TT,h)                                    % incrementation of the loading bar (just for fun)
    end

close(h)                                                   % close the loading bar

%% Transformation of the misfit M into likelihood chi and normalization

chi        = 1./exp(0.5*M);

max_chi    = max(chi(:));
norm_chi   = chi./max_chi;

%% Resampling the likelihood by comparing to a random value between 0 and 1

jt=0; 
for it=1:TT
    for j= 1:TT
    R=rand;
    if (norm_chi(it,j)>R)
        jt=jt+1;
        s_chi(jt)= norm_chi(it,j); 
        s_SP_t(jt) = SP_t_matrix(it,j);
        s_mu(jt) = mu_matrix(it,j);
    end
    end 
end

%% Extract 1d PDFs and confidence intervals for MBTP8

nbin                 = 20;
% For SP
[n_SP_t,xout_SP_t]   =   hist(log10(s_SP_t),nbin);
xwork_SP_t           =   cumsum(n_SP_t/sum(n_SP_t));

% -1 sigma 
ix_SP_t              =   find(xwork_SP_t>0.175,1);
SP_t_1sd               =   xout_SP_t(ix_SP_t);

% +1 sigma 
ix_SP_t                =   find(xwork_SP_t>0.825,1);
SP_t_1su               =   xout_SP_t(ix_SP_t);

% -2 sigma 
ix_SP_t                =   find(xwork_SP_t>0.025,1);
SP_t_2sd               =   xout_SP_t(ix_SP_t);

% +2 sigma 
ix_SP_t                =   find(xwork_SP_t>0.925,1);
SP_t_2su               =   xout_SP_t(ix_SP_t);

% Median
ix_SP_t                =   find(xwork_SP_t>0.50,1);
SP_t_M                 =   xout_SP_t(ix_SP_t);

%% For mu

[n_mu,xout_mu]           =   hist(s_mu,nbin);
xwork_mu                 =   cumsum(n_mu/sum(n_mu));

% -1 sigma 
ix_mu                    =   find(xwork_mu>0.175,1);
mu_1sd                   =   xout_mu(ix_mu);

% +1 sigma 
ix_mu                    =   find(xwork_mu>0.825,1);
mu_1su                   =   xout_mu(ix_mu);

% +2 sigma 
ix_mu                    =   find(xwork_mu>0.025,1);
mu_2sd                   =   xout_mu(ix_mu);

% +2 sigma 
ix_mu                    =   find(xwork_mu>0.925,1);
mu_2su                   =   xout_mu(ix_mu);

% Median
ix_mu                    =   find(xwork_mu>0.50,1);
mu_M                     =   xout_mu(ix_mu);

%% Print the results

fprintf('\nResult for the calibration with only MBTP8 \n')

SP_t_M_ok   = 10^(SP_t_M);
SP_t_1su_ok = 10^(SP_t_1su);
SP_t_1sd_ok = 10^(SP_t_1sd);
SP_t_2su_ok = 10^(SP_t_2su);
SP_t_2sd_ok = 10^(SP_t_2sd);

disp(['SPxt Median     = ' num2str(SP_t_M_ok,3) ]);
disp(['SPxt 1sigma sup = ' num2str(abs(SP_t_M_ok-SP_t_1su_ok),3)]);
disp(['SPxt 1sigma inf = ' num2str(abs(SP_t_M_ok-SP_t_1sd_ok),3)]);
disp(['SPxt 2sigma sup = ' num2str(abs(SP_t_M_ok-SP_t_2su_ok),3)]);
disp(['SPxt 2sigma inf = ' num2str(abs(SP_t_M_ok-SP_t_2sd_ok),3)]);

disp(['mu Median     = ' num2str(mu_M,3) ' mm-1']);
disp(['mu 1sigma sup = ' num2str(abs(mu_M-mu_1su),3) ' mm-1']);
disp(['mu 1sigma inf = ' num2str(abs(mu_M-mu_1sd),3) ' mm-1']);
disp(['mu 2sigma sup = ' num2str(abs(mu_M-mu_2su),3) ' mm-1']);
disp(['mu 2sigma inf = ' num2str(abs(mu_M-mu_2sd),3) ' mm-1']);


%% Creating output models

LxTx_M   = exp(-SP_t_M_ok*exp(-mu_M*x_model));

%% Plotting

figure('NumberTitle','off','Position',[00 00 1200 350])
sl =5;

subplot(1,2,1)
errorbar(x(1:nC1)     ,LxTx(1:nC1)     ,LxTx_e(1:nC1)     ,'go' ,'markersize', sl)
hold on
errorbar(x(nC1+1:nC2),LxTx(nC1+1:nC2),LxTx_e(nC1+1:nC2),'b>' ,'markersize', sl)
errorbar(x(nC2+1:n)  ,LxTx(nC2+1:n)  ,LxTx_e(nC2+1:n)  ,'rd' ,'markersize', sl)
plot(x_model,LxTx_M ,'k', 'LineWidth',1)

xlabel('Depth [mm]')
ylabel('IRSL Normalized Intensity')
legend('Core 1','Core 2','Core 3','Infered model with Median','Location','Southeast')
title('Evolution of the luminescence signal with depth')
axis([0 28 0 1.3])

subplot(1,2,2)
surface(r_SP_t,r_mu,norm_chi); axis square; colorbar;title(colorbar,'Likelihood');shading interp;set(gca,'XScale','log');
hold on

SP_vec_plot      = linspace(mu_min,mu_max,TT);
SP_M_line_plot   = SP_t_M_ok.*ones(size(SP_vec_plot));

mu_vec_plot      = linspace(SP_t_min,SP_t_max,TT); 
mu_M_line_plot   = mu_M.*ones(size(mu_vec_plot));

plot(mu_vec_plot,mu_M_line_plot,'w','LineWidth',1)
plot(SP_M_line_plot,SP_vec_plot,'w','LineWidth',1)

xlabel('Bleaching rate x time [Dimensionless]')
ylabel('Attenuation coeff. [mm-1]')
axis([SP_t_min SP_t_max mu_min mu_max])
title('Probability distribution')

