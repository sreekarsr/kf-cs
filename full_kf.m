clear;
close all;
%% Full Kalman Filter
% Course Project : EE6110 Adaptive Signal Processing
% EE18B154 Sreekar Sai Ranganathan

addpath ADMDS1.1; % for Dantzig Selector
% Credits : YhongZhangAI

%% Parameter Initialisation
% Model Dimensions
m = 256; % length of signal (x_t)
n = 72; % length of observation (y_t)

% Model matrix
A = randn(n,m); % n × m i.i.d. Gaussian entries
A = normc(A);% normalise each column of A

% Maximum Sparsity of x_t
S_max = 8;% maximum sparsity run for 8, 16, 25

% Noise Variances
noisevar_obs = ((1/3)*sqrt(S_max/n))^2;
noisevar_init = 9;
noisevar_sys = 1;

Niter =100 ; % no of monte carlo simulations

tvec = 1:1:10;

X = NaN(m,length(tvec));
Y = NaN(n,length(tvec));

MSE_vec = zeros(1,length(tvec));

for k = 1:Niter
    % Support sets
    T1 = sort(randperm(m, S_max - 2)'); % initial support set (till T4)
    T1c = setdiff((1:m)',T1);
    T5 = sort([T1; T1c(randperm(length(T1c),2)')]); % final support set (for 5 onwards)

    fprintf('\nIteration:%d\n',k);
    x = zeros(m,1);
    T = [];% initialise support set
    
    %% Simulation (ground truth)
    for t=tvec
        Tlast = T;
        if t==1
            T = T1;
        elseif t==5
            T = T5;
        end
        Delta = setdiff(T,Tlast);

        v = zeros(m,1);
        v(Delta) = sqrt(noisevar_init)*randn(length(Delta),1);
        v(Tlast) = sqrt(noisevar_sys)*randn(length(Tlast),1);
        % NOTE : NOT SIMULATING DELETION

        x = x + v;
        w = sqrt(noisevar_obs)*randn(n,1);
        y = A(:,T)*x(T) + w;

        % store the pair
        X(:,t) = x;
        Y(:,t) = y;
    end
    clear x y w v;

    %% Full KF
    P = zeros(m,m); % P=0 (initialising largest possible matrix)
    P_prior = NaN(m,m);
    K = NaN(m,n); % 

    T = 1:m;
    xcap = zeros(m,1);

    for t=tvec
        %% KF prediction and update
        %     eq (4) and (5)
        % TODO : TAKE CARE OF ANY HANGING INDICES (not updated during these steps)

        xcap_prior = xcap;
        P_prior(T,T) = P(T,T) + noisevar_sys*eye(length(T));

        R_ie = A(:,T)*P_prior(T,T)*(A(:,T)') + noisevar_obs*eye(n); % innovation error covariance - nxn matrix
        K(T,:) = P_prior(T,T)*(A(:,T)')*(inv(R_ie));

        Tc = setdiff((1:m)',T);
        % assign new estimates
        xcap(T) = xcap_prior(T) +  K(T,:)*(Y(:,t) - A*xcap_prior);
        % xcap(Tc) = xcap(Tc); %redundant! for representational purposes
        P(T,T) = (eye(length(T)) - K(T,:)*A(:,T))*P_prior(T,T);


        MSE_vec(t) = MSE_vec(t) + norm(X(:,t)-xcap)^2; % should change for Monte Carlo
    end
 
end
MSE_vec = MSE_vec/Niter;
save("MSE_fullKF_8.mat",'MSE_vec')
figure;
plot(tvec,MSE_vec);
% ylim([0,30]);
