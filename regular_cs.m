clear;
close all;
%% Regular Compressed Sensing
% Course Project : EE6110 Adaptive Signal Processing
% EE18B154 Sreekar Sai Ranganathan

addpath ADMDS1.1; % for Dantzig Selector
% Credits : YongZhangAI

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

% DS parameters TODO : have to change these for this one
lambda_m = sqrt(2*log(m));
delta = lambda_m*sqrt(noisevar_obs);
eps = 1e-3;
maxiter = 100; % usually around 50 iterations sufficient

Niter =100 ; % no of monte carlo simulations

tvec = 1:1:10;

X = NaN(m,length(tvec));
Y = NaN(n,length(tvec));

MSE_vec = zeros(Niter,length(tvec));

for k = 1:Niter
    % Support sets
    T1 = sort(randperm(m, S_max - 2)'); % initial support set (till T4)
    T1c = setdiff((1:m)',T1);
    T5 = sort([T1; T1c(randperm(length(T1c),2)')]); % final support set (for 5 onwards)

    fprintf('\nIteration:%d\n',k);
    x = zeros(m,1);
    T=[];
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
        
        x = x + v;
        w = sqrt(noisevar_obs)*randn(n,1);
        y = A(:,T)*x(T) + w;

        % store the pair
        X(:,t) = x;
        Y(:,t) = y;
    end
    clear x y w v;

    %%  Regular CS (Dantzig Selector)
    for t=tvec
            % Regular CS (Dantzig Selector)
            [xcap,iter,dval,time] = selector(A,ones(m,1),' ',Y(:,t),delta,eps,maxiter);
            nzidx = sum(abs([X(:,t) xcap]),2)~=0;
            disp('idx x xcap');disp([find(nzidx) X(nzidx,t) xcap(nzidx)]);
            MSE_vec(k,t) = norm(X(:,t)-xcap)^2; % should change for Monte Carlo
    end
 
end

figure;
plot(tvec,MSE_vec);

MSE_vec_avg = sum(MSE_vec,1)./Niter;
figure;
plot(tvec,MSE_vec_avg);
% ylim([0,0.4]);

save(sprintf("MSE_vec_regcs_%d.mat",S_max),'MSE_vec_avg');
