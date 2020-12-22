clear;
close all;
%% Kalman Filtered Compressed Sensing
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

% algorithm parameters
FEN_thresh = 7e-2; % have to change this based on genie-aided KF (different for each Smax), and maybe some theoretical perspective?
lambda_m = sqrt(2*log(m));
delta = lambda_m*sqrt(noisevar_obs);
alpha_a = 1.0;%threshold for addition (based on obtained values)

% DS parameters
eps = 1e-3;%TODO : see what this value means - tolerance for alternating direction method
maxiter = 1e3; %TODO : check reasonable estimate for this

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

    %% KF-CS Algorithm (known T1)
    % Smax order filter, with T1 and T5 known
    P = zeros(m,m); % P=0 (initialising largest possible matrix)
    P_prior = NaN(m,m);
    K = NaN(m,n);

    T = T1;
    xcap = zeros(m,1);

    for t=tvec
        %% KF prediction and update
        %     eq (4) and (5)
        xcap_prior = xcap;
        P_prior(T,T) = P(T,T) + noisevar_sys*eye(length(T));

        R_ie = A(:,T)*P_prior(T,T)*(A(:,T)') + noisevar_obs*eye(n); % innovation error covariance - nxn matrix
        K(T,:) = P_prior(T,T)*(A(:,T)')*(inv(R_ie));

        Tc = setdiff((1:m)',T);
        % assign new estimates
        xcap(T) = xcap_prior(T) +  K(T,:)*(Y(:,t) - A*xcap_prior);
        % xcap(Tc) = xcap(Tc); %redundant! for representational purposes
        P(T,T) = (eye(length(T)) - K(T,:)*A(:,T))*P_prior(T,T);

        %% Check for addition (non-zero mean for filtering error)
        % calculate filtering error norm
        filt_error = Y(:,t) - A*xcap;
        R_fe = (eye(n) - A(:,T)*K(T,:))*R_ie*(eye(n) - A(:,T)*K(T,:));
        FEN = filt_error'*R_fe*filt_error; % Filtering error norm
        fprintf('\nt = %d FEN : %1.5e',t,FEN);

        if(length(T)==S_max)
            disp('Reached S_max... Not updating support');

            if(FEN > FEN_thresh)
                warning('FEN above threshold while Smax already reached!');
            end

        elseif(FEN > FEN_thresh)
            %% Addition step
            
            % (a) Run CS (Dantzig Selector)
            [betacap,iter,dval,time] = selector(A(:,Tc),ones(length(Tc),1),' ',filt_error,delta,eps,maxiter);

            % Threshold to estimate increase in support
            nz = find(betacap.^2 > alpha_a);%careful here
            Deltacap = Tc(nz);

            % b. compute T_new
            Tnew = sort(union(T,Deltacap));

            % set T = Tnew, expand P_prior for additional supports
            T = Tnew;
            P_prior(Deltacap,Deltacap) = noisevar_init*eye(length(Deltacap));
            P_prior(Deltacap,sort(setdiff(T,Deltacap))) = 0;
            P_prior(sort(setdiff(T,Deltacap)),Deltacap) = 0;

            % Rerun KF update (5)         
            R_ie = A(:,T)*P_prior(T,T)*(A(:,T)') + noisevar_obs*eye(n); % innovation error covariance - nxn matrix
            K(T,:) = P_prior(T,T)*(A(:,T)')*(inv(R_ie));

            Tc = setdiff((1:m)',T);
            % assign new estimates
            xcap(T) = xcap_prior(T) +  K(T,:)*(Y(:,t) - A*xcap_prior);
            % xcap(Tc) = xcap(Tc); %redundant! for representational purposes
            P(T,T) = (eye(length(T)) - K(T,:)*A(:,T))*P_prior(T,T);

            % Post-addition FEN     
            filt_error = Y(:,t) - A*xcap;
            R_fe = (eye(n) - A(:,T)*K(T,:))*R_ie*(eye(n) - A(:,T)*K(T,:));
            FEN = filt_error'*R_fe*filt_error; % Filtering error norm
            fprintf('\nt = %d FEN : %1.5e  (new support)',t,FEN);
            
            % Deletion step (can try implementing as an exercise to show you did something extra)

        end

        MSE_vec(t) = MSE_vec(t) + norm(X(:,t)-xcap)^2; % should change for Monte Carlo
    end
 
end
MSE_vec = MSE_vec/Niter;
save("MSE_kfcs_known_8.mat",'MSE_vec')
figure;
plot(tvec,MSE_vec);
ylim([0,0.4]);

%% Normal CS (Using DS?) Threshold values?


%% Full KF