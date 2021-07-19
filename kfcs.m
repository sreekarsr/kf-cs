clear;
close all;
delete('log.txt');
diary('log.txt');
diary on;
%% Kalman Filtered Compressed Sensing
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
S_max = 25;% maximum sparsity run for 8, 16, 25
type = "unknown";
% algorithm decision parameters
FEN_thresh =100; % have to change this based on genie-aided KF (different for each Smax), and maybe some theoretical perspective?
alpha_a = 0.1;%threshold for addition (based on obtained values)

% Noise Variances
noisevar_obs = ((1/3)*sqrt(S_max/n))^2;
noisevar_init = 9;
noisevar_sys = 1;


% DS parameters
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
    T = [];% initialise support set for simulation

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

    if(type=="known")
        T=T1;
    elseif(type=="unknown")
        T = [];% initialise support set
    end
    xcap = zeros(m,1);
    disp('T1:');disp(T1);
    disp('T5:');disp(T5);
    disp('Delta:');disp(setdiff(T5,T1));
    for t=tvec
        %% KF prediction and update
        %     eq (4) and (5)
        xcap_prior = xcap;
        P_prior(T,T) = P(T,T) + noisevar_sys*eye(length(T));

        R_ie = A(:,T)*P_prior(T,T)*(A(:,T)') + noisevar_obs*eye(n); % innovation error covariance - nxn matrix
        K(T,:) = P_prior(T,T)*(A(:,T)')*(inv(R_ie));

        % assign new estimates
        xcap(T) = xcap_prior(T) +  K(T,:)*(Y(:,t) - A*xcap_prior);
        % xcap(Tc) = xcap(Tc); %redundant! for representational purposes
        P(T,T) = (eye(length(T)) - K(T,:)*A(:,T))*P_prior(T,T);

        %% Check for addition (non-zero mean for filtering error)
        % calculate filtering error norm
        filt_error = Y(:,t) - A*xcap;
        R_fe = (eye(n) - A(:,T)*K(T,:))*R_ie*(eye(n) - A(:,T)*K(T,:));
        FEN = filt_error'*inv(R_fe)*filt_error; % Filtering error norm
        fprintf('t = %d FEN : %1.5e\n',t,FEN);

        if(FEN > FEN_thresh)
            disp('FEN is greater than threshold');
            %% Addition step
            Tc = setdiff((1:m)',T);

            % (a) Run CS (Dantzig Selector)
            [betacap,iter,dval,time] = selector(A(:,Tc),ones(length(Tc),1),' ',filt_error,delta,eps,maxiter);

            % Threshold to estimate increase in support
            nz = find(betacap.^2 > alpha_a);%careful here
            Deltacap = Tc(nz);

            % b. compute T_new
            Tnew = sort(union(T,Deltacap));

            % set T = Tnew, expand P_prior for additional supports
            T = Tnew;
            disp('added to support:');
            disp(Deltacap);
            if(length(T)==length(T5))
                if(T~=T5)
                    disp('Support of same size but not matching');
                    disp('Not added yet: ');
                    if(t<5)
                        disp(setdiff(T1,T));
                    else
                        disp(setdiff(T5,T));
                    end
                    disp('Extra added:');
                    if(t<5)
                        disp(setdiff(T,T1));
                    else
                        disp(setdiff(T,T5));
                    end
                     disp('beta values')
                     disp([Tc(betacap~=0) betacap(betacap~=0).^2]);
                end
            else
                 disp('Size not matching'); 
                    disp('Not added yet: ');
                    if(t<5)
                        disp(setdiff(T1,T));
                    else
                        disp(setdiff(T5,T));
                    end
                    disp('Extra added:');
                    if(t<5)
                        disp(setdiff(T,T1));
                    else
                        disp(setdiff(T,T5));
                    end
                    disp('beta values')
                    disp([Tc(betacap~=0) betacap(betacap~=0).^2]);
            end
            
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
            FEN = filt_error'*inv(R_fe)*filt_error; % Filtering error norm
            fprintf('\nt = %d FEN : %1.5e  (new support)\n',t,FEN);
            
            % Deletion step (can try implementing as an exercise to show you did something extra)

        end
           
        MSE_vec(k,t) = norm(X(:,t)-xcap)^2; % should change for Monte Carlo
    end
     disp('10 steps complete')
            if(length(T)==length(T5))
                if(T~=T5)
                    disp('Support of same size but not matching');
                     disp('Not added yet: ');
                    if(t<5)
                        disp(setdiff(T1,T));
                    else
                        disp(setdiff(T5,T));
                    end
                    disp('Extra added:');
                    if(t<5)
                        disp(setdiff(T,T1));
                    else
                        disp(setdiff(T,T5));
                    end
                end
            else
                 disp('Size not matching'); 
                   disp('Not added yet: ');
                    if(t<5)
                        disp(setdiff(T1,T));
                    else
                        disp(setdiff(T5,T));
                    end
                    disp('Extra added:');
                    if(t<5)
                        disp(setdiff(T,T1));
                    else
                        disp(setdiff(T,T5));
                    end
            end
 
end

figure;
plot(tvec,MSE_vec);

MSE_vec_avg = sum(MSE_vec,1)./Niter;
figure;
plot(tvec,MSE_vec_avg);
ylim([0,30]);
% ylim([0 0.4]);
% ylim([0 5]);
xlim([1 10]);

save(sprintf("MSE_vec_kfcs_%s_%d.mat",type,S_max),'MSE_vec_avg')

diary off;
