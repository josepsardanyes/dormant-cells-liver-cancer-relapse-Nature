%% =========================================================================
%  ABC-SMC (Approximate Bayesian Computation - Sequential Monte Carlo)
%  Distance metric: Euclidean distance between simulated and observed data
%
%  Algorithm: Toni et al. (2009) J. R. Soc. Interface
%
%  Structure:
%    1. User-defined section  – ODE, priors, tolerances, data
%    2. ABC-SMC engine        – population iterations
%    3. Results & diagnostics – plots and parameter estimates
% =========================================================================

clear; clc; close all;
rng(42);   % reproducibility

%% =========================================================================
%  SECTION 1 – USER CONFIGURATION
% =========================================================================

% ---- 1a. Experimental data -----------------------------------------------
%  t_obs  : observation times  (column vector, length T)
%  Y_obs  : observed tumour volume (column vector, length T)

script_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(script_dir);
sample_id = 3; % Select experimental sample: 1, 2, 3, or 4
if ~ismember(sample_id, 1:4)
    error('sample_id must be 1, 2, 3, or 4.');
end
sample_label = sprintf('sample%d', sample_id);
data_file = fullfile(repo_dir, 'data', ...
    sprintf('volume_treated_%d_NoZeros.dat', sample_id));
observed_data = load(data_file);
t_obs = observed_data(:,1);
Y_obs = observed_data(:,4);
results_dir = fullfile(repo_dir, 'results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

[T, n_observed_states] = size(Y_obs); %#ok<ASGLU> % one observed variable

% ---- 1b. Parameter bounds (uniform prior) --------------------------------
%  Rows = parameters 1..9.  Columns = [lower, upper].
%  Adjust to your biological / physical knowledge.

param_names = {'p1','p2','p3','p4','p5','p6','p7','p8','p9'};
prior_lb = [0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0];
%prior_ub = [0.3, 0.02, 0.1, 0.3, 50, 0.1, 50, 0.3, 0.5];
%prior_ub = [0.3, 0.02, 0.1, 0.3, 0.1, 0.1, 0.1, 0.3, 0.5];
prior_ub = [0.3, 0.02, 0.1, 0.3, 100, 0.1, 100, 0.3, 0.5];

n_params  = 9;

% ---- 1c. ABC-SMC settings ------------------------------------------------
N          = 300;           % particles per population
%n_pop      = 20;             % number of SMC populations (iterations)
n_pop      = 25;             % number of SMC populations (iterations)
p_acc      = 0.5;           % target acceptance rate for adaptive epsilon
epsilon_0  = 1e0;           % first threshold (accept all, auto-set below)
max_tries  = 1e6;           % safety cap on simulation attempts per particle

% Perturbation kernel bandwidth (component-wise Gaussian, scaled each pop)
kernel_factor = 0.5;        % multiply component std by this factor

% ODE solver options
ode_opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.5);

% Initial conditions for the ODE  (length = n_states = 5)
y0 = [80; 0.0; 0.0; 0.0; 0.0];   % <-- set to your system's ICs

%% =========================================================================
%  SECTION 2 – ABC-SMC ENGINE
% =========================================================================

fprintf('=== ABC-SMC  |  %d particles  |  %d populations ===\n\n', N, n_pop);

% Storage
particles  = zeros(N, n_params, n_pop);   % accepted particles
weights    = zeros(N, n_pop);             % importance weights
distances  = zeros(N, n_pop);            % Euclidean distances
epsilons   = zeros(1, n_pop);            % thresholds used

% ------------------------------------------------------------------
%  Population 0  – sample directly from prior, keep best N distances
% ------------------------------------------------------------------
fprintf('Population 1/%d  (sampling from prior) ...\n', n_pop);

raw_theta = zeros(max_tries, n_params);
raw_dist  = zeros(max_tries, 1);
n_sim     = 0;

for k = 1:max_tries
    %fprintf('k = %d\n',k);
    theta_k = sample_prior(prior_lb, prior_ub, n_params);
    [d,TotalQ,RootMeanSqErr,volume] = simulate_and_distance(theta_k, t_obs, y0, Y_obs, ode_opts);
    %fprintf('k = %d distance = %d \n',k,d);
    if d<epsilon_0
        n_sim = n_sim + 1;
        raw_theta(n_sim,:) = theta_k;
        raw_dist(n_sim)    = d;
        if n_sim >= 5*N, break; end   % collect pool, then trim
        %figure(1)
        %plot(t_obs,volume)
        %hold on
    end
    fprintf('n_sim = %d\n',n_sim);
end

raw_theta = raw_theta(1:n_sim,:);
raw_dist  = raw_dist(1:n_sim);

% Sort and keep N best
[raw_dist_sorted, idx] = sort(raw_dist);
keep = idx(1:min(N, n_sim));

particles(:,:,1) = raw_theta(keep,:);
distances(:,1)   = raw_dist(keep);
weights(:,1)     = 1/N;
epsilons(1)      = raw_dist_sorted(min(N, n_sim));

fprintf('  epsilon_1 = %.4f  |  acceptance rate ~ %.1f%%\n\n', ...
        epsilons(1), 100*n_sim/max_tries);

% ------------------------------------------------------------------
%  Populations 1..n_pop-1  – SMC iterations
% ------------------------------------------------------------------
for pop = 2:n_pop

    prev_particles = particles(:,:,pop-1);
    prev_weights   = weights(:,pop-1);

    % Compute perturbation kernel bandwidth (2 * weighted variance)
    wmu  = sum(prev_weights .* prev_particles, 1);
    wvar = sum(prev_weights .* (prev_particles - wmu).^2, 1);
    sigma_k = kernel_factor * sqrt(2 * wvar + 1e-12);   % component std

    % Adaptive epsilon: target p_acc quantile of previous distances
    epsilon = quantile(distances(:,pop-1), p_acc);
    epsilons(pop) = epsilon;

    fprintf('Population %d/%d  |  epsilon = %.4f ...\n', pop, n_pop, epsilon);

    new_particles = zeros(N, n_params);
    new_distances = zeros(N, 1);
    n_accepted    = 0;
    n_tried       = 0;
    
    QuiescentCells=zeros(1,N);
    RootMeanSquaredError=zeros(1,N);
    while n_accepted < N && n_tried < max_tries

        % 1. Sample ancestor from previous population
        anc_idx = randsample(N, 1, true, prev_weights);
        theta_star = prev_particles(anc_idx,:);

        % 2. Perturb with Gaussian kernel
        theta_prop = theta_star + sigma_k .* randn(1, n_params);

        % 3. Enforce prior bounds (reject if outside)
        if any(theta_prop < prior_lb) || any(theta_prop > prior_ub)
            n_tried = n_tried + 1;
            continue
        end

        % 4. Simulate and compute distance
        [d,TotalQ,RootMeanSqErr,volume] = simulate_and_distance(theta_prop, t_obs, y0, Y_obs, ode_opts);
        n_tried = n_tried + 1;

        if isfinite(d) && d <= epsilon
            %fprintf('epsilon = %d \n',epsilon)
            n_accepted = n_accepted + 1;
            if pop == n_pop && n_accepted > 0
               QuiescentCells(1,n_accepted) = TotalQ;
               RootMeanSquaredError(1,n_accepted) = RootMeanSqErr;
            end 
            new_particles(n_accepted,:) = theta_prop;
            new_distances(n_accepted)   = d;
            if pop == n_pop
              figure(1)
              plot(t_obs,volume)
              hold on
            end
        end
        fprintf('n_accepted in pop %d = %d\n',pop,n_accepted);        
    end

    if n_accepted < N
        warning('Only %d/%d particles accepted in population %d.', ...
                n_accepted, N, pop);
        new_particles = new_particles(1:n_accepted,:);
        new_distances = new_distances(1:n_accepted);
        % Pad if needed (rare edge case)
        extra = N - n_accepted;
        new_particles = [new_particles; new_particles(end-extra+1:end,:)];
        new_distances = [new_distances; new_distances(end-extra+1:end)];
    end

    % 5. Compute importance weights
    log_prior   = zeros(N, 1);   % log p(theta) – uniform => constant, cancels
    log_kern_sum = zeros(N, 1);

    for i = 1:N
        log_k = zeros(N,1);
        for j = 1:N
            log_k(j) = log(prev_weights(j) + 1e-300) + ...
                       sum(-0.5*((new_particles(i,:) - prev_particles(j,:))./sigma_k).^2 ...
                           - log(sigma_k*sqrt(2*pi)));
        end
        log_kern_sum(i) = logsumexp_vec(log_k);
    end

    log_w = log_prior - log_kern_sum;
    log_w = log_w - max(log_w);
    w     = exp(log_w);
    w     = w / sum(w);

    particles(:,:,pop) = new_particles;
    distances(:,pop)   = new_distances;
    weights(:,pop)     = w;

    eff_N = 1 / sum(w.^2);
    fprintf('epsilon %d  |  Accepted %d  |  attempts %d  |  acc_rate %.1f%%  |  ESS %.0f\n\n',...
            epsilon, n_accepted, n_tried, 100*n_accepted/n_tried, eff_N);
       
end

QuiescentCells

figure(1)
gscatter(t_obs,Y_obs)
hold on

figure(2)
boxplot(QuiescentCells)

figure(3)
boxplot(distances(:,n_pop))

figure(4)
boxplot(RootMeanSquaredError)

%% =========================================================================
%  SECTION 3 – RESULTS & DIAGNOSTICS
% =========================================================================

format shortEng;

final_particles = particles(:,:,end);
final_weights   = weights(:,end);

fprintf('=== Posterior summary (final population) ===\n');
fprintf('%-6s  %8s  %8s  %8s  %8s\n','Param','Mean','Std','2.5%','97.5%');
for k = 1:n_params
    vals = final_particles(:,k);
    w    = final_weights;
    pmean = sum(w .* vals);
    pstd  = sqrt(sum(w .* (vals - pmean).^2));
    sorted_vals = sort(vals);
    p025 = sorted_vals(round(0.025*N));
    p975 = sorted_vals(round(0.975*N));
    fprintf('%-6s  %8.4f  %8.4f  %8.4f  %8.4f\n', ...
            param_names{k}, pmean, pstd, p025, p975);
end

% ---- Plot 1: Posterior marginals -----------------------------------------
figure('Name','Posterior Marginals','Position',[50 50 1200 600]);
for k = 1:n_params
    subplot(3,4,k);
    %histogram(final_particles(:,k), 30, 'Normalization','pdf', ...
    %          'FaceColor',[0.2 0.5 0.8], 'EdgeColor','none','FaceAlpha',0.7);
    ecdf(final_particles(:,k),'Bounds','on');
    hold on;
    PriorSample = prior_lb(k)+(prior_ub(k)-prior_lb(k))*rand(1,1000);
    ecdf(PriorSample);
    xlabel(param_names{k});
    ylabel('Empirical CDF');
    title(['Posterior: ' param_names{k}]);
    legend('Posterior','Prior')
    box off;
end
sgtitle('ABC-SMC Posterior Marginals (final population)','FontWeight','bold');

% ---- Plot 2: Epsilon schedule --------------------------------------------
figure('Name','Epsilon Schedule','Position',[100 100 500 350]);
plot(1:n_pop, epsilons, 'o-','LineWidth',2,'MarkerSize',7,'Color',[0.8 0.2 0.2]);
xlabel('Population'); ylabel('\epsilon');
title('Tolerance schedule'); grid on; box off;

% ---- Plot 3: Posterior predictive check ----------------------------------
figure('Name','Posterior Predictive','Position',[150 150 1100 700]);
t_fine = linspace(t_obs(1), t_obs(end), 200)';
state_names = {'State 1','State 2','State 3','State 4','State 5'};

n_samples = min(200, N);
sample_idx = randsample(N, n_samples, true, final_weights);

Y_pred = zeros(length(t_fine), 5, n_samples);
for s = 1:n_samples
    theta_s = final_particles(sample_idx(s),:);
    try
        [~, Ys] = ode45(@(t,y) ode_system(t,y,theta_s), t_fine, y0, ode_opts);
        Y_pred(:,:,s) = Ys;
    catch
        Y_pred(:,:,s) = NaN;
    end
end

for st = 1:n_observed_states
    %subplot(2,3,st);
    pred_st = squeeze(Y_pred(:,1,:)+Y_pred(:,2,:)+Y_pred(:,3,:));
    q05 = quantile(pred_st, 0.05, 2);
    q50 = quantile(pred_st, 0.50, 2);
    q95 = quantile(pred_st, 0.95, 2);

    fill([t_fine; flipud(t_fine)], [q05; flipud(q95)], ...
         [0.6 0.8 1.0], 'EdgeColor','none','FaceAlpha',0.5); hold on;
    plot(t_fine, q50, 'b-', 'LineWidth',1.5);
    plot(t_obs,  Y_obs(:), 'ro', 'MarkerSize',7, 'LineWidth',1.5);
    xlabel('Time'); ylabel('Value');
    title(state_names{st});
    legend('90% CI','Median','Data','Location','best');
    box off;
end
sgtitle('Posterior Predictive Check','FontWeight','bold');

% ---- Plot 4: Pairwise scatter (first 4 params) ---------------------------
figure('Name','Pairwise Posterior','Position',[200 200 800 700]);
%n_show = min(4, n_params);
n_show = n_params;
idx_s  = randsample(N, min(300,N), true, final_weights);
for r = 1:n_show
    for c = 1:n_show
        subplot(n_show, n_show, (r-1)*n_show+c);
        if r == c
            histogram(final_particles(idx_s,r), 20, ...
                      'FaceColor',[0.2 0.5 0.8],'EdgeColor','none');
        else
            scatter(final_particles(idx_s,c), final_particles(idx_s,r), ...
                    8, [0.2 0.5 0.8], 'filled', 'MarkerFaceAlpha', 0.5);
        end
        if r == n_show, xlabel(param_names{c}); end
        if c == 1,      ylabel(param_names{r}); end
        box off;
    end
end
sgtitle('Pairwise Posterior Scatter (first 4 params)','FontWeight','bold');

fprintf('\nDone. All figures displayed.\n');

% ---- Plot 5: Correlation with quiescence ---------------------------------
figure('Name','Correlation with quiescence','Position',[50 50 1200 600]);
for k = 1:n_params
    subplot(3,4,k);
    %histogram(final_particles(:,k), 30, 'Normalization','pdf', ...
    %          'FaceColor',[0.2 0.5 0.8], 'EdgeColor','none','FaceAlpha',0.7);
    gscatter(final_particles(:,k),QuiescentCells(:))
    xlabel(param_names{k});
    ylabel('Quiescence');
    title(['Quiescence: ' param_names{k}]);
    box off;
end
sgtitle('ABC-SMC Quiescence vs params. (final population)','FontWeight','bold');

writematrix(QuiescentCells,fullfile(results_dir,sprintf('quiescentcells_%s.dat',sample_label)));
writematrix(distances(:,n_pop),fullfile(results_dir,sprintf('distances_%s.dat',sample_label)));
writematrix(RootMeanSquaredError,fullfile(results_dir,sprintf('RootMeanSquaredError_%s.dat',sample_label)));
writematrix(final_particles,fullfile(results_dir,sprintf('final_particles_%s.dat',sample_label)));

% ---- Plot 6: Time series ---------------------------------

figure('Name','Time course with average params.','Position',[50 50 1200 600]);
for k = 1:n_params
   average(k)=mean(final_particles(:,k)); 
end
[~, Ys] = ode15s(@(t,y) ode_system(t,y,average), t_obs, y0, ode_opts);
plot(t_obs,Ys(:,1)+Ys(:,2)+Ys(:,3))
hold on
gscatter(t_obs,Y_obs)
xlabel('Time');
ylabel('Volume');
title('Volume (best-fitting parameters)');
box off;

figure('Name','Time course with best params.','Position',[50 50 1200 600]);
subplot(3,1,1)
[minepsilon,minepsilonindex]=min(distances(:,n_pop));
[~, Ys] = ode15s(@(t,y) ode_system(t,y,final_particles(minepsilonindex,:)), t_obs, y0, ode_opts);
plot(t_obs,Ys(:,1)+Ys(:,2)+Ys(:,3))
hold on
gscatter(t_obs,Y_obs)
xlabel('Time');
ylabel('Volume');
title('Volume (average parameters)');
box off;
subplot(3,1,2)
plot(t_obs,Ys(:,1:3))
subplot(3,1,3)
plot(t_obs,Ys(:,2))

writematrix(final_particles(minepsilonindex,:),fullfile(results_dir,sprintf('BestParsSample%d.dat',sample_id)));

figure('Name','Time course for minum quiescence','Position',[50 50 1200 600]);

[minQ,minQindex]=min(QuiescentCells);
[~, Ys] = ode15s(@(t,y) ode_system(t,y,final_particles(minQindex,:)), t_obs, y0, ode_opts);
subplot(2,1,1)
plot(t_obs,Ys(:,1:3)./(Ys(:,1)+Ys(:,2)+Ys(:,3)))
subplot(2,1,2)
plot(t_obs,Ys(:,2)./(Ys(:,1)+Ys(:,2)+Ys(:,3)))

figure('Name','Time course for minum quiescence (absolute)','Position',[50 50 1200 600]);

[minQ,minQindex]=min(QuiescentCells);
[~, Ys] = ode15s(@(t,y) ode_system(t,y,final_particles(minQindex,:)), t_obs, y0, ode_opts);
subplot(2,1,1)
plot(t_obs,Ys(:,1:3))
subplot(2,1,2)
plot(t_obs,Ys(:,2))

writematrix(final_particles(minQindex,:),fullfile(results_dir,sprintf('MinQuiescSample%d.dat',sample_id)));




%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

% ---- ODE system (5 states, 9 parameters) ----------------------------------
%  CUSTOMISE THIS to your actual model equations.
%  y(1..5) = state variables;  p(1..9) = free parameters
function dydt = ode_system(t, y, p)
    dydt = zeros(5,1);
    
    % p(1)=gS p(2)=deltaS p(3)=muSR p(4)=muSP p(5)=muPS p(6)=muRP p(7)=muPR p(8)=dP p(9)=gR
    
    % Drug parameters fixed
    sigmaD = 0.9;deltaD = 0.00075;deltaI = 0.0025;
    
    dydt(1) = p(1)*y(1) - p(2)*y(5)*y(1) - (p(3) + p(4))*y(1) + p(5)*y(2);
    dydt(2) = p(4)*y(1) + p(6)*y(3) - (p(5)+p(7))*y(2) - p(8)*y(2);
    dydt(3) = p(9)*y(3) - p(6)*y(3) + p(3)*y(1) + p(7)*y(2);
    dydt(4) = sigmaD*heaviside(26-t) - deltaD*y(4)*(y(1)+y(2)+y(3));
    dydt(5) = deltaD*y(4)*y(1)-deltaI*y(5);
    
end

% ---- Simulate ODE and compute Euclidean distance -------------------------
function [d,TotalQ,RootMeanSqErr,volume] = simulate_and_distance(theta, t_obs, y0, Y_obs, opts)
   % try
        [t, Y_sim] = ode15s(@(t,y) ode_system(t,y,theta), t_obs, y0, opts);
        if size(Y_sim,1) ~= length(t_obs)
            d = Inf; return;
        end
        volume = Y_sim(:,1)+Y_sim(:,2)+Y_sim(:,3);
        time_span = max(t_obs,[],1) - min(t_obs,[],1);
        TotalQ = trapz(t_obs(:),(Y_sim(:,2)./volume(:)))/time_span;
        % Normalise each state by its data range to avoid scale dominance
        data_range = max(Y_obs,[],1) - min(Y_obs,[],1);
        data_range(data_range < 1e-12) = 1;
        diff = (volume - Y_obs) ./ data_range;
        d = sqrt(sum(diff(:).^2));
        RootMeanSqErr = rmse(volume,Y_obs);
   % catch
   %     d = Inf;
   %     TotalQ = NaN;
   % end
   %fprintf('distance = %d TotalQ = %f rmse = %f \n',d,TotalQ,RootMeanSqErr)
end

% ---- Sample from uniform prior -------------------------------------------
function theta = sample_prior(lb, ub, n)
    theta = lb + (ub - lb) .* rand(1, n);
end

% ---- Log-sum-exp (numerically stable) ------------------------------------
function s = logsumexp_vec(v)
    m = max(v);
    s = m + log(sum(exp(v - m)));
end
