
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                 Source Code Author Information                  %%%%%
%%%%%                     Mostafa Mehdipour Ghazi                     %%%%%
%%%%%                   mostafa.mehdipour@gmail.com                   %%%%%
%%%%%                      Created on 01/08/2017                      %%%%%
%%%%%                      Updated on 01/08/2020                      %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  This code is an implementation of the algorithm published in:  %%%%%
%%%%%  Robust parametric modeling of Alzheimer's disease progression  %%%%%
%%%%%                https://arxiv.org/abs/1908.05338                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Training, validation, and testing the proposed model for biomarker value prediction and clinical status classification using a simulated data 

%% Data Preparation

tic

restoredefaultpath
close all
clear
clc

% Reading the longitudinal data
data = readtable('./data/simulation_data.csv');
addpath('./rpdpm'); % add path to the source codes

% Extracting desired fields from the longitudinal data
subjects_vec = data{:, 'SubjectID'}; % subject IDs
labels_vec = data{:, 'Label'}; % visiting status of subjects
ages_vec = data{:, 'Age'}; % visiting ages of subjects
biomarkers = setdiff(data.Properties.VariableNames, {'SubjectID', 'Label', 'Age'}); % biomarker names
y_vec = data{:, ismember(data.Properties.VariableNames, biomarkers)}; % biomarker values

% Removing data lacking visiting age information and vice versa
y_vec(isnan(ages_vec), :) = NaN;
ages_vec(~sum(~isnan(y_vec), 2)) = [];
labels_vec(~sum(~isnan(y_vec), 2)) = [];
y_vec(~sum(~isnan(y_vec), 2), :) = [];

% Arranging data in a desired form of 2D/3D arrays
subjects = unique(subjects_vec, 'sorted'); % unique subjects
I = numel(subjects); % number of unique subjects
K = numel(biomarkers); % number of biomarkers
J = histcounts(subjects_vec, I); % number of available visits per subject
y = NaN(K, I, max(J)); % biomarker values
ages = NaN(I, max(J)); % visiting ages of subjects
labels = cell(I, max(J)); % visiting status of subjects
for i = 1 : I
    idx_i = find(subjects_vec == subjects(i));
    [val_sort, idx_sort] = sort(ages_vec(idx_i), 'ascend');
    ages(i, 1 : numel(idx_i)) = val_sort;
    y(:, i, 1 : numel(idx_i)) = y_vec(idx_i(idx_sort), :)';
    labels(i, 1 : numel(idx_i)) = labels_vec(idx_i(idx_sort));
end

% Finding clinical status of subjects in their first and last available visits
idx_labels = ~cellfun(@isempty, labels); % logical indices of available labels
labels_exist = labels(sum(idx_labels, 2) > 0, :); % labels having at least one label per subject
idx_labels_exist = idx_labels(sum(idx_labels, 2) > 0, :); % logical indices of available labels having at least one label per subject
idx_labels_first = arrayfun(@(x) find(idx_labels_exist(x, :), 1, 'first'), 1 : size(idx_labels_exist, 1)); % index of first label per subject
idx_labels_first = sub2ind(size(idx_labels_exist), 1 : size(idx_labels_exist, 1), idx_labels_first); % linear array indices of first available labels
idx_labels_last = arrayfun(@(x) find(idx_labels_exist(x, :), 1, 'last'), 1 : size(idx_labels_exist, 1)); % index of last label per subject
idx_labels_last = sub2ind(size(idx_labels_exist), 1 : size(idx_labels_exist, 1), idx_labels_last); % linear array indices of last available labels
diagnose = cell(2, I); % first and last diagnosis of subjects
diagnose(:, sum(idx_labels, 2) > 0) = [labels_exist(idx_labels_first); labels_exist(idx_labels_last)];

% Filtering data by rejecting outliers
ranges = [0 30; 0 Inf; -Inf Inf]; % ranges of biomarkers
for k = 1 : K
    idx_ij = find(~isnan(y(k, :)));
    y_k = y(k, idx_ij);
    idx_out = find(y_k < ranges(k, 1) | y_k > ranges(k, 2));
    y(k, idx_ij(idx_out)) = NaN; % reject outliers
end

% Removing subjects with a few number of distinct visits
idx_rmv = find(sum(~isnan(nanmean(y, 1)), 3) < 2 | sum(~cellfun(@isempty, labels), 2)' < 2);
y(:, idx_rmv, :) = [];
ages(idx_rmv, :) = [];
labels(idx_rmv, :) = [];
diagnose(:, idx_rmv) = [];
subjects(idx_rmv, :) = [];

% Splitting data to training and test subsets
ratio_test = 0.2; % proportion of test subjects to entire data
classes = {'CN', 'MCI', 'AD'}; % distinct class labels ordered w.r.t. disease progression
num_points = sum(sum(~isnan(y), 3), 1); % number of available points per subject
idx_test = []; % indices of test subjects
rng(0); % random number generation seed
for u = 1 : numel(classes)
    for uu = u : numel(classes)
        idx_u = find(strcmp(diagnose(1, :), classes(u)) & strcmp(diagnose(2, :), classes(uu))); % first and last diagnoses
        med_u = median(num_points(idx_u)); % within class median of available points
        idx_many = find(num_points(idx_u) > med_u); % subjects with many points
        if ratio_test * numel(idx_many) > 1
            idx_rnd = randperm(numel(idx_many), floor(ratio_test * numel(idx_many)));
            idx_test = cat(2, idx_test, idx_u(idx_many(idx_rnd)));
        end
        idx_few = find(num_points(idx_u) <= med_u); % subjects with few points
        if ratio_test * numel(idx_few) > 1
            idx_rnd = randperm(numel(idx_few), floor(ratio_test * numel(idx_few)));
            idx_test = cat(2, idx_test, idx_u(idx_few(idx_rnd)));
        end
    end
end
idx_train = find(~ismember(1 : I, idx_test)); % indices of training subjects
y_train = y(:, idx_train, :);
ages_train = ages(idx_train, :);
labels_train = labels(idx_train, :);
diagnose_train = diagnose(:, idx_train);
subjects_train = subjects(idx_train, :);
y_test = y(:, idx_test, :);
ages_test = ages(idx_test, :);
labels_test = labels(idx_test, :);
diagnose_test = diagnose(:, idx_test);
subjects_test = subjects(idx_test, :);

time_data = toc;

%% Robust Parametric DPM Training

tic

% Number of biomarkers, subjects, and visits
[K1, I1, J1] = size(y_train);

% Optimization options
optim_options = optimoptions('lsqnonlin', 'Display', 'none', 'Jacobian', 'off', 'DerivativeCheck', 'off', 'MaxIter', 1e+3, 'TolFun', 1e-6, 'TolX', 1e-6, 'Algorithm', 'trust-region-reflective');
loss_type = 'logistic'; % robust estimation loss function type
fit_type = 'proposed'; % fitting function type
L_max = 50; % maximum number of alternating iterations
L_min = 10; % minimum number of alternating iterations
btstrp = 1; % number of bootstraps

% Parameter initialization for different bootstraps
Idx_train = zeros(I1, btstrp);
Alpha_train = zeros(I1, btstrp);
Beta_train = zeros(I1, btstrp);
Sigma = zeros(K1, btstrp);
A = zeros(K1, btstrp);
B = zeros(K1, btstrp);
C = zeros(K1, btstrp);
D = zeros(K1, btstrp);
G = zeros(K1, btstrp);
BIC = zeros(L_max, btstrp);
Iters_opt = zeros(1, btstrp);
DPS_train = zeros(I1, J1, btstrp);
Loss_train = zeros(L_max, btstrp);
Loss_valid = zeros(L_max, btstrp);
NMAE_valid = zeros(L_max, btstrp);

% Biomarker modeling
num_points = sum(sum(~isnan(y_train), 3), 1); % number of available points per training subject
for n = 1 : btstrp
    
    fprintf('bootstrap %i> ', n); % display bootstrap number
    
    % Randomly sampling a group of within-class subjects with replacement for training and the rest for validation
    idx_btstrp = []; % indices of bootstrapped subjects
    rng(n); % random number generation seed
    for u = 1 : numel(classes)
        for uu = u : numel(classes)
            idx_u = find(strcmp(diagnose_train(1, :), classes(u)) & strcmp(diagnose_train(2, :), classes(uu))); % first and last diagnoses
            med_u = median(num_points(idx_u)); % within class median of available points
            idx_many = find(num_points(idx_u) > med_u); % subjects with many points
            if numel(idx_many)
                idx_rnd = randi(numel(idx_many), numel(idx_many), 1);
                idx_btstrp = cat(2, idx_btstrp, idx_u(idx_many(idx_rnd)));
            end
            idx_few = find(num_points(idx_u) <= med_u); % subjects with few points
            if numel(idx_few)
                idx_rnd = randi(numel(idx_few), numel(idx_few), 1);
                idx_btstrp = cat(2, idx_btstrp, idx_u(idx_few(idx_rnd)));
            end
        end
    end
    Idx_train(:, n) = idx_btstrp; % indices of training subjects
    idx_valid = find(~ismember(1 : I1, idx_btstrp)); % indices of validation subjects
    
    % Fitting sigmoids to biomarkers
    [idx_train, idx_btstrp1, idx_btstrp2] = unique(idx_btstrp); % unique training indices
    [alpha, beta, Sigma(:, n), a, b, c, d, g, Loss_train(:, n), BIC(:, n)] = ...
        rpdpm_train(y_train(:, idx_train, :), ages_train(idx_train, :), labels_train(idx_train, :), classes, ranges, loss_type, fit_type, L_max, optim_options);
    
    % Validation set modeling performance
    w = zeros(numel(idx_valid), J1, K1);
    for i = 1 : numel(idx_valid)
        w(i, :) = 1 / sum(sum(~isnan(y_train(:, idx_valid(i), :)), 3), 1); % normalized w.r.t. the number of available points per subject
    end
    alpha_valid = zeros(numel(idx_valid), 1);
    beta_valid = zeros(numel(idx_valid), 1);
    for l = 1 : L_max
        theta = [a(:, l + 1), b(:, l + 1), c(:, l + 1), d(:, l + 1), g(:, l + 1)];
        [alpha_valid, beta_valid, ~, y_valid_fit] = rpdpm_test(y_train(:, idx_valid, :), ages_train(idx_valid, :), ...
            alpha_valid, beta_valid, Sigma(:, n), theta(:, 1), theta(:, 2), theta(:, 3), theta(:, 4), theta(:, 5), loss_type, fit_type, optim_options);
        Loss_valid(l, n) = sum(objective_markers(theta(:), y_train(:, idx_valid, :), ages_train(idx_valid, :), w(:), Sigma(:, n), ...
            repmat(alpha_valid, 1, J1), repmat(beta_valid, 1, J1), loss_type, fit_type) .^ 2);
        ae_valid = abs(y_train(:, idx_valid, :) - y_valid_fit); % absolute errors
        ae_valid(1 : K1, :) = ae_valid(1 : K1, :) ./ Sigma(1 : K1, n); % normalized (scaled) absolute errors
        NMAE_valid(l, n) = nanmean(ae_valid(:)); % normalized mean absolute errors
    end
    
    % Optimal parameters
    [~, L_opt] = min(flip(Loss_valid(L_min : end, n))); % the latest minimum
    Iters_opt(n) = L_max - L_opt + 1;
    alpha = alpha(:, Iters_opt(n) + 1);
    beta = beta(:, Iters_opt(n) + 1);
    A(:, n) = a(:, Iters_opt(n) + 1);
    B(:, n) = b(:, Iters_opt(n) + 1);
    C(:, n) = c(:, Iters_opt(n) + 1);
    D(:, n) = d(:, Iters_opt(n) + 1);
    G(:, n) = g(:, Iters_opt(n) + 1);
    dps = repmat(beta, 1, J1) + repmat(alpha, 1, J1) .* ages_train(idx_train, :); % linear DPS
    
    % Standardizing DPS and related parameters w.r.t. the normal group
    dps_cn = dps(cellfun(@(x) strcmp(x, classes{1}), labels_train(idx_train, :))); % normal visits
    mu_cn = nanmean(dps_cn(:)); % mean
    sigma_cn = nanstd(dps_cn(:)); % standard deviation
    alpha = alpha / sigma_cn; % normalized alpha
    beta = (beta - mu_cn) / sigma_cn; % normalized beta
    dps = (dps - mu_cn) / sigma_cn; % normalized DPS
    C(:, n) = (C(:, n) - mu_cn) / sigma_cn; % normalized c
    B(:, n) = B(:, n) * sigma_cn; % normalized b
    
    % Replicated parameter values for all bootstrapped subjects
    Alpha_train(:, n) = alpha(idx_btstrp2);
    Beta_train(:, n) = beta(idx_btstrp2);
    DPS_train(:, :, n) = dps(idx_btstrp2, :);
    
end
bic_avg = mean(BIC(logical(full(sparse(Iters_opt, 1 : btstrp, ones(1, btstrp)))))); % average BIC
bic_std = std(BIC(logical(full(sparse(Iters_opt, 1 : btstrp, ones(1, btstrp)))))); % standard deviation of BICs
nmae_valid_avg = mean(NMAE_valid(logical(full(sparse(Iters_opt, 1 : btstrp, ones(1, btstrp)))))); % average of normalized mean absolute errors
nmae_valid_std = std(NMAE_valid(logical(full(sparse(Iters_opt, 1 : btstrp, ones(1, btstrp)))))); % standard deviation of normalized mean absolute errors
fprintf('Training BIC = %4.4f \x00B1 %4.4f \n', [bic_avg, bic_std]); % display training modeling performance
fprintf('Validation NMAE = %4.4f \x00B1 %4.4f \n', [nmae_valid_avg, nmae_valid_std]); % display validation modeling performance

% Visualizing the DPM results
plot_marker_order(C, biomarkers);
plot_dps_pdf(DPS_train, C, labels_train, classes, biomarkers, Idx_train);
plot_trajectory(A, B, C, D, G, DPS_train, y_train, labels_train, classes, biomarkers, Idx_train, fit_type);

time_rpdpm_train = toc;

%% Robust Parametric DPM Testing

tic

% Estimating test subject parameters using biomaker values and ages
[K2, I2, J2] = size(y_test); % number of biomarkers, subjects, and visits
y_test_fit = zeros(K2, I2, J2, btstrp); % estimated values of biomarkers
DPS_test = zeros(I2, J2, btstrp);
Alpha_test = zeros(I2, btstrp);
Beta_test = zeros(I2, btstrp);
for n = 1 : btstrp
    [Alpha_test(:, n), Beta_test(:, n), DPS_test(:, :, n), y_test_fit(:, :, :, n)] = rpdpm_test(y_test, ages_test, ...
        Alpha_test(:, n), Beta_test(:, n), Sigma(:, n), A(:, n), B(:, n), C(:, n), D(:, n), G(:, n), loss_type, fit_type, optim_options);
end

% Test set modeling performance
NMAE_test = zeros(btstrp, 1); % overall test error per bootstrap
MAE_test = zeros(K2, btstrp); % test error per biomarker per bootstrap
for n = 1 : btstrp
    ae_test = abs(y_test - y_test_fit(:, :, :, n)); % absolute errors
    MAE_test(1 : K2, n) = nanmean(ae_test(1 : K2, :), 2); % mean absolute error
    ae_test = ae_test(1 : K2, :) ./ Sigma(1 : K2, n); % normalized (scaled) absolute errors
    NMAE_test(n) = nanmean(ae_test(:)); % normalized mean absolute error
end
nmae_test_avg = mean(NMAE_test); % average of normalized mean absolute errors
nmae_test_std = std(NMAE_test); % standard deviation of normalized mean absolute errors
mae_test_avg = mean(MAE_test, 2); % average of mean absolute errors
mae_test_std = std(MAE_test, 0, 2); % standard deviation of mean absolute errors
fprintf('Test NMAE = %4.4f \x00B1 %4.4f \n', [nmae_test_avg, nmae_test_std]); % display test modeling performance

time_rpdpm_test = toc;

%% Robust Parametric DPM Classification Training

tic

% Training Bayesian classifiers per bootstrap
UC = numel(classes); % number of unique class labels
dist_models = repmat({'kernel'}, 1, UC); % distribution models for likelihoods
posterior_bayes = cell(UC, btstrp);
likelihood_bayes = cell(UC, btstrp);
evidence_bayes = cell(1, btstrp);
prior_bayes = zeros(UC, btstrp);
for n = 1 : btstrp
    [posterior_bayes(:, n), likelihood_bayes(:, n), prior_bayes(:, n), evidence_bayes{n}] = bayes_train(DPS_train(:, :, n), labels_train(Idx_train(:, n), :), classes, dist_models);
end

time_classify_train = toc;

%% Robust Parametric DPM Classification Testing

tic

% Classification performance analysis for test samples
predict_scores_test = nan(UC, I2 * J2, btstrp); % posterior probabilities for test scores
auc_test = zeros(1, btstrp); % test multiclass area under the curve per bootstrap
for n = 1 : btstrp
    dps_test_n = DPS_test(:, :, n);
    for cc = 1 : UC
        predict_scores_test(cc, :, n) = posterior_bayes{cc, n}(dps_test_n(:)');
    end
    auc_test(n) = multiclass_auc(predict_scores_test(:, :, n)', labels_test(:), classes);
end
auc_test_avg = mean(auc_test); % average multiclass area under the curve
auc_test_std = std(auc_test); % standard deviation of multiclass areas under the curves
fprintf('Test AUC = %4.4f \x00B1 %4.4f \n', [auc_test_avg, auc_test_std]); % display test classification performance

time_classify_test = toc;
