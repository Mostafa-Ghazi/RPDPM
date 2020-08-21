function [alpha, beta, sigma, a, b, c, d, g, loss, bic] = rpdpm_train(y, t, labels, classes, ranges, loss_type, fit_type, iters, optims)

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

% Training the robust parametric disease progression model

% I: number of subjects
% K: number of biomarkers
% J: number of visits (time points)
% y: biomarker measurements [K x I x J]
% t: visiting age of subjects [I x J]
% labels: array of labels [I x J]
% classes: distinct class labels ordered w.r.t. disease progression
% ranges: actual ranges of biomarkers [K x 2]
% loss_type: robust estimation loss function type
% fit_type: fitting function type
% iters: number of alternating itersations
% optims: optimization options for fitting biomarkers
% [alpha, beta]: subject-specific parameters each of which [I x (iters + 1)]
% sigma: scale factor of biomarker errors [K x 1]
% [a, b, c, d, g]: biomarker-specific parameters each of which [K x (iters + 1)]
% loss: training loss [iters x 1]
% bic: Bayesian information critersion [iters x 1]

% Initializing the model parameters
[K, I, J] = size(y);
sigma = ones(K, 1);
alpha = zeros(I, iters + 1);
beta = zeros(I, iters + 1);
a = ones(K, iters + 1);
b = ones(K, iters + 1);
c = zeros(K, iters + 1);
d = zeros(K, iters + 1);
g = ones(K, iters + 1);

% Number of optimization parameters
if strcmpi(fit_type, 'proposed') || strcmpi(fit_type, 'richards')
    num_params = 2 * I + 5 * K - 2 * sum(~sum(isinf(ranges), 2));
elseif strcmpi(fit_type, 'verhulst') || strcmpi(fit_type, 'gompertz')
    num_params = 2 * I + 4 * K - 2 * sum(~sum(isinf(ranges), 2));
end

% Updating initial values of the parameters
for k = 1 : K
    sigma(k) = nanstd(y(k, cellfun(@(x) strcmp(x, classes(1)), labels(:)))); % standard deviation as a dispersion measure
    if nanmean(y(k, cellfun(@(x) strcmp(x, classes(1)), labels(:)))) < nanmean(y(k, cellfun(@(x) strcmp(x, classes(end)), labels(:)))) % positive slope (ascending)
        d(k, 1) = nanmin(y(k, :)); % a > d
        a(k, 1) = nanmax(y(k, :)); % a > d
        b(k, 1) = 4 / (a(k, 1) - d(k, 1)); % b > 0
    else % negative slope (descending)
        d(k, 1) = nanmax(y(k, :)); % a < d
        a(k, 1) = nanmin(y(k, :)); % a < d
        b(k, 1) = - 4 / (a(k, 1) - d(k, 1)); % b > 0
    end
end

% Assigning rangess of the model parameters
lb = [- Inf(K, 1); zeros(K, 1); zeros(K, 1); - Inf(K, 1); 0.1 * ones(K, 1); zeros(I, 1); - Inf(I, 1)]; % lower bounds of parameters [a, b, c, d, g, alpha, beta]
ub = [Inf(K, 1); Inf(K, 1); Inf(K, 1); Inf(K, 1); 100 * ones(K, 1); Inf(I, 1); Inf(I, 1)]; % upper bounds of parameters [a, b, c, d, g, alpha, beta]
for k = 1 : K
    if ~sum(isinf(ranges(k, :)), 2) % cognitive tests with known ranges
        if a(k, 1) > d(k, 1) % positive slope (ascending)
            lb(k) = ranges(k, 2); % a
            ub(k) = ranges(k, 2); % a
            lb(3 * K + k) = ranges(k, 1); % d
            ub(3 * K + k) = ranges(k, 1); % d
        else % negative slope (descending)
            lb(k) = ranges(k, 1); % a
            ub(k) = ranges(k, 1); % a
            lb(3 * K + k) = ranges(k, 2); % d
            ub(3 * K + k) = ranges(k, 2); % d
        end
    else % other biomarkers
        lb(k) = ranges(k, 1); % a
        ub(k) = ranges(k, 2); % a
        lb(3 * K + k) = ranges(k, 1); % d
        ub(3 * K + k) = ranges(k, 2); % d
    end
end

% Weights for normalizing the residuals
w = zeros(I, J, K);
for i = 1 : I
    w(i, :) = 1 / sum(sum(~isnan(y(:, i, :)), 3), 1); % normalized w.r.t. the number of available points per subject
end
w_i = permute(w, [1, 3, 2]);

% Fitting the model
fprintf('iteration #'); % display itersations
loss = zeros(iters, 1); % training loss
bic = zeros(iters, 1); % Bayesian information critersion
for l = 1 : iters
    
    fprintf(' %i', l); % displays current itersation's number
    
    % Replicating subject-specific parameters for all time-points
    alpha_k = repmat(alpha(:, l), 1, J);
    beta_k = repmat(beta(:, l), 1, J);
    
    % Learning biomarker-specific parameters by fitting sigmoids to biomarkers
    F = @(x) objective_markers(x, y, t, w(:), sigma, alpha_k, beta_k, loss_type, fit_type);
    x0 = [a(:, l); b(:, l); c(:, l); d(:, l); g(:, l)];
    params = lsqnonlin(F, x0, lb(1 : end - 2 * I), ub(1 : end - 2 * I), optims);
    
    % Updating biomarker-specific parameters
    a(:, l + 1) = params(1 : K);
    b(:, l + 1) = params(K + 1 : 2 * K);
    c(:, l + 1) = params(2 * K + 1 : 3 * K);
    d(:, l + 1) = params(3 * K + 1 : 4 * K);
    g(:, l + 1) = params(4 * K + 1 : 5 * K);
    
    % Replicating biomarker-specific parameters for all time-points
    a_i = repmat(a(:, l + 1), 1, J);
    b_i = repmat(b(:, l + 1), 1, J);
    c_i = repmat(c(:, l + 1), 1, J);
    d_i = repmat(d(:, l + 1), 1, J);
    g_i = repmat(g(:, l + 1), 1, J);
    sigma_i = repmat(sigma, 1, J);
    
    % Learning subject-specific parameters
    for i = 1 : I
        
        % Finding measures with available values
        y_i = squeeze(y(:, i, :));
        t_i = repmat(t(i, :), K, 1);
        idx = find(~isnan(t_i(:)) & ~isnan(y_i(:)));
        
        % Fitting sigmoids to biomarkers
        F = @(x) objective_subject(x, y_i(idx), t_i(idx), w_i(i, idx)', sigma_i(idx), a_i(idx), b_i(idx), c_i(idx), d_i(idx), g_i(idx), loss_type, fit_type);
        x0 = [alpha(i, l); beta(i, l)];
        params = lsqnonlin(F, x0, lb(end - 2 * I + i : I : end), ub(end - 2 * I + i : I : end), optims);
        
        % Updating subject-specific parameters
        alpha(i, l + 1) = params(1);
        beta(i, l + 1) = params(2);
        
    end
    
    % Training set performance and goodness of fit
    theta = [a(:, l + 1); b(:, l + 1); c(:, l + 1); d(:, l + 1); g(:, l + 1)];
    res = objective_markers(theta, y, t, w(:), sigma, repmat(alpha(:, l + 1), 1, J), repmat(beta(:, l + 1), 1, J), loss_type, fit_type);
    loss(l) = sum(res .^ 2);
    bic(l) = model_selection(loss(l), numel(res), num_params);
    
end
fprintf('\n'); % go to the next line

end