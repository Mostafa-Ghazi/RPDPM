function [posterior, likelihood, prior, evidence] = bayes_train(x, labels, unique_labels, models)

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

% Training Bayesian classifiers using kernel density estimation

% x: vector of distributions of x-values in data scatter
% labels: vector of true labels for a multiclass problem
% unique_labels: array of unique class labels ordered w.r.t. x-values progression
% models: array of distribution names to fit likelihoods ordered w.r.t. x-values progression
% posterior: array of functions of posterior probabilities
% likelihood: array of functions of class-conditional probabilities
% prior: vector of prior probabilities of classes
% evidence: function of prior probabilities of data

% Finding indices of available points
idx = find(~cellfun(@isempty, labels(:)) & ~isnan(x(:)));
labels = labels(idx);
x = x(idx);

% Finding prior probabilities of classes and fitting models to likelihood distributions
C = numel(unique_labels); % number of unique class labels
prior = zeros(C, 1);
likelihood = cell(C, 1);
for c = 1 : C
    prior(c) = sum(strcmp(labels, unique_labels(c))) / numel(idx);
    theta = fitdist(x(strcmp(labels, unique_labels(c))), models{c});
    likelihood{c} = @(x) pdf(theta, x);
end

% Finding prior probabilities of data
evidence = @(x) 0;
for c = 1 : C
    evidence = @(x) (evidence(x) + prior(c) * likelihood{c}(x));
end

% Finding posterior probabilities
posterior = cell(C, 1);
for c = 1 : C
    posterior{c} = @(x) (prior(c) * likelihood{c}(x) ./ evidence(x));
end

end