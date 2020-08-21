function [bicr, aicr] = model_selection(loss, n, p)

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

% Model selection in robust regression using AIC and BIC

% loss: sum of residual errors (residual sum of squares in l2 case)
% n: number of observations
% p: number of model parameters
% bicr: Bayesian information criterion for robust estimators
% aicr: Akaike information criterion for robust estimators

% Calculating robust aic and bic
% bic = n * log(loss / n) + log(n) * p; % l2
% aic = n * log(loss / n) + 2 * p; % l2
bicr = 2 * loss + log(n) * p;
aicr = 2 * loss + 2 * p;
if n / p < 40 % second-order bias adjustment for small sample size
    aicr = aicr + (2 * p * (p + 1)) ./ (n - p - 1);
end

end