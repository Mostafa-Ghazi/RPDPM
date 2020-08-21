function F = objective_markers(params, y, t, w, sigma, alpha, beta, loss_type, fit_type)

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

% The objective function for estimating the biomarker-specific parameters

[K, I, J] = size(y);
F = zeros(I * J * K, 1); % residuals
s = beta + alpha .* t; % linear DPS
for k = 1 : K
    F((k - 1) * I * J + 1 : k * I * J) = (y(k, :)' - sigmoid_function(fit_type, s(:), params(k), params(K + k), params(2 * K + k), params(3 * K + k), params(4 * K + k))) ./ sigma(k); % l2
end

if ~strcmp(loss_type, 'l2')
    switch loss_type
        case 'l1'
            G = sqrt(abs(F)); % l1
        case 'l1-l2'
            G = sqrt(2 * (sqrt(1 + F .^ 2) - 1)); % l1-l2
        case 'logistic'
            G = sqrt(1.205 ^ 2 * log(cosh(F / 1.205))); % logistic
        case 'welsch'
            G = sqrt(2.9846 ^ 2 * (1 - exp(- (F / 2.9846) .^ 2))); % welsch
        case 'fair'
            G = sqrt(1.3998 ^ 2 * (abs(F / 1.3998) - log(1 + abs(F / 1.3998)))); % fair
        case 'cauchy-lorentz'
            G = sqrt(2.3849 ^ 2 * log(1 + (F / 2.3849) .^ 2)); % cauchy-lorentz
        case 'talwar'
            G = sqrt(F .^ 2 .* (abs(F) <= 2.795) + 2.795 ^ 2 * (abs(F) > 2.795)); % talwar
        case 'huber'
            G = sqrt(F .^ 2 .* (abs(F) <= 1.345) + 1.345 ^ 2 * (2 * abs(F / 1.345) - 1) .* (abs(F) > 1.345)); % huber
        case 'andrew'
            G = sqrt(1.339 ^ 2 * ((1 - cos(F / 1.339)) .* (abs(F) <= (pi * 1.339)) + 2 * (abs(F) > (pi * 1.339)))); % andrew
        case 'tukey'
            G = sqrt(4.6851 ^ 2 / 3 * ((1 - (1 - (F / 4.6851) .^ 2) .^ 3) .* (abs(F) <= 4.6851) + 1 * (abs(F) > 4.6851))); % tukey
        case 'modified-huber'
            G = sqrt(1.2107 ^ 2 * ((1 - cos(abs(F / 1.2107))) .* (abs(F) <= (1.2107 * pi / 2)) + (abs(F / 1.2107) + 1 - pi / 2) .* (abs(F) > (1.2107 * pi / 2)))); % modified huber
    end
    F(~(isnan(G) | isinf(G))) = G(~(isnan(G) | isinf(G)));
end

F = sqrt(w) .* F; % weighted residuals
F = F(~(isnan(F) | isinf(F)));

end