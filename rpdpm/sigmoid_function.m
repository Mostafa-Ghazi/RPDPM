function f = sigmoid_function(fit_type, s, a, b, c, d, g)

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

% The logistic sigmoid function for fitting biomarker trajectories

% fit_type: fitting function type
% [a, d]: minimum or maximum (asymptotes) of the fit depending on sign(b)
% c: inflection point
% g: symmetry
% s: disease progression score (DPS)
% f: sigmoidal function for fitting trajectory of biomarker

f = (a - d) .* (1 + exp(- b .* (s - c) ./ g) ./ g) .^ (- g) + d; % Proposed logistic function (modified Stannard)
% (a - d) * b * (1 + 1 / g) ^ (- g - 1) / g: slope at c

if ~strcmp(fit_type, 'proposed')
    switch fit_type
        case 'richards'
            f = (a - d) .* (1 + g .* exp(- b .* (s - c))) .^ (- 1 ./ g) + d; % Richards logistic function
            % (a - d) * b * (1 + g) ^ (- 1 / g - 1): slope at c
        case 'verhulst'
            f = (a - d) ./ (1 + exp(- b .* (s - c))) + d; % Verhulst logistic function
            % (a - d) * b / 4: slope at c
        case 'gompertz'
            f = (a - d) .* exp(- exp(- b .* (s - c))) + d; % Gompertz logistic function
            % (a - d) * b / e: slope at c
    end
end

end