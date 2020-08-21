function plot_trajectory(A, B, C, D, G, S, y, labels, classes, biomarkers, idx, fit_type)

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

% Plotting trajectories (biomarker values versus disease progression scores)

% I: number of subjects
% K: number of biomarkers
% J: number of visits (time points)
% btstrp: number of bootstraps
% [A, B, C, D, G]: biomarker-specific parameters each of which [K x btstrp]
% S: disease progression scores (DPS) [I x J x btstrp]
% y: biomarker measurements [K x I x J]
% labels: array of labels [I x J]
% classes: distinct class labels ordered w.r.t. disease progression
% biomarkers: array of biomarker names [1 x K]
% idx: indices of randomly sampled subjects [I x btstrp]
% fit_type: fitting function type

[K, I, J] = size(y); % number of biomarkers, subjects, and visits
btstrp = size(S, 3); % number of bootstraps
classes_extend = classes; % appended class labels with missing labels

% Calculating average dps per subject visit across different bootstraps
S = permute(S, [2, 1, 3]);
S_avg = zeros(J, I);
for i = 1 : I
    S_avg(:, i) = nanmean(S(:, find(idx == i)), 2);
end
S = S_avg';

% DPS resolution
% s = [min(S(:)) - 0.05 * (max(S(:)) - min(S(:))) : 0.001 :  max(S(:)) + 0.05 * (max(S(:)) - min(S(:)))];
s = [-1 : 0.001 : 2];

CM = {'k'; 'b'; 'g'; 'r'; 'm'; 'c'; 'y'}; % marker and line colors
M = {'.'; 'o'; 'x'; '+'; '*'; 's'; 'p'}; % markers
LS = {'-'; ':'; '-.'; '--'}; % line styles

% Plotting biomarker trajectories individually
for k = 1 : K
    
    h = figure('Visible', 'off'); clf;
    hold on;
    
    % Plotting data scatter
    for u = 1 : numel(classes)
        idx_u = strcmp(labels, classes(u));
        y_k_u = y(k, idx_u);
        S_u = S(idx_u);
        h_data(u) = plot(S_u(:), y_k_u(:), M{1}, 'Color', CM{u + 1}, 'MarkerSize', 8); % plot labeled biomarker data
    end
    idx_nu = cellfun(@isempty, labels);
    if sum(idx_nu(:))
        y_k_nu = y(k, idx_nu);
        S_nu = S(idx_nu);
        h_data(u + 1) = plot(S_nu(:), y_k_nu(:), M{1}, 'Color', [0.5 0.5 0.5], 'MarkerSize', 8); % plot unlabeled biomarker data
        classes_extend = cat(2, classes, 'Missing');
    end
    
    % Plotting logistic curves
    f_avg = 0; % average function
    for n = 1 : btstrp
        f = sigmoid_function(fit_type, s, A(k, n), B(k, n), C(k, n), D(k, n), G(k, n)); % fit logistic curve
        plot(s, f, 'Color', [0.7 0.7 0.7], 'LineStyle', LS{1}, 'LineWidth', 0.65);
        f_avg = f_avg + f / btstrp;
    end
    plot(s, f_avg, 'Color', CM{1}, 'LineStyle', LS{1}, 'LineWidth', 1.5);
    hold off;
    
    % Saving the figure
    ylim([min(y(k, :)) - 0.05 * (max(y(k, :)) - min(y(k, :))),  max(y(k, :)) + 0.05 * (max(y(k, :)) - min(y(k, :)))]);
    xlim([s(1), s(end)]), axis square, title({biomarkers{k}, ''}), xlabel('DPS'), ylabel('Biomarker value');
    set(gca, 'FontName', 'times', 'FontSize', 13, 'TickLength', [0.01, 0.01]);
    h_leg = legend(h_data, classes_extend, 'Location', 'northeastoutside', 'Orientation', 'vertical', 'FontSize', 14, 'TextColor', 'black'); legend boxoff;
    title(h_leg, 'Diagnosis');
    print(h, '-dpng', ['fig_trajectory' num2str(k) '.png']); close; % save the figure
    
end

% Plotting all normalized biomarker trajectories in the one figure
h = figure('Visible', 'off'); clf;
hold on;
for k = 1 : K
    f_avg = 0; % average function
    for n = 1 : btstrp
        f = sigmoid_function(fit_type, s, A(k, n), B(k, n), C(k, n), D(k, n), G(k, n)); % fit logistic curve
        f = (f - D(k, n)) / (A(k, n) - D(k, n)); % scale to [0 1]
        f_avg = f_avg + f / btstrp;
    end
    cm = circshift(CM, - k); % circular shifting of colors
    h_plot(k) = plot(s, f_avg, 'Color', cm{1}, 'LineStyle', LS{ceil(k/numel(cm))}, 'LineWidth', 1.2);
end
hold off;
ylim([0 1]), axis square, xlabel('DPS'), ylabel('Normalized biomarker value');
set(gca, 'FontName', 'times', 'FontSize', 13, 'TickLength', [0.01, 0.01]);
legend(h_plot, biomarkers, 'Location', 'northeastoutside', 'FontSize', 11, 'TextColor', 'black'); legend boxoff;
print(h, '-dpng', 'fig_trajectories.png'); close; % save the figure

end