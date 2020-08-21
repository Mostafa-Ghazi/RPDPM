function plot_dps_pdf(S, C, labels, classes, biomarkers, idx)

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

% Plotting within-class distributions of DPSs and ordered biomarkers

% I: number of subjects
% J: number of visits (time points)
% btstrp: number of bootstraps
% S: disease progression scores (DPS) [I x J x btstrp]
% C: inflection points [K x btstrp]
% labels: array of labels [I x J]
% classes: distinct class labels ordered w.r.t. disease progression
% biomarkers: array of biomarker names [1 x K]
% idx: indices of randomly sampled subjects [I x btstrp]

% Calculating average dps per subject visit across different bootstraps
[I, J, btstrp] = size(S); % number of subjects, visits, and bootstraps
S = permute(S, [2, 1, 3]);
S_avg = zeros(J, I);
for i = 1 : I
    S_avg(:, i) = nanmean(S(:, find(idx == i)), 2);
end
s = S_avg';

% DPS range
% s_range = [min(s(:)) - 0.05 * (max(s(:)) - min(s(:))), max(s(:)) + 0.05 * (max(s(:)) - min(s(:)))];
s_range = [-1, 2];

% Marker colors
CM = {'k'; 'b'; 'g'; 'r'; 'm'; 'c'; 'y'};

% Plotting within-class DPS distributions
h = figure('Visible', 'off'); clf;
h1 = subplot(2, 1, 1); % displays pdf
hold on;
for u = 1 : numel(classes)
    s_u = s(strcmp(labels, classes(u)));
    [pd, s_u] = ksdensity(s_u);
    h_pdf(u) = area(s_u, pd, 'FaceColor', CM{u + 1}, 'LineStyle', 'none', 'FaceAlpha', 0.45);
end
idx_nu = cellfun(@isempty, labels);
if sum(idx_nu(:))
    s_nu = s(idx_nu);
    [pd, s_nu] = ksdensity(s_nu);
    h_pdf(u + 1) = area(s_nu, pd, 'FaceColor', CM{1}, 'LineStyle', 'none', 'FaceAlpha', 0.45);
    classes = cat(2, classes, 'Missing');
end
hold off;
xlim(s_range), xlabel('DPS'), ylabel('Probability'), set(gca, 'Ticklength', [0.006 0], 'XTickLabel', [], 'FontName', 'times', 'FontSize', 8);
h1_leg = legend(h_pdf, classes, 'Location', 'northeast', 'FontSize', 9, 'TextColor', 'black'); title(h1_leg, 'Diagnosis'); legend boxoff;

% Plotting ordered biomarker variations
h2 = subplot(2, 1, 2);
if size(C, 2) == 1
    C = repmat(C, 1, 2);
end
[~, idx_sort] = sort(median(C, 2), 'descend'); % sort markers based on median inflection points
boxplot(C(idx_sort, :)', 'orientation', 'horizontal', 'PlotStyle', 'traditional', 'BoxStyle', 'outline', 'MedianStyle', 'line', 'Widths', 0.4, 'Symbol', '');
xlim(s_range), xlabel('DPS'), set(gca, 'Ticklength', [0.006 0], 'YTickLabel', biomarkers(idx_sort), 'box', 'off', 'FontName', 'times', 'FontSize', 7.5);
set(findobj(gcf, 'tag', 'Median'), 'Color', 'k'); set(findobj(gcf, 'tag', 'Box'), 'Color', 'm'); % change box plot color

% Stacking the two plots, aligning them, and saving the figure
p1 = get(h1, 'Position'); p2 = get(h2, 'Position'); p1(2) = p2(2) + p2(4) + 0.012; p1(1) = p2(1); set(h1, 'pos', p1);
print(h, '-dpng', 'fig_dps_pdf.png'); close; % save the figure

end