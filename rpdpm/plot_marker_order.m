function plot_marker_order(C, biomarker)

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

% Plotting temporal ordering of biomarkers in the disease course

% K: number of biomarkers
% btstrp: number of bootstraps
% C: inflection points [K x btstrp]
% biomarkers: array of biomarker names [1 x K]

[K, btstrp] = size(C); % number of biomarkers and bootstraps

% Ordering probabilites
[~, order] = sort(C, 1, 'ascend');
order_freq = zeros(numel(biomarker));
for k = 1 : K
    for n = 1 : K
        order_freq(k, n) = sum(order(n, :) == k);
    end
end
order_freq = order_freq / btstrp;

% Sorting biomarkers to obtain a descending diagonal-like array
order_freq_sum = cumsum(order_freq, 2, 'forward'); % cdf per biomarker
for i = 1 : K - 1
    for j = i + 1 : K
        order_freq_sum_diff = order_freq_sum(j, :) - order_freq_sum(i, :); % difference between consecutive cdfs
        idx_order_freq_sum_diff_first = find(order_freq_sum_diff ~= 0, 1, 'first'); % first difference between cdfs
        if order_freq_sum_diff(idx_order_freq_sum_diff_first) > 0
            order_freq([i, j], :) = order_freq([j, i], :);
            order_freq_sum([i, j], :) = order_freq_sum([j, i], :);
            biomarker([i, j]) = biomarker([j, i]);
        end
    end
end

% Plotting and saving the figure
h = figure('Visible', 'off'); clf;
imagesc(order_freq), colormap(flipud(gray)); % change the colormap to gray
matrix_string = strtrim(cellstr(num2str(order_freq(:), '%0.2f'))); % convert matrix values to strings
[x, y] = meshgrid(1 : length(order_freq)); % create coordinates for the strings
h_string = text(x(:), y(:), matrix_string(:), 'HorizontalAlignment', 'center', 'FontSize', 6.5); % plot the strings
mean_color = mean(get(gca, 'CLim')); % mean value of the color range
string_color = repmat(order_freq(:) > mean_color, 1, 3); % choose white or black for the text color based on mean color value
set(h_string, {'Color'}, num2cell(string_color, 2)); % change the text colors
set(gca, 'Ticklength', [0 0], 'XTick', [1 : 1 : numel(biomarker)], 'YTick', [1 : 1 : numel(biomarker)], 'YTickLabel', biomarker, 'FontSize', 10, 'FontName', 'times');
xlabel('Event position'), ylabel('Biomarker'), title({'Temporal Ordering', ''});
print(h, '-dpng', 'fig_marker_order.png'); close; % save the figure

end