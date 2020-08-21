function mauc = multiclass_auc(scores, labels, unique_labels)

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

% Classification performance based on multiclass area under the ROC curve

% N: number of labels
% C: number of unique class labels
% scores: matrix of posterior probabilities [N x C]
% labels: vector of true labels for a multiclass problem [N x 1]
% unique_labels: unique class labels corresponding to the scores columns [1 x C]
% mauc: average area under the curve (one-vs-one)

C = numel(unique_labels); % number of unique class labels
CC = C * (C - 1) / 2; % number of combinations of class pairs

% Calculating the accuracy
mauc = 0; % average auc
for c = 1 : C - 1
    
    idx_c = find(strcmp(labels, unique_labels(c)));
    n_c = numel(idx_c); % number of within class labels
    
    for cc = c + 1 : C
        
        idx_cc = find(strcmp(labels, unique_labels(cc)));
        n_cc = numel(idx_cc); % number of within class labels
        
        s_c = scores([idx_c; idx_cc], c);
        s_cc = scores([idx_c; idx_cc], cc);
        
        [~, idx_sort_c] = sort(s_c, 'ascend'); % sorts scores of both classes
        [~, rank_c] = sort(idx_sort_c, 'ascend'); % ranks of classes
        
        [~, idx_sort_cc] = sort(s_cc, 'ascend'); % sorts scores of both classes
        [~, rank_cc] = sort(idx_sort_cc, 'ascend'); % ranks of classes
        
        if (n_c ~= 0 && n_cc ~= 0)
            auc_c = (sum(rank_c(1 : n_c)) - n_c * (n_c + 1) / 2) / (n_c * n_cc);
            auc_cc = (sum(rank_cc(n_c + 1 : n_c + n_cc)) - n_cc * (n_cc + 1) / 2) / (n_c * n_cc);
            mauc = mauc + (auc_c + auc_cc) / 2;
        end
        
    end
    
end
mauc = mauc / CC;

end