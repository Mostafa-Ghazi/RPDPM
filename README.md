# RPDPM
Robust Parametric Disease Progression Modeling
<br />

# Description
The robust parametric disease progression modeling (RPDPM) can be applied to time series regression, prediction, and classification tasks to jointly model the trajectories of dynamic features or biomarkers, to predict the temporal developments, and to classify the clinical labels in sequential data with missing values. This toolbox is an implementation of the algorithm proposed in [1].
<br />

# Algorithm
The algorithm, which is based on alternating M-estimation to address potential curve-fitting problems such as outliers, linearly transforms ages to disease progression scores and jointly fits logistic sigmoids to the longitudinal dynamics of biomarkers. The estimated parameters are then used to temporally order the biomarkers in the disease course and to predict biomarker values as well as to classify the clinical status per time point.
<br />

# Dependencies
MATLAB (9.8), Statistics and Machine Learning Toolbox (11.7), Optimization Toolbox (8.5).
<br />

# Functions
<br />

# Inputs
rpdpm_demo.m
<br />
Line 28: a CSV file path containing columns with variable names 'SubjectID', 'Label', 'Age', and biomarkers or features. Missing labels must be assigned as empty cells and missing values need to be specified as NaN.
<br />
Line 72: ranges, the possible minimum and maximum values of each feature. If a feature range is unknown, it should be set as [-Inf, Inf].
<br />
Line 89: ratio_test, proportion of test subjects to all available subjects in data partitioning.
<br />
Line 90: classes, distinct class labels ordered w.r.t. disease progression.
<br />
Line 133: loss_type, robust estimation loss function type.
<br />
Line 134: fit_type, fitting function type.
<br />
Line 137: btstrp, number of bootstraps.
<br />

# Outputs
<br />

# Citation
When you publish your research using this toolbox, please cite [1] as
<br />
<br />
@article{Ghazi2019,
<br />
  title = {Robust parametric modeling of {A}lzheimer's disease progression},
  <br />
  author = {Mehdipour Ghazi, Mostafa and Nielsen, Mads and Pai, Akshay and Modat, Marc and Cardoso, M Jorge and Ourselin, S{\'e}bastien and S{\o}rensen, Lauge},
  <br />
  journal = {arXiv preprint arXiv:1908.05338},
  <br />
  volume = {},
  <br />
  pages = {},
  <br />
  year = {2019},
  <br />
  publisher = {},
  <br />
}
<br />

# References
[1] Mehdipour Ghazi, M., Nielsen, M., Pai, A., Modat, M., Cardoso, M.J., Ourselin, S., and Sørensen, L., 2019. Robust parametric modeling of Alzheimer’s disease progression. arXiv preprint arXiv:1908.05338.
<br />
[2] Jedynak, B.M., Lang, A., Liu, B., Katz, E., Zhang, Y., Wyman, B.T., Raunig, D., Jedynak, C.P., Caffo, B., Prince, J.L., 2012. A computational neurodegenerative disease progression score: method and results with the Alzheimer’s Disease Neuroimaging Initiative cohort. NeuroImage 63, 1478–1486.
<br />

Contact: mostafa.mehdipour@gmail.com
<br />
