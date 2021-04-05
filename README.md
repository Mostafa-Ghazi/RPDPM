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
MATLAB (tested with v9.8), Statistics and Machine Learning Toolbox (tested with v11.7), Optimization Toolbox (tested with v8.5).
<br />

# Inputs
•	A CSV file containing age information, labels, and measurements in columns under variable names 'SubjectID', 'Label', 'Age', and 'Features'. Missing labels and missing values need to be assigned as empty cell and NaN, respectively.
<br />
•	The minimum and maximum values of each feature. If a feature range is unknown, it should be set as [-Inf, Inf].
<br />
•	Proportion of test subjects to all available subjects in data partitioning.
<br />
•	Distinct class labels ordered w.r.t. disease progression.
<br />
•	Robust estimation loss function type.
<br />
•	Fitting function type.
<br />
•	The minimum and maximum number of alternating iterations.
<br />
•	Number of bootstraps.
<br />

# Outputs
•	Training performance across all bootstraps (BIC) printed to the command window.
<br />
•	Validation performance across all bootstraps (NMAE) printed to the command window.
<br />
•	Testing performance across all bootstraps (NMAE and AUC) printed to the command window.
<br />
•	A saved figure displaying the estimated class-conditional likelihoods using the DPSs.
<br />
•	A saved figure displaying the temporal ordering of biomarkers in the disease course.
<br />
•	Saved figures displaying the estimated trajectories per biomarker.
<br />

# Citation
When you publish your research using this toolbox, please cite [1] as
<br />
<br />
@article{Ghazi2020,
<br />
  title = {Robust parametric modeling of {A}lzheimer's disease progression},
  <br />
  author = {Mehdipour Ghazi, Mostafa and Nielsen, Mads and Pai, Akshay and Modat, Marc and Cardoso, M Jorge and Ourselin, S{\'e}bastien and S{\o}rensen, Lauge},
  <br />
  journal = {NeuroImage},
  <br />
  volume = {225},
  <br />
  pages = {117460},
  <br />
  year = {2021},
  <br />
  publisher = {Elsevier},
  <br />
}
<br />

# References
[1] Mehdipour Ghazi, M., Nielsen, M., Pai, A., Modat, M., Cardoso, M.J., Ourselin, S., and Sørensen, L., 2019. Robust parametric modeling of Alzheimer’s disease progression. NeuroImage 225, 117460.
<br />
[2] Jedynak, B.M., Lang, A., Liu, B., Katz, E., Zhang, Y., Wyman, B.T., Raunig, D., Jedynak, C.P., Caffo, B., Prince, J.L., 2012. A computational neurodegenerative disease progression score: method and results with the Alzheimer’s Disease Neuroimaging Initiative cohort. NeuroImage 63, 1478–1486.
<br />

Contact: mostafa.mehdipour@gmail.com
<br />
