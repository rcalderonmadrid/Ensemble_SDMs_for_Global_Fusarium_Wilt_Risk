# Ensemble SDMs for Global Fusarium Wilt Risk Assessment

Species Distribution Modeling (SDM) workflow for predicting global Fusarium wilt risk using ensemble methods, remote sensing and earth system modeling data.

---

## 📂 Project Structure

```
Ensemble_SDMs_for_Global_Fusarium_Wilt/
├── scripts/
│   ├── 00_setup.R                                    # Package installation
│   ├── download_RS_variables.R                       # Download pre-optimized data
│   ├── 01_variable_selection_correlation.R           # Step 1: Correlation filtering
│   ├── 02_variable_optimization_contribution.R       # Step 2: Contribution filtering
│   └── 03_final_modeling_ensemble_projection.R       # Step 3: Final models
├── R/
│   ├── install_requirements.R                        # Package installation script
│   └── utils.R                                       # Utility functions
├── Data/
│   ├── RS_variables/                                 # Environmental predictors
│   ├── occurrences/                                  # Species occurrence data
│   ├── Step 1/                                       # Step 1 outputs
│   ├── Step 2/                                       # Step 2 outputs
│   └── Step 3/                                       # Step 3 outputs
├── requirements.txt                                  # R package dependencies
└── README.md                                         # This file
```

---

## 🚀 Quick Start

### Option A: Full Workflow (Custom Variable Selection)

```r
# 1. Setup environment
source("scripts/00_setup.R")

# 2. Step 1: Correlation-based variable selection (~2-4 hours)
source("scripts/01_variable_selection_correlation.R")

# 3. Step 2: Contribution-based optimization (~3-6 hours)
source("scripts/02_variable_optimization_contribution.R")

# 4. Step 3: Final modeling and ensemble (~8-12 hours)
source("scripts/03_final_modeling_ensemble_projection.R")
```

### Option B: Use Pre-Optimized Variables (Faster)

```r
# 1. Setup environment
source("scripts/00_setup.R")

# 2. Download pre-optimized variables (~10 minutes)
source("scripts/download_RS_variables.R")

# 3. Skip to Step 3: Final modeling
source("scripts/03_final_modeling_ensemble_projection.R")
```

---

## 📋 Workflow Description

### 00. Setup (`00_setup.R`)
- Reads package list from `requirements.txt`
- Calls `R/install_requirements.R` to install all required R packages
- Configures R environment
- Loads utility functions from `R/utils.R`

**Duration:** 10-30 minutes  
**Resources:** 1 core, 4GB RAM

---

### Download Data (`download_RS_variables.R`)
**Purpose:** Download pre-optimized RS variables from Figshare (optional)

**Process:**
- Downloads `.7z` archive from Figshare
- Extracts environmental raster layers
- Validates file integrity

**Use case:** Skip Steps 1-2 if you want to use pre-processed variables

**Outputs:**
- `Data/RS_variables/*.tif` - Environmental predictor rasters

**Duration:** 5-15 minutes  
**Resources:** 1 core, 2GB RAM, 500MB disk space

---

### Step 1: Correlation-Based Variable Selection
**Script:** `01_variable_selection_correlation.R`

**Purpose:** Remove highly correlated variables to reduce multicollinearity

**Process:**
1. **Occurrence filtering**
   - Spatial thinning (1 km minimum distance)
   - Geographic filtering to reduce clustering
   
2. **Spatial block cross-validation**
   - Calculates spatial autocorrelation range
   - Creates 5 spatial blocks for independent CV
   
3. **Pseudo-absence generation**
   - 10 replicates with environmental/geographic constraints
   - Equal distribution across spatial blocks
   
4. **Initial model calibration**
   - Trains 6 algorithms: ANN, GBM, MAXNET, RF, GAM, MARS
   - Uses default parameters ("bigboss" strategy)
   
5. **Correlation filtering**
   - Computes variable correlations (Spearman, threshold = 0.7)
   - Uses jackknife test to determine which correlated variable contributes less
   - Iteratively removes correlated variables

**Outputs:**
- `Data/Step 1/vars_corr_removed.csv` - Uncorrelated variable set
- `Data/Step 1/filtered_occ.csv` - Filtered occurrences
- `Data/Step 1/Model_initial.rds` - Initial model object
- `Data/Step 1/Pseudoabsences/psa_*.csv` - Pseudo-absence datasets

**Typical results:** ~50 variables → ~30 uncorrelated variables

**Duration:** 2-4 hours  
**Resources:** 15-25 cores, 32GB RAM, 5GB disk space

---

### Step 2: Contribution-Based Variable Optimization
**Script:** `02_variable_optimization_contribution.R`

**Purpose:** Remove low-contribution variables to improve model parsimony

**Process:**
1. **Uses Step 1 outputs**
   - Loads uncorrelated variables from Step 1
   - Applies same occurrence filtering workflow
   
2. **New pseudo-absence generation**
   - Fresh spatial blocks and pseudo-absences
   - Independent from Step 1 for robust evaluation
   
3. **Model calibration**
   - Trains models with uncorrelated variables
   - Calculates variable importance (permutation-based)
   
4. **Contribution optimization**
   - Identifies variables with contribution ≤ 2%
   - Uses jackknife test to evaluate removal impact
   - Removes variable if validation performance is maintained/improved
   - Iterates until all remaining variables are important

**Outputs:**
- `Data/Step 2/vars_optimized.csv` - Final optimized variable set
- `Data/Step 2/filtered_occ.csv` - Filtered occurrences
- `Data/Step 2/Model_initial.rds` - Model with optimized variables
- `Data/Step 2/Pseudoabsences/psa_*.csv` - Pseudo-absence datasets

**Typical results:** ~30 variables → ~15-20 optimized variables

**Duration:** 3-6 hours  
**Resources:** 15-25 cores, 32GB RAM, 8GB disk space

---

### Step 3: Final Modeling, Ensemble & Projection
**Script:** `03_final_modeling_ensemble_projection.R`

**Purpose:** Final model calibration with tuned parameters and ensemble predictions

**Process:**
1. **Uses Step 2 outputs**
   - Loads optimized variables
   - Applies occurrence filtering workflow
   
2. **Enhanced pseudo-absence generation**
   - 11 replicates: 10 for training + 1 for independent evaluation
   - Provides unbiased model assessment
   
3. **Hyperparameter tuning**
   - Grid search for each algorithm:
     - **ANN:** Hidden units, weight decay
     - **GAM:** Term selection, smoothing method
     - **GBM:** Trees, depth, learning rate
     - **MARS:** Interaction degree, pruning
     - **RF:** Number of predictors (mtry)
     - **MAXNET:** Uses optimized defaults
   - Evaluation metric: TSS (True Skill Statistic)
   
4. **Final model calibration**
   - Trains models with tuned hyperparameters
   - Spatial block cross-validation (5 folds × 10 replicates)
   - Extended evaluation metrics: ROC, TSS, FAR, BIAS, POD, POFD, SR, BOYCE
   
5. **Ensemble modeling**
   - Filters models by TSS ≥ 0.8
   - 6 ensemble strategies:
     - **EMmean:** Simple mean
     - **EMmedian:** Median prediction
     - **EMcv:** CV-weighted mean
     - **EMci:** Confidence interval
     - **EMca:** Committee averaging
     - **EMwmean:** Performance-weighted mean
   
6. **Spatial projections**
   - Projects individual models onto current conditions
   - Projects ensemble models
   - Outputs as GeoTIFF rasters

**Outputs:**
- `Data/Step 3/Tuning/tuned_*.rds` - Tuned hyperparameters
- `Data/Step 3/Tuned.ModelOut.rds` - Final tuned models
- `Data/Step 3/Ensemble_model.rds` - Ensemble model object
- `Data/Step 3/model_scores.csv` - Individual model evaluation
- `Data/Step 3/EM_scores.csv` - Ensemble evaluation
- `Data/Step 3/Single_models_proj.rds` - Individual projections
- `Data/Step 3/Ensemble_models_proj.rds` - Ensemble projections
- `Data/Step 3/Foxy/proj_Present_RS/*.tif` - Projection rasters

**Duration:** 8-12 hours  
**Resources:** 20-30 cores, 64GB RAM, 15GB disk space

---

## 📊 Variable Selection Pipeline

```
All RS Variables (n = ~50)
         ↓
┌────────────────────────────────────┐
│ Step 1: Correlation Filtering     │
│ - Removes |r| ≥ 0.7               │
│ - Jackknife test for selection    │
└────────────────────────────────────┘
         ↓
Uncorrelated Variables (n = ~30)
         ↓
┌────────────────────────────────────┐
│ Step 2: Contribution Optimization │
│ - Removes variables ≤ 2%          │
│ - Maintains/improves performance  │
└────────────────────────────────────┘
         ↓
Optimized Variables (n = ~15-20)
         ↓
┌────────────────────────────────────┐
│ Step 3: Final Modeling            │
│ - Tuned hyperparameters           │
│ - Ensemble predictions            │
└────────────────────────────────────┘
```

---

## 💻 Computational Requirements

| Step | CPU Cores | RAM | Time | Disk Space |
|------|-----------|-----|------|------------|
| **00. Setup** | 1 | 4GB | 10-30 min | 2GB |
| **Download** | 1 | 2GB | 5-15 min | 500MB |
| **Step 1** | 15-25 | 32GB | 2-4 hours | 5GB |
| **Step 2** | 15-25 | 32GB | 3-6 hours | 8GB |
| **Step 3** | 20-30 | 64GB | 8-12 hours | 15GB |

**Recommendations:**
- Use HPC cluster for Steps 1-3
- Adjust `nb.cpu` and `mc.cores` parameters based on available resources
- Monitor memory usage during spatial operations

---

## 📦 Required Data

### Input Data

1. **Environmental Variables** (`Data/RS_variables/`)
   - Remote sensing predictor rasters (.tif format)
   - Global coverage, 1 km resolution
   - **Note:** Not included in repository due to file size
   - **Download:** Available on Figshare (DOI: 10.6084/m9.figshare.30416476)
   - Alternatively, use `download_RS_variables.R` script

2. **Occurrence Data** (`Data/occurrences/occurrences.csv`)
   - Species presence records
   - Required columns: `x` (longitude), `y` (latitude)
   - Coordinate system: WGS84 (EPSG:4326)
   - **Download:** Available on Figshare (DOI: 10.6084/m9.figshare.30416476)

### Output Data Structure

```
Data/
├── Step 1/
│   ├── vars_corr_removed.csv           # Uncorrelated variables
│   ├── filtered_occ.csv                # Filtered occurrences
│   ├── spatial_autocorrelation.rds     # Autocorrelation metrics
│   ├── block_partition.tif             # Spatial CV blocks
│   ├── Model_initial.rds               # Initial models
│   └── Pseudoabsences/
│       └── psa_*.csv                   # Pseudo-absence replicates
│
├── Step 2/
│   ├── vars_optimized.csv              # Optimized variables
│   ├── filtered_occ.csv                # Filtered occurrences
│   ├── Model_initial.rds               # Models with optimized vars
│   └── Pseudoabsences/
│       └── psa_*.csv                   # Pseudo-absence replicates
│
└── Step 3/
    ├── Tuned.ModelOut.rds              # Tuned models
    ├── Ensemble_model.rds              # Ensemble models
    ├── model_scores.csv                # Evaluation metrics
    ├── EM_scores.csv                   # Ensemble metrics
    ├── Tuning/
    │   └── tuned_*.rds                 # Hyperparameters
    ├── Pseudoabsences/
    │   └── psa_*.csv                   # 11 replicates
    └── Foxy/
        └── proj_*/*.tif                # Spatial projections
```

---

## 🔧 Configuration

### Adjusting Computational Resources

Edit these parameters in each script:

```r
# Parallel processing cores
nb.cpu = 15        # For biomod2 functions
mc.cores = 10      # For mclapply functions
workers = 2        # For future/foreach

# Memory management (terra)
terra::terraOptions(
  memmax = 100,    # Maximum memory in GB
  steps = 500      # Processing chunks
)

# Future package (for large objects)
options(future.globals.maxSize = 25000 * 1024^2)  # 25 GB
```

### Modifying Filtering Thresholds

**Step 1 - Correlation threshold:**
```r
RemoveCorrVar(..., cor_th = 0.7)  # Default: 0.7
```

**Step 2 - Contribution threshold:**
```r
OptimizeVar(..., th = 2)  # Default: 2%
```

**Step 3 - Model selection threshold:**
```r
BIOMOD_EnsembleModeling(..., metric.select.thresh = 0.8)  # Default: TSS ≥ 0.8
```

---

## 📚 Key References

- **Occurrence filtering:** Velazco et al. (2021) - Environmental filtering methodology
- **Spatial blocks:** Valavi et al. (2019) - Spatial cross-validation for SDMs
- **Pseudo-absences:** Barbet-Massin et al. (2012) - Pseudo-absence generation strategies
- **Ensemble models:** Araújo & New (2007) - Ensemble forecasting in species distribution modeling

---

## 👥 Contact

**Rocío Calderón Madrid**  
Email: r.calderon@csic.es

For questions about:
- RS variables or data access
- Methodology implementation
- Results interpretation

---

## 📄 License

[Add your license information here]

---

## 🙏 Acknowledgments

[Add acknowledgments here]

---

## Installation and Setup

### Prerequisites

- **R** (≥ 4.3.0)
- **Operating System**: Linux/macOS recommended for HPC; Windows possible with limitations
- **Disk space**: ~50 GB for complete workflow (input data + outputs)
- **Computational resources**: See requirements table below

### Step 1: Clone the Repository

```bash
git clone https://github.com/[your-username]/Ensemble_SDMs_for_Global_Fusarium_Wilt_Risk.git
cd Ensemble_SDMs_for_Global_Fusarium_Wilt_Risk
```

### Step 2: Install Required R Packages

The repository includes a `requirements.txt` file listing all necessary R packages.

**Option A - Direct installation:**
```r
source("R/install_requirements.R")
```

**Option B - Through setup script:**
```r
source("scripts/00_setup.R")
```

Both options read from `requirements.txt` and install the following packages:
- **SDM frameworks**: `biomod2`, `flexsdm`, `terra`
- **Machine learning**: `randomForest`, `gbm`, `mgcv`, `earth`, `nnet`, `maxnet`
- **Parallel processing**: `foreach`, `doFuture`, `future`, `parallel`
- **Data manipulation**: `dplyr`, `tidyr`
- **Spatial analysis**: `ape`, `sf`
- **Utilities**: `devtools`, `ggplot2`, `e1071`, `automap`

**Installation time**: 10–30 minutes depending on your system

### Step 3: Download Input Data from Figshare

1. Visit the Figshare repository: [https://doi.org/10.6084/m9.figshare.30416476](https://doi.org/10.6084/m9.figshare.30416476)
2. Download the following datasets:
   - **Occurrence records** → place in `Data/occurrences/`
   - **Environmental variables** (RS predictors) → place in `Data/RS_variables/`
3. Ensure file paths match those expected by the scripts:
   - Occurrence file: `Data/occurrences/occurrences.csv`
   - RS variables: `Data/RS_variables/*.tif`

**Download size**: ~2–5 GB compressed

---

## Workflow Execution

The analysis follows a **sequential three-step workflow**. Run scripts in numerical order:

### Step 1: Correlation-Based Variable Selection

**Script**: `01_variable_selection_correlation.R`

**Purpose**: Remove highly correlated environmental predictors to reduce multicollinearity.

```r
source("scripts/01_variable_selection_correlation.R")
```

**Key processes**:
1. **Occurrence filtering**:
   - Spatial thinning (1 km minimum distance between points)
   - Geographic filtering to reduce sampling bias
   
2. **Spatial block cross-validation**:
   - Calculate spatial autocorrelation range in environmental space
   - Generate 5 spatial blocks ensuring independence between folds
   
3. **Pseudo-absence generation**:
   - Create 10 replicate sets with environmental and geographic constraints
   - Equal distribution across spatial blocks
   
4. **Initial model calibration**:
   - Train 6 algorithms (ANN, GBM, MAXNET, RF, GAM, MARS) with default parameters
   - Evaluate using spatial cross-validation
   
5. **Correlation filtering**:
   - Identify variable pairs with |r| ≥ 0.7 (Spearman correlation)
   - Use jackknife test to determine which variable contributes less
   - Iteratively remove correlated variables until threshold met

**Outputs**:
- `vars_corr_removed.csv`: Uncorrelated variable set (~30 variables)
- `filtered_occ.csv`: Spatially and environmentally filtered occurrences
- `block_partition.tif`: Spatial blocks for cross-validation
- `Model_initial.rds`: Initial model object for evaluation
- `Pseudoabsences/psa_*.csv`: 10 pseudo-absence replicates

**Duration**: 2–4 hours  
**Resources**: 15–25 CPU cores, 32 GB RAM, 5 GB disk space

---

### Step 2: Contribution-Based Variable Optimization

**Script**: `02_variable_optimization_contribution.R`

**Purpose**: Remove low-contribution variables to improve model parsimony and predictive performance.

```r
source("scripts/02_variable_optimization_contribution.R")
```

**Key processes**:
1. **Load Step 1 outputs**:
   - Use uncorrelated variable set from Step 1
   - Apply same occurrence filtering workflow (independent random seed)
   
2. **New spatial structure**:
   - Generate fresh spatial blocks and pseudo-absences
   - Ensures independent evaluation from Step 1
   
3. **Model calibration**:
   - Train models using uncorrelated variables
   - Calculate variable importance via permutation
   
4. **Contribution optimization**:
   - Identify variables with ≤2% contribution (permutation importance)
   - Jackknife test: remove variable if performance maintained/improved
   - Iterate until all remaining variables are ecologically important

**Outputs**:
- `vars_optimized.csv`: Final optimized variable set (~15–20 variables)
- `filtered_occ.csv`: Filtered occurrences (new random filtering)
- `block_partition.tif`: New spatial blocks
- `Model_initial.rds`: Model trained with optimized variables
- `Pseudoabsences/psa_*.csv`: 10 new pseudo-absence replicates

**Duration**: 3–6 hours  
**Resources**: 15–25 CPU cores, 32 GB RAM, 8 GB disk space

---

### Step 3: Final Modeling, Ensemble & Spatial Projection

**Script**: `03_final_modeling_ensemble_projection.R`

**Purpose**: Final model calibration with tuned hyperparameters, ensemble modeling, and spatial projections.

```r
source("scripts/03_final_modeling_ensemble_projection.R")
```

**Key processes**:

**Part 1–4: Data preparation**
- Load optimized variables from Step 2
- Apply occurrence filtering (spatial + environmental with Moran's I optimization)
- Generate 11 pseudo-absence replicates (10 training + 1 independent evaluation)
- Create spatial blocks for cross-validation

**Part 5–6: Data formatting**
- Format data for biomod2 with evaluation set
- Set up spatial block cross-validation structure

**Part 7: Hyperparameter tuning**
- Grid search for each algorithm:
  - **ANN**: Hidden units (2–35), weight decay (0.01–0.3)
  - **GAM**: Term selection, smoothing methods (GCV.Cp, REML, ML, etc.)
  - **GBM**: Trees (500–2500), depth (2–10), learning rate (0.001–0.1)
  - **MARS**: Interaction degree (1–2), pruning parameters
  - **RF**: Number of predictors to sample (mtry)
  - **MAXNET**: Uses optimized defaults
- Evaluation metric: TSS (True Skill Statistic)

**Part 8: Model calibration**
- Train all algorithms with tuned hyperparameters
- Spatial cross-validation: 5 folds × 10 replicates = 50 models per algorithm
- Calculate 8 evaluation metrics: ROC, TSS, FAR, BIAS, POD, POFD, SR, BOYCE
- Compute variable importance (10 permutations)

**Part 9: Ensemble modeling**
- Filter models by TSS ≥ 0.8 on validation data
- Build 6 ensemble strategies:
  - **EMmean**: Simple mean of predictions
  - **EMmedian**: Median of predictions
  - **EMcv**: Coefficient of variation weighted mean
  - **EMci**: Confidence interval (α = 0.05)
  - **EMca**: Committee averaging
  - **EMwmean**: Performance-weighted mean (proportional decay)

**Part 10: Spatial projections**
- Project individual models onto present-day environmental conditions
- Project ensemble models
- Output continuous suitability and binary presence/absence maps (TSS threshold)

**Outputs**:
- `spatial_autocorrelation.rds`: Autocorrelation metrics for block size
- `filtered_occ.csv`: Final filtered occurrences
- `FoxyBiomodData.rds`: Formatted biomod2 object with evaluation data
- `Tuning/tuned_*.rds`: Optimized hyperparameters for each algorithm
- `Tuned.ModelOut.rds`: Final calibrated models (300 models total)
- `Ensemble_model.rds`: Ensemble model object
- `model_scores.csv`: Evaluation metrics for all individual models
- `summary_model_scores.csv`: Summary statistics by algorithm
- `var_imp.csv`: Variable importance scores
- `summary_models_var_imp.csv`: Normalized variable importance summary
- `EM_scores.csv`: Ensemble model evaluation metrics
- `EM_var_imp.csv`: Ensemble variable importance
- `Single_models_proj.rds`: Individual model projections
- `Ensemble_models_proj.rds`: Ensemble projections
- `Foxy/proj_Present_RS/*.tif`: Spatial projection rasters (GeoTIFF)
- `Foxy/proj_EM_Present_RS/*.tif`: Ensemble projection rasters

**Duration**: 8–12 hours  
**Resources**: 20–30 CPU cores, 64 GB RAM, 15 GB disk space

---

## Computational Requirements

| Step | CPU Cores | RAM | Time | Disk Space |
|------|-----------|-----|------|------------|
| **00. Setup** | 1 | 4 GB | 10–30 min | 2 GB |
| **Step 1** | 15–25 | 32 GB | 2–4 hours | 5 GB |
| **Step 2** | 15–25 | 32 GB | 3–6 hours | 8 GB |
| **Step 3** | 20–30 | 64 GB | 8–12 hours | 15 GB |
| **Total** | — | — | **13–22 hours** | **~30 GB** |

**Recommendations**:
- Use HPC cluster or high-performance workstation
- Adjust `nb.cpu` and `mc.cores` parameters based on available resources
- Monitor memory during raster operations (use `terra::terraOptions()`)
- Steps can be run separately; intermediate results are saved

**Memory optimization**:
```r
# In each script, adjust based on your system:
terra::terraOptions(
  memmax = 100,      # Maximum memory in GB
  steps = 500        # Processing chunks for large rasters
)

# For parallel processing with large objects:
options(future.globals.maxSize = 25000 * 1024^2)  # 25 GB
```

---

## Variable Selection Pipeline Summary

```
Initial RS Variables (n ≈ 50)
         ↓
┌─────────────────────────────────────────┐
│ Step 1: Correlation Filtering          │
│ • Remove |r| ≥ 0.7 (Spearman)          │
│ • Jackknife test for variable selection│
│ • Maintain model performance           │
└─────────────────────────────────────────┘
         ↓
Uncorrelated Variables (n ≈ 30)
         ↓
┌─────────────────────────────────────────┐
│ Step 2: Contribution Optimization      │
│ • Remove variables ≤2% importance      │
│ • Jackknife test                       │
│ • Improve parsimony & generalization   │
└─────────────────────────────────────────┘
         ↓
Optimized Variables (n ≈ 15–20)
         ↓
┌─────────────────────────────────────────┐
│ Step 3: Final Modeling                 │
│ • Hyperparameter tuning                │
│ • Spatial cross-validation             │
│ • Ensemble predictions                 │
│ • Present-day projections              │
└─────────────────────────────────────────┘
```

---

## Future Climate Projections

**Note**: Scripts for CMIP6 ESM projections use the trained ensemble models from Step 3. Future projection outputs are available on Figshare.

**Projection setup**:
- **Time periods**: 2041–2060 (mid-century), 2081–2100 (end-century)
- **Scenarios**: SSP2-4.5 (moderate mitigation), SSP5-8.5 (high emissions)
- **ESMs**: Multi-model ensemble from CMIP6 archive (bias-corrected)
- **Output**: Individual ESM projections + multi-model ensemble mean
- **Analysis**: Range change maps (stable, expansion, contraction zones)

---

## Citation

If you use this workflow or data in your research, please cite:

> [Author list]. ([Year]). Global environmental suitability and potential establishment of *Fusarium oxysporum* f. sp. *cubense* TR4 under present and future climate conditions. *[Journal Name]*, *[Volume]*([Issue]), [Pages]. DOI: [manuscript DOI]

**Data repository:**
> [Author list]. ([Year]). Dataset: Ensemble SDMs for global Fusarium wilt risk assessment. *Figshare*. DOI: [10.6084/m9.figshare.30416476](https://doi.org/10.6084/m9.figshare.30416476)

---

## License

This project is licensed under the [MIT License](LICENSE) (or specify your license). Data on Figshare are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

---

## Contact

For questions or issues, please open an issue on GitHub or contact:

- **Lead author**: [Name] ([email@institution.edu])
- **Repository maintainer**: [Name] ([email@institution.edu])

---

## Acknowledgments

This research was supported by [funding sources]. We acknowledge the World Climate Research Programme's Working Group on Coupled Modelling for coordinating CMIP6 data, and the climate modeling centers for producing and sharing their model output.

Remote sensing data were obtained from [NASA MODIS/Landsat archives, etc.]. Occurrence data were compiled from [GBIF, CABI, literature sources, etc.].

---
