# Ensemble SDMs for Global Fusarium Wilt Risk Assessment

Species Distribution Modeling (SDM) workflow for predicting global Fusarium wilt risk using ensemble methods and remote sensing data.

---

## 📂 Project Structure

```
Ensemble_SDMs_for_Global_Fusarium_Wilt_Risk/
├── scripts/
│   ├── 00_setup.R                                    # Package installation
│   ├── download_RS_variables.R                       # Download pre-optimized data
│   ├── 01_variable_selection_correlation.R           # Step 1: Correlation filtering
│   ├── 02_variable_optimization_contribution.R       # Step 2: Contribution filtering
│   └── 03_final_modeling_ensemble_projection.R       # Step 3: Final models
├── R/
│   └── utils.R                                       # Utility functions
├── Data/
│   ├── RS_variables/                                 # Environmental predictors
│   ├── occurrences/                                  # Species occurrence data
│   ├── Step 1/                                       # Step 1 outputs
│   ├── Step 2/                                       # Step 2 outputs
│   └── Step 3/                                       # Step 3 outputs
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
- Installs required R packages from `requirements.txt`
- Configures R environment
- Loads utility functions

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
   - Global coverage, consistent resolution
   - **Note:** Not included in repository due to file size
   - **Download:** Use `download_RS_variables.R` or contact r.calderon@csic.es

2. **Occurrence Data** (`Data/occurrences/occurrences.csv`)
   - Species presence records
   - Required columns: `x` (longitude), `y` (latitude)
   - Coordinate system: WGS84 (EPSG:4326)

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
