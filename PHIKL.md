# PHIKL Pharmacokinetic Modeling System Documentation

This document provides detailed information about using NONMEM and MONOLIX control streams for pharmacokinetic parameter estimation using the SAEM and FOCE algorithms.

## Table of Contents

1. [Overview](#overview)
2. [NONMEM Control Stream Format](#nonmem-control-stream-format)
3. [MONOLIX Control Stream Format](#monolix-control-stream-format)
4. [Supported Model Types](#supported-model-types)
5. [Estimation Methods](#estimation-methods)
6. [Data File Format](#data-file-format)
7. [Complete Examples](#complete-examples)
8. [Parsing Details](#parsing-details)
9. [Best Practices](#best-practices)
10. Machine Learning Models

---

## Overview

This system provides a Rust-based pharmacokinetic modeling framework that can parse and execute both NONMEM and MONOLIX control streams. It supports:

- **Model Types**: One, two, and three-compartment pharmacokinetic models, plus machine learning models (Random Forest, XGBoost)
- **Estimation Methods**: FOCE (First-Order Conditional Estimation), FOCE-I (with interaction), and SAEM (Stochastic Approximation Expectation-Maximization)
- **Machine Learning Methods**: Random Forest and XGBoost (Extreme Gradient Boosting) for data-driven modeling
- **Control Stream Formats**: Both NONMEM (.ctl, .mod) and MONOLIX (.mlxtran) formats

---

## NONMEM Control Stream Format

NONMEM control streams use a block-based structure with dollar signs ($) denoting each section.

### Basic Structure

```
$PROBLEM [Description]
$INPUT [Column definitions]
$DATA [Data file path]
$SUBROUTINES [ADVAN routine]
$PK [Pharmacokinetic model]
$THETA [Fixed effects initial estimates]
$OMEGA [Random effects variance]
$SIGMA [Residual error variance]
$ERROR [Error model]
$ESTIMATION [Estimation method]
$COVARIANCE [Covariance calculation]
$TABLE [Output tables]
```

### Key Sections

#### $DATA Section

Specifies the dataset file to use:

```
$DATA example_dataset.csv IGNORE=C
$DATA /path/to/data.csv
```

- `IGNORE=C`: Skip rows starting with 'C' (typically headers/comments)

#### $SUBROUTINES Section

Defines the compartmental model structure:

```
$SUBROUTINES ADVAN1 TRANS2    ; One-compartment
$SUBROUTINES ADVAN3 TRANS4    ; Two-compartment
$SUBROUTINES ADVAN11 TRANS4   ; Three-compartment
```

**Common ADVAN routines:**
- `ADVAN1`: One-compartment model
- `ADVAN2`: One-compartment with first-order absorption
- `ADVAN3`: Two-compartment model
- `ADVAN4`: Two-compartment with first-order absorption
- `ADVAN11`: Three-compartment model
- `ADVAN12`: Three-compartment with first-order absorption

#### $PK Section

Defines the pharmacokinetic parameters:

```
$PK
  CL = THETA(1) * EXP(ETA(1))
  V  = THETA(2) * EXP(ETA(2))
  S1 = V
```

- `THETA(n)`: Population (fixed effect) parameters
- `ETA(n)`: Individual random effects
- `S1`: Scaling factor for compartment 1

For two-compartment models:

```
$PK
  CL = THETA(1) * EXP(ETA(1))
  V1 = THETA(2) * EXP(ETA(2))
  Q  = THETA(3)
  V2 = THETA(4)
  S1 = V1
```

#### $THETA Section

Initial estimates for population parameters:

**Format 1: Lower bound, initial, upper bound**
```
$THETA
  (0, 2.0)      ; CL (L/h)
  (0, 20.0)     ; V (L)
  (0, 0.1, 1)   ; Ka (1/h) with upper bound
```

**Format 2: Single value**
```
$THETA
  2.0           ; CL (L/h)
  20.0          ; V (L)
```

#### $OMEGA Section

Between-subject variability (variance of random effects):

**Diagonal matrix (independent random effects):**
```
$OMEGA
  0.09          ; IIV on CL (variance)
  0.04          ; IIV on V (variance)
```

**Block matrix (correlated random effects):**
```
$OMEGA BLOCK(2)
  0.09          ; IIV on CL
  0.01 0.04     ; Correlation CL-V, IIV on V
```

#### $SIGMA Section

Residual error variance:

```
$SIGMA
  0.01          ; Proportional error variance
  0.5           ; Additive error variance
```

#### $ERROR Section

Defines the error model:

```
$ERROR
  IPRED = F
  IRES = DV - IPRED
  W = IPRED * THETA(3) + THETA(4)
  IF (W.EQ.0) W = 1
  IWRES = IRES / W
  Y = IPRED + W * ERR(1)
```

- `F`: Model prediction
- `IPRED`: Individual prediction
- `DV`: Dependent variable (observed value)
- `W`: Weighting factor
- `ERR(1)`: Residual error

#### $ESTIMATION Section

Specifies the estimation method:

**FOCE (First-Order Conditional Estimation):**
```
$ESTIMATION METHOD=1 MAXEVAL=9999 PRINT=5
```

**FOCE with interaction:**
```
$ESTIMATION METHOD=1 MAXEVAL=9999 INTERACTION PRINT=5
```

**SAEM:**
```
$ESTIMATION METHOD=SAEM MAXEVAL=9999 PRINT=5
```

Options:
- `METHOD=1`: FOCE
- `INTERACTION`: Include eta-epsilon interaction
- `MAXEVAL`: Maximum evaluations
- `PRINT`: Print frequency

---

## MONOLIX Control Stream Format

MONOLIX uses an XML-like structure with sections enclosed in angle brackets.

### Basic Structure

```
<DATAFILE>
[File information and content mapping]

<MODEL>
[Individual and longitudinal models]

<FIT>
[Data to model mapping]

<PARAMETER>
[Initial parameter estimates]

<MONOLIX>
[Tasks and settings]
```

### Key Sections

#### <DATAFILE> Section

Defines the dataset structure:

```
<DATAFILE>

[FILEINFO]
file = 'example_dataset.csv'
delimiter = comma
header = {ID, TIME, DV, AMT, EVID, CMT, WEIGHT, AGE}

[CONTENT]
ID = {use=identifier}
TIME = {use=time}
DV = {use=observation, name=concentration, type=continuous}
AMT = {use=amount}
EVID = {use=eventidentifier}
CMT = {use=compartment}
WEIGHT = {use=covariate, type=continuous}
AGE = {use=covariate, type=continuous}
```

**Column Uses:**
- `identifier`: Subject ID
- `time`: Time of observation/dose
- `observation`: Dependent variable
- `amount`: Dose amount
- `eventidentifier`: Event type (0=observation, 1=dose)
- `compartment`: Compartment number
- `covariate`: Covariate variable

#### <MODEL> Section

Defines individual and longitudinal models:

**[INDIVIDUAL] subsection:**
```
[INDIVIDUAL]
input = {CL_pop, omega_CL, V_pop, omega_V}

DEFINITION:
CL = {distribution=logNormal, typical=CL_pop, sd=omega_CL}
V = {distribution=logNormal, typical=V_pop, sd=omega_V}
```

**Distribution types:**
- `logNormal`: Log-normal distribution (common for PK parameters)
- `normal`: Normal distribution

**[LONGITUDINAL] subsection:**
```
[LONGITUDINAL]
input = {a, b}

file = 'lib:oral1_1cpt_kaVCl.txt'

DEFINITION:
concentration = {distribution=normal, prediction=Cc, errorModel=combined1(a, b)}
```

**Library model files:**
- `lib:oral1_1cpt_kaVCl.txt`: One-compartment with first-order absorption
- `lib:oral2_2cpt_kaClV1QV2.txt`: Two-compartment with absorption
- `lib:bolus1_1cpt_VCl.txt`: One-compartment IV bolus
- `lib:infusion2_2cpt_ClV1QV2.txt`: Two-compartment IV infusion

**Error models:**
- `constant`: Additive error (a)
- `proportional`: Proportional error (a)
- `combined1`: Combined error (a*prediction + b)
- `combined2`: Combined error (sqrt(a^2 + (b*prediction)^2))

#### <PARAMETER> Section

Initial parameter values:

```
<PARAMETER>
CL_pop = {value=2.0, method=MLE}
V_pop = {value=20.0, method=MLE}
omega_CL = {value=0.3, method=MLE}
omega_V = {value=0.2, method=MLE}
a = {value=0.1, method=MLE}
b = {value=0.5, method=MLE}
```

- `value`: Initial estimate
- `method=MLE`: Maximum likelihood estimation

#### <MONOLIX> Section

Tasks and settings:

```
<MONOLIX>

[TASKS]
populationParameters()
individualParameters(method = {conditionalMean, conditionalMode})
fim(method = StochasticApproximation)
logLikelihood(method = ImportanceSampling)

[SETTINGS]
GLOBAL:
exportpath = 'output_directory'

POPULATION:
exploratoryautostop = no
smoothingautostop = no
```

---

## Supported Model Types

### One-Compartment Model

**Differential Equations:**
```
dA/dt = -CL/V * A
C = A/V
```

**Parameters:**
- `CL`: Clearance (L/h)
- `V`: Volume of distribution (L)
- `Ka`: Absorption rate constant (1/h) - for oral dosing

**NONMEM Example:**
```
$SUBROUTINES ADVAN1 TRANS2

$PK
  CL = THETA(1) * EXP(ETA(1))
  V  = THETA(2) * EXP(ETA(2))
  S1 = V
```

**MONOLIX Example:**
```
file = 'lib:oral1_1cpt_kaVCl.txt'

CL = {distribution=logNormal, typical=CL_pop, sd=omega_CL}
V = {distribution=logNormal, typical=V_pop, sd=omega_V}
```

### Two-Compartment Model

**Differential Equations:**
```
dA1/dt = -CL/V1 * A1 - Q/V1 * A1 + Q/V2 * A2
dA2/dt = Q/V1 * A1 - Q/V2 * A2
C = A1/V1
```

**Parameters:**
- `CL`: Clearance (L/h)
- `V1`: Central volume (L)
- `Q`: Inter-compartmental clearance (L/h)
- `V2`: Peripheral volume (L)
- `Ka`: Absorption rate constant (1/h) - for oral dosing

**NONMEM Example:**
```
$SUBROUTINES ADVAN3 TRANS4

$PK
  CL = THETA(1) * EXP(ETA(1))
  V1 = THETA(2) * EXP(ETA(2))
  Q  = THETA(3)
  V2 = THETA(4)
  S1 = V1
```

**MONOLIX Example:**
```
file = 'lib:oral2_2cpt_kaClV1QV2.txt'

CL = {distribution=logNormal, typical=CL_pop, sd=omega_CL}
V1 = {distribution=logNormal, typical=V1_pop, sd=omega_V1}
Q = {distribution=logNormal, typical=Q_pop, sd=0.1}
V2 = {distribution=logNormal, typical=V2_pop, sd=0.1}
```

### Three-Compartment Model

**Differential Equations:**
```
dA1/dt = -CL/V1 * A1 - Q2/V1 * A1 + Q2/V2 * A2 - Q3/V1 * A1 + Q3/V3 * A3
dA2/dt = Q2/V1 * A1 - Q2/V2 * A2
dA3/dt = Q3/V1 * A1 - Q3/V3 * A3
C = A1/V1
```

**Parameters:**
- `CL`: Clearance (L/h)
- `V1`: Central volume (L)
- `Q2`: Inter-compartmental clearance to peripheral 1 (L/h)
- `V2`: Peripheral volume 1 (L)
- `Q3`: Inter-compartmental clearance to peripheral 2 (L/h)
- `V3`: Peripheral volume 2 (L)

**NONMEM Example:**
```
$SUBROUTINES ADVAN11 TRANS4

$PK
  CL = THETA(1) * EXP(ETA(1))
  V1 = THETA(2) * EXP(ETA(2))
  Q2 = THETA(3)
  V2 = THETA(4)
  Q3 = THETA(5)
  V3 = THETA(6)
  S1 = V1
```

---

## Estimation Methods

### FOCE (First-Order Conditional Estimation)

A linearization method that approximates the likelihood by expanding around conditional estimates of random effects.

**Characteristics:**
- Deterministic algorithm
- Fast for simple models
- May struggle with non-linear models or heavy tailed distributions

**NONMEM Syntax:**
```
$ESTIMATION METHOD=1 MAXEVAL=9999 PRINT=5
```

### FOCE-I (FOCE with Interaction)

Extends FOCE by including eta-epsilon interaction in the error model.

**Characteristics:**
- More accurate for models with proportional error
- Slightly slower than standard FOCE

**NONMEM Syntax:**
```
$ESTIMATION METHOD=1 MAXEVAL=9999 INTERACTION PRINT=5
```

### SAEM (Stochastic Approximation Expectation-Maximization)

A stochastic algorithm that uses Monte Carlo sampling to handle random effects.

**Characteristics:**
- Robust for complex models
- Better handles non-linearities
- Standard in MONOLIX (default)
- Requires more iterations but typically more reliable

**NONMEM Syntax:**
```
$ESTIMATION METHOD=SAEM MAXEVAL=9999 NBURN=1000 NITER=500
```

**MONOLIX:**
SAEM is the default method in MONOLIX and doesn't require explicit specification.

---

## Data File Format

Dataset files should be in CSV format with the following columns:

### Required Columns

| Column | Description | Values |
|--------|-------------|--------|
| ID | Subject identifier | Unique integer for each subject |
| TIME | Time of observation/dose | Numeric (hours) |
| DV | Dependent variable | Concentration or response |
| AMT | Dose amount | Dose amount (0 for observations) |
| EVID | Event identifier | 0=observation, 1=dose |
| CMT | Compartment | 1=central, 2=peripheral, etc. |

### Optional Columns

| Column | Description | Example |
|--------|-------------|---------|
| WEIGHT | Body weight | kg |
| AGE | Age | years |
| SEX | Sex | 0=female, 1=male |
| RACE | Race | 1=Caucasian, 2=Black, etc. |
| CRCL | Creatinine clearance | mL/min |
| BMI | Body mass index | kg/m² |

### Example Dataset

```csv
ID,TIME,DV,AMT,EVID,CMT,WEIGHT,AGE
1,0.0,.,100,1,1,70,45
1,0.5,8.2,0,0,1,70,45
1,1.0,6.5,0,0,1,70,45
1,2.0,4.1,0,0,1,70,45
1,4.0,2.3,0,0,1,70,45
2,0.0,.,150,1,1,85,52
2,0.5,10.5,0,0,1,85,52
2,1.0,8.1,0,0,1,85,52
```

**Key Points:**
- Use `.` or blank for missing observations
- EVID=1 rows (doses) should have AMT > 0 and DV missing
- EVID=0 rows (observations) should have AMT=0 and DV present
- TIME is typically in hours
- Concentrations (DV) in appropriate units (e.g., ng/mL, µg/mL)

---

## Complete Examples

### Example 1: One-Compartment NONMEM Model

**File: `example_1comp.ctl`**

```
$PROBLEM One-compartment PK model example

$INPUT ID TIME DV AMT EVID CMT WEIGHT AGE SEX

$DATA example_dataset.csv IGNORE=C

$SUBROUTINES ADVAN1 TRANS2

$PK
  CL = THETA(1) * EXP(ETA(1))
  V  = THETA(2) * EXP(ETA(2))
  S1 = V

$THETA
  (0, 2.0)    ; CL (L/h)
  (0, 20.0)   ; V (L)

$OMEGA
  0.09        ; IIV-CL (30% CV)
  0.04        ; IIV-V (20% CV)

$SIGMA
  0.01        ; Residual error

$ERROR
  IPRED = F
  W = IPRED
  Y = IPRED + W * ERR(1)

$ESTIMATION METHOD=1 MAXEVAL=9999 PRINT=5

$COVARIANCE PRINT=E

$TABLE ID TIME IPRED CWRES NOPRINT ONEHEADER FILE=sdtab001
```

### Example 2: One-Compartment MONOLIX Model

**File: `example_1comp.mlxtran`**

```
<DATAFILE>

[FILEINFO]
file = 'example_dataset.csv'
delimiter = comma
header = {ID, TIME, DV, AMT, EVID, CMT, WEIGHT, AGE}

[CONTENT]
ID = {use=identifier}
TIME = {use=time}
DV = {use=observation, name=concentration, type=continuous}
AMT = {use=amount}
EVID = {use=eventidentifier}
CMT = {use=compartment}
WEIGHT = {use=covariate, type=continuous}
AGE = {use=covariate, type=continuous}

<MODEL>

[INDIVIDUAL]
input = {CL_pop, omega_CL, V_pop, omega_V}

DEFINITION:
CL = {distribution=logNormal, typical=CL_pop, sd=omega_CL}
V = {distribution=logNormal, typical=V_pop, sd=omega_V}

[LONGITUDINAL]
input = {a}

file = 'lib:oral1_1cpt_kaVCl.txt'

DEFINITION:
concentration = {distribution=normal, prediction=Cc, errorModel=proportional(a)}

<FIT>
data = concentration
model = concentration

<PARAMETER>
CL_pop = {value=2.0, method=MLE}
V_pop = {value=20.0, method=MLE}
omega_CL = {value=0.3, method=MLE}
omega_V = {value=0.2, method=MLE}
a = {value=0.1, method=MLE}

<MONOLIX>

[TASKS]
populationParameters()
individualParameters(method = {conditionalMean, conditionalMode})
fim(method = StochasticApproximation)
logLikelihood(method = ImportanceSampling)

[SETTINGS]
GLOBAL:
exportpath = 'example_1comp'
```

### Example 3: Two-Compartment NONMEM Model

**File: `example_2comp.ctl`**

```
$PROBLEM Two-compartment PK model example

$INPUT ID TIME DV AMT EVID CMT WEIGHT AGE

$DATA two_compartment_dataset.csv IGNORE=C

$SUBROUTINES ADVAN3 TRANS4

$PK
  CL = THETA(1) * EXP(ETA(1))
  V1 = THETA(2) * EXP(ETA(2))
  Q  = THETA(3)
  V2 = THETA(4)
  S1 = V1

$THETA
  (0, 2.0)    ; CL (L/h)
  (0, 20.0)   ; V1 (L)
  (0, 1.0)    ; Q (L/h)
  (0, 30.0)   ; V2 (L)

$OMEGA BLOCK(2)
  0.09        ; IIV-CL
  0.01 0.04   ; IIV-V1

$SIGMA
  0.01        ; Residual error

$ERROR
  IPRED = F
  W = IPRED * THETA(5) + THETA(6)
  IF (W.EQ.0) W = 1
  Y = IPRED + W * ERR(1)

$THETA
  (0, 0.1)    ; Proportional error
  (0, 0.5)    ; Additive error

$ESTIMATION METHOD=1 MAXEVAL=9999 FOCE INTERACTION PRINT=5

$COVARIANCE PRINT=E

$TABLE ID TIME IPRED CWRES NOPRINT ONEHEADER FILE=sdtab002
```

---

## Parsing Details

### NONMEM Parser Behavior

The parser (`src/parsers/nonmem.rs`) extracts:

1. **Model Type Detection:**
   - Automatically detects from ADVAN routine or compartment counts
   - Looks for keywords: `NCOMP`, `Q`, `Q2`, `Q3`, `V2`, `V3`

2. **Estimation Method:**
   - `METHOD=1` → FOCE
   - `METHOD=1 INTERACTION` → FOCE-I
   - `METHOD=SAEM` → SAEM

3. **Parameters:**
   - `$THETA`: Extracts initial values (handles both `(lower, init, upper)` and single value formats)
   - `$OMEGA`: Extracts diagonal variance values
   - `$SIGMA`: Extracts residual error variances

4. **Data File:**
   - Extracts filename from `$DATA` line

### MONOLIX Parser Behavior

The parser (`src/parsers/monolix.rs`) extracts:

1. **Model Type Detection:**
   - From library file name (e.g., `lib:oral1_1cpt` → one-compartment)
   - From keywords: `peripheral`, `k12`, `k13`, `Q`, `Q2`, `Q3`

2. **Estimation Method:**
   - Always defaults to SAEM (MONOLIX standard)

3. **Parameters:**
   - Extracts `pop_*` parameters as initial estimates
   - Extracts `omega_*` as random effect variances
   - Extracts `a` and `b` as residual error parameters

4. **Data File:**
   - Extracts from `file =` line in `<DATAFILE>` section

---

## Best Practices

### Model Building

1. **Start Simple**: Begin with a one-compartment model before adding complexity
2. **Visual Inspection**: Plot data before modeling to understand PK profile
3. **Initial Estimates**: Use literature values or simple calculations for initial estimates
4. **Covariate Analysis**: Start with base model, then add covariates systematically

### Parameter Initialization

**Good Initial Estimates:**
- `CL`: 2-10 L/h for typical adult
- `V`: 10-50 L for typical adult
- `Ka`: 0.5-2 1/h for oral absorption
- `IIV (omega)`: 0.04-0.25 (20-50% CV)
- `Residual error`: 0.01-0.1 (10-30% error)

**Avoid:**
- Initial values of 0 or negative values
- Unrealistic bounds (e.g., CL > 1000 L/h)
- Overly restrictive bounds that prevent convergence

### Error Models

**Choose Based on Data:**
- **Additive**: `Y = IPRED + a*ERR(1)` - Constant error across concentration range
- **Proportional**: `Y = IPRED * (1 + b*ERR(1))` - Error increases with concentration
- **Combined**: `Y = IPRED + (a + b*IPRED)*ERR(1)` - Most flexible, use when both patterns exist

### Estimation Method Selection

**Use FOCE when:**
- Simple models (one or two compartments)
- Limited random effects
- Need fast estimation
- Model is well-behaved

**Use SAEM when:**
- Complex models (three+ compartments)
- Many random effects
- Convergence issues with FOCE
- Non-linear models
- Sparse data

### Convergence Troubleshooting

If estimation doesn't converge:

1. **Check data**: Ensure no missing values, correct units, proper EVID coding
2. **Simplify model**: Remove random effects or compartments
3. **Adjust initial estimates**: Use values closer to expected results
4. **Increase iterations**: MAXEVAL for NONMEM
5. **Try different method**: Switch between FOCE and SAEM
6. **Check bounds**: Ensure parameter bounds are reasonable

### Validation

After estimation:

1. **Check convergence**: Look for successful minimization message
2. **Evaluate precision**: Check standard errors (RSE < 30% preferred)
3. **Visual diagnostics**: Goodness-of-fit plots, VPC
4. **Residual analysis**: Check for trends or patterns
5. **Individual fits**: Verify predictions match observations

---

## Usage Examples

### Running NONMEM Parser

```rust
use pharmacokinetics::parsers::parse_nonmem_control_stream;

let control_data = parse_nonmem_control_stream("example_1comp.ctl")?;

println!("Model: {:?}", control_data.model_type);
println!("Method: {:?}", control_data.estimation_method);
println!("Data file: {:?}", control_data.data_file);
```

### Running MONOLIX Parser

```rust
use pharmacokinetics::parsers::parse_monolix_control_stream;

let control_data = parse_monolix_control_stream("example_1comp.mlxtran")?;

println!("Model: {:?}", control_data.model_type);
println!("Method: {:?}", control_data.estimation_method);
println!("Initial estimates: {:?}", control_data.initial_estimates);
```

---

## Additional Resources

### NONMEM Documentation
- [NONMEM Users Guide](https://nonmem.iconplc.com/)
- ADVAN subroutines reference
- $ESTIMATION options guide

### MONOLIX Documentation
- [Lixoft MONOLIX](https://lixoft.com/products/monolix/)
- Model library reference
- Mlxtran language guide

### Pharmacokinetic Resources
- Gabrielsson & Weiner: "Pharmacokinetic and Pharmacodynamic Data Analysis"
- Bauer: "NONMEM Users Guide"
- FDA Guidance on Population PK

---

## Troubleshooting Common Issues

### Issue: "File not found" Error

**Solution:**
- Ensure data file path in control stream is correct
- Use relative paths from control stream location
- Check file permissions

### Issue: Parser Returns Wrong Model Type

**Solution:**
- Explicitly specify compartment count using `NCOMP=` (NONMEM)
- Use correct library file (MONOLIX)
- Verify model equations include all compartments

### Issue: Missing Initial Estimates

**Solution:**
- Ensure `$THETA`, `$OMEGA`, `$SIGMA` blocks are present (NONMEM)
- Verify `<PARAMETER>` section is complete (MONOLIX)
- Check for parsing errors in parameter values

### Issue: Estimation Doesn't Start

**Solution:**
- Verify dataset format matches control stream specifications
- Check that all required columns are present
- Ensure EVID and AMT columns are correctly coded
- Verify TIME values are numeric and properly ordered

---

## Machine Learning Models

In addition to traditional compartmental models, this system supports machine learning approaches for pharmacokinetic modeling.

### Random Forest

Random Forest is an ensemble learning method that uses multiple decision trees for regression.

**Advantages:**
- Non-parametric: No assumptions about underlying PK model structure
- Handles non-linear relationships naturally
- Robust to outliers
- Provides feature importance rankings
- Fast training and prediction

**Usage:**

```bash
phikl -d dataset.csv -m rf -o ./output
```

**Parameters:**
- Number of trees: 100 (default)
- Maximum depth: 10 (default)
- Features used: time, dose amount, time since dose, cumulative dose, subject ID

**Output:**
- Feature importance scores
- RMSE (Root Mean Square Error)
- R-squared
- MAE (Mean Absolute Error)
- Predictions by individual and time point

**When to Use:**
- Exploratory analysis
- Complex dose-response relationships
- High-dimensional covariate spaces
- Quick modeling without mechanistic assumptions
- As a benchmark against mechanistic models

### XGBoost (Extreme Gradient Boosting)

XGBoost is a gradient boosting algorithm that builds an ensemble of decision trees sequentially.

**Advantages:**
- High predictive accuracy
- Handles missing data
- Built-in regularization to prevent overfitting
- Efficient computation
- Feature importance analysis

**Usage:**

```bash
phikl -d dataset.csv -m xgb -o ./output
```

**Parameters:**
- Number of boosting rounds: 100 (default)
- Maximum tree depth: 6 (default)
- Learning rate: 0.1 (default)
- Features used: time, dose amount, time since dose, cumulative dose, subject ID

**Output:**
- Feature importance scores
- RMSE (Root Mean Square Error)
- R-squared
- MAE (Mean Absolute Error)
- Predictions by individual and time point
- Training time

**When to Use:**
- Need for high predictive accuracy
- Complex non-linear patterns
- Alternative to mechanistic models
- Model comparison and validation
- Feature selection and importance analysis

### ML Model Features

Both Random Forest and XGBoost models automatically extract and use the following features from your dataset:

1. **Time**: Observation time point
2. **Dose Amount**: Most recent dose amount
3. **Time Since Dose**: Time elapsed since last dose
4. **Cumulative Dose**: Total cumulative dose administered
5. **Subject ID**: Individual subject identifier

### Comparing ML Models with Mechanistic Models

You can compare ML models directly with compartmental models:

```bash
phikl -d dataset.csv -m 1comp -m 2comp -m rf -m xgb --compare -o ./output
```

This generates a comprehensive comparison report including:
- AIC and BIC for all models
- RMSE and R-squared comparisons
- Model rankings
- Recommendations based on fit criteria

### ML Model Interpretation

**Feature Importance:**
Feature importance scores indicate which variables contribute most to predictions. Higher scores mean greater importance.

Example interpretation:
- If "time_since_dose" has high importance: Predictions strongly depend on time after dosing
- If "cumulative_dose" has high importance: Total exposure matters more than individual doses
- If "id" has low importance: Population-level predictions are consistent across individuals

**Model Diagnostics:**
- **RMSE**: Lower values indicate better fit (same units as observations)
- **R-squared**: Proportion of variance explained (0-1, higher is better)
- **MAE**: Average absolute prediction error (same units as observations)

### Limitations of ML Models

1. **No Mechanistic Interpretation**: ML models don't provide PK parameters (CL, V, etc.)
2. **Extrapolation Risk**: Poor performance outside training data range
3. **No Covariate Effects**: Can't directly estimate covariate relationships
4. **Individual Parameters**: Don't provide individual PK parameter estimates
5. **Simulation**: Can't be used for traditional PK simulations

### Best Practices for ML Models

1. **Use as Complementary Tools**: ML models work best alongside mechanistic models
2. **Check Feature Importance**: Validate that important features make pharmacokinetic sense
3. **Cross-Validation**: Consider model performance on held-out data
4. **Compare with Mechanistic Models**: Use AIC/BIC to compare with traditional approaches
5. **Interpretation**: Focus on predictive performance rather than mechanistic understanding
6. **Data Requirements**: ML models need sufficient data (>100 observations recommended)

### Example Workflows

**Exploratory Analysis:**
```bash
# Quick exploratory modeling with ML
phikl -d new_compound.csv -m rf -m xgb -o ./exploration

# Then build mechanistic model
phikl -d new_compound.csv -m 1comp -m 2comp -e saem -o ./mechanistic

# Compare all approaches
phikl -d new_compound.csv -m 1comp -m 2comp -m rf -m xgb --compare -o ./comparison
```

**Feature Importance Analysis:**
```bash
# Use RF to identify important features
phikl -d dataset.csv -m rf -o ./features

# Check ml_summary_report.txt for feature importance rankings
# Use insights to inform covariate selection in mechanistic models
```

**Validation:**
```bash
# Use ML as benchmark
phikl -d training_data.csv -m 1comp -m xgb -o ./validation
# Compare predictive performance
```

---

## Command Line Interface

### Running Different Model Types

**Compartmental Models:**
```bash
# Single compartment with SAEM
phikl -d dataset.csv -m 1comp -e saem -o ./output

# Multiple compartments
phikl -d dataset.csv -m 1comp -m 2comp -m 3comp -o ./output

# All compartmental models
phikl -d dataset.csv -m all -o ./output
```

**Machine Learning Models:**
```bash
# Random Forest
phikl -d dataset.csv -m rf -o ./output

# XGBoost
phikl -d dataset.csv -m xgb -o ./output

# Both ML models
phikl -d dataset.csv -m rf -m xgb -o ./output
```

**Combined Analysis:**
```bash
# All models (compartmental + ML)
phikl -d dataset.csv -m all -o ./output

# Custom combination
phikl -d dataset.csv -m 1comp -m 2comp -m rf -m xgb --compare -o ./output
```

### Available Options

- `-d, --dataset`: Path to dataset CSV file
- `-C, --control-stream`: Path to NONMEM/MONOLIX control stream
- `-m, --model`: Model type(s): 1comp, 2comp, 3comp, rf (random-forest), xgb (xgboost), all
- `-e, --method`: Estimation method(s): saem, foce, foce-i, all
- `-o, --output`: Output directory
- `-i, --iterations`: Number of iterations (default: 1000)
- `-b, --burn-in`: Burn-in iterations (default: 200)
- `-c, --chains`: Number of MCMC chains (default: 4)
- `--compare`: Generate comparison report

---

## Version Compatibility

This system is compatible with:

- **NONMEM**: Versions 7.x and later
- **MONOLIX**: Versions 2018R1 and later
- **Dataset formats**: Standard NONMEM/MONOLIX CSV format
- **Machine Learning**: Random Forest and XGBoost implementations via linfa and smartcore libraries

---

## Contributing

For issues, feature requests, or contributions, please refer to the project repository.

## License

Refer to the main project license file.
