[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/PWP-models&file=pwp.m)

# PWP-models

## 📘 Overview
pwp is a MATLAB function that fits Prentice–Williams–Peterson (PWP) Cox models for recurrent events, together with the Andersen–Gill (AG) model for comparison. Its interface and output mirror R's survival::coxph with strata() and cluster().

When a patient is hospitalized twice for heart failure, an asthmatic has six exacerbations or a device fails a dozen times, a classical Cox model is no longer adequate: repeated events of the same subject are correlated, and the risk of the next event often depends on how many events have already occurred. The PWP model is the classical conditional approach to this setting:

- the partial likelihood is stratified by event number, so that each recurrence has its own baseline hazard λ0k(t)
- a subject is at risk for the k-th event only after having experienced the (k−1)-th event
- covariate effects are shared across strata (or, optionally, allowed to differ)
- robust (sandwich) standard errors clustered on subject account for within-subject correlation

Two clocks are available:

- PWP-CP (counting process / total time): λk(t | Z) = λ0k(t)·exp(βᵀZ), t > t(k−1), with t = time since study entry
- PWP-GT (gap time): λk(s | Z) = λ0k(s)·exp(βᵀZ), s = t − t(k−1), with the clock reset to zero after every event

The Andersen–Gill model uses the same counting-process data with a single common baseline hazard: λ(t | Z) = λ0(t)·exp(βᵀZ).

If the model is not specified, pwp summarizes the data and asks a few questions to choose the appropriate one.

## ✨ Features
- Three models: PWP-CP (total time), PWP-GT (gap time) and Andersen–Gill
- Interactive model selection: when 'Model' is omitted, pwp shows rows and events per event number and asks:
  - whether the baseline risk is expected to change with event order (no → AG)
  - whether the risk is driven by time since study entry (→ PWP-CP) or by time since the previous event (→ PWP-GT)
  - from which event number high-order strata should be collapsed (a value is suggested)
  - whether some covariate effects should differ across event numbers
  - at the end, the equivalent non-interactive call is printed, for reproducibility
- Accepts:
  - a table in counting-process format (one row per at-risk interval: id, start, stop, event), as in R, or
  - separate vectors id, tstart, tstop, event and a covariate matrix X
- Event numbers (strata) and gap times are computed internally from each subject's event history, so subjects may have several rows between two events (e.g. time-varying covariates)
- Categorical, string and cellstr covariates are dummy coded (first level = reference)
- Collapsing of rare high-order event numbers into one stratum (R: pmin(enum, K))
- Stratum-specific covariate effects (R: x:strata(enum))
- Efron (default, as in R) or Breslow handling of tied event times
- Delayed entry handled through (start, stop] risk sets
- Computes:
  - coefficients, hazard ratios and (1 − α)·100% confidence intervals
  - model-based and robust standard errors clustered on subject (R: cluster(id))
  - likelihood ratio, Wald and score (log-rank) tests
  - stratum-specific baseline cumulative hazards
- Prints an R-like summary, including the equivalent coxph call
- Optional plot of the baseline cumulative hazards by event number
- Base MATLAB only: no toolbox required

## 📥 Installation
1. Download or clone the repository:
   https://github.com/dnafinder/PWP-models

2. Add the folder to your MATLAB path:
      addpath('path_to_PWP-models_folder')

3. Verify that MATLAB can see the function:
      which pwp

If MATLAB returns the path to pwp.m, the installation is successful.

## ⚙️ Requirements
- MATLAB R2016b or later (implicit expansion)
- No additional toolboxes are required

## 📈 Usage

Data must be in counting-process format, one row per at-risk interval per subject:

    %  id  start  stop  event  treat  age
    %   1    0     5.2    1      1     58
    %   1   5.2   11.7    1      1     58
    %   1  11.7   24.0    0      1     58
    %   2    0    24.0    0      0     66

PWP counting process (total time), collapsing event numbers ≥ 5:

    out = pwp(dat, 'Covariates', {'treat','age'}, 'Model', 'cp', 'MaxStrata', 5);

PWP gap time:

    out = pwp(dat, 'Covariates', {'treat','age'}, 'Model', 'gt', 'MaxStrata', 5);

Andersen–Gill model:

    out = pwp(dat, 'Covariates', {'treat','age'}, 'Model', 'ag');

Let pwp ask the questions and choose the model:

    out = pwp(dat, 'Covariates', {'treat','age'});

Treatment effect allowed to differ by event number:

    out = pwp(dat, 'Covariates', {'treat','age'}, 'Model', 'cp', ...
              'MaxStrata', 5, 'ByStratum', 'treat');

Matrix syntax, with custom variable names in the table:

    out = pwp(id, tstart, tstop, event, X, 'Model', 'cp', 'VarNames', {'treat','age'});
    out = pwp(tbl, 'ID', 'patient', 'Start', 't0', 'Stop', 't1', 'Event', 'relapse', ...
              'Covariates', {'arm'}, 'Model', 'gt');

Plot the stratum-specific baseline cumulative hazards:

    out = pwp(dat, 'Covariates', {'treat','age'}, 'Model', 'cp', 'MaxStrata', 5, 'Plot', true);

## 🔢 Inputs

pwp(tbl, 'Covariates', covs)  
pwp(tbl, 'Covariates', covs, Name, Value, ...)  
pwp(id, tstart, tstop, event, X, Name, Value, ...)

- tbl  
  - Type: table
  - Description: counting-process data, one row per at-risk interval per subject. Rows of the same subject must not overlap.

- id, tstart, tstop, event, X (matrix syntax)  
  - id: subject identifier (numeric, string, cellstr or categorical)
  - tstart, tstop: interval limits, as time since study entry (tstop > tstart ≥ 0)
  - event: 1 if the interval ends with an event, 0 if censored
  - X: n-by-p numeric matrix of covariates

- Name-Value pairs
  - Model: 'ask' (default, interactive), 'cp' (or 'tt', 'total'), 'gt' (or 'gap'), 'ag'
  - Covariates: (table syntax) names of the covariate variables
  - ID, Start, Stop, Event: (table syntax) variable names; defaults 'id', 'start', 'stop', 'event'
  - VarNames: (matrix syntax) names of the columns of X
  - MaxStrata: event numbers ≥ MaxStrata are collapsed into one stratum; default Inf
  - ByStratum: covariates whose effect may differ across strata (names, or column indices for the matrix syntax); default none
  - Ties: 'efron' (default) or 'breslow'
  - Robust: true (default) uses robust SE for z tests and confidence intervals; false uses model-based SE
  - Alpha: significance level for confidence intervals; default 0.05
  - Display: print the summary; default true
  - Plot: plot the baseline cumulative hazards; default false
  - MaxIter: maximum Newton–Raphson iterations; default 30
  - Tol: relative convergence tolerance on the log-likelihood; default 1e-9

## 📤 Outputs

    out = pwp(...);

- out.model, out.modelCode: model description and code ('cp', 'gt', 'ag')
- out.call: equivalent pwp call (useful after the interactive mode)
- out.Rcall: equivalent R survival::coxph call
- out.coef: table with coef, HR, se, robust_se, z, p and HR confidence limits
- out.beta, out.V, out.Vrobust: estimates, model-based and robust covariance matrices
- out.loglik: partial log-likelihood of the null and fitted models
- out.tests: LR, Wald and Score tests, each as [statistic df p-value]
- out.strata: subjects, rows and events per stratum
- out.basehaz: baseline cumulative hazard per stratum, at the covariate means (as R basehaz/survfit)
- out.data: analysed intervals (id, enum, stratum, start, stop, event) on the time scale of the model
- out.n, out.nsubjects, out.nevents, out.nmissing, out.iterations, out.converged

## 🧠 Interpretation
- PWP is conditional: exp(β) is the hazard ratio for the next event, given that the previous events have occurred. AG estimates a marginal rate ratio. The two models answer different questions and their hazard ratios are not expected to coincide.
- PWP-CP answers "how does the covariate change the risk along the study timeline, given the event count?"; it is suited to ageing, disease progression or cumulative exposure.
- PWP-GT answers "how does the covariate change the time to the next event after each recurrence?"; it is suited to renewal-like, rest-and-relapse processes.
- If PWP-CP and PWP-GT hazard ratios disagree noticeably, that is itself a finding: the time-to-next-event dynamics differ from the calendar-time dynamics, and a single AG hazard ratio would hide it.
- Baseline cumulative hazards that rise faster for higher event numbers indicate state dependence (each event makes the next one more likely), which AG cannot represent.

## 📌 Notes
- Always use robust standard errors (the default): without them, standard errors assume independent rows and are anticonservative. With clustered data the likelihood ratio test is reported for reference only; use the robust Wald test.
- Cap the event number: a subject with 17 events creates 17 strata, most with tiny risk sets. In the interactive mode pwp suggests the last event number with at least 10 events.
- PWP-GT requires the clock to restart after every event. pwp computes gap times from each subject's event history, not as stop − start of each row, so it remains correct when a subject has several rows between two events.
- Proportional hazards are assumed within each event-number stratum.
- If the worst-off subjects drop out (e.g. death stops further hospitalizations), censoring is informative: consider joint frailty or competing-risks models.
- Numerical results were cross-checked against an independent implementation on the bladder2 data of R's survival package (AG, PWP-CP and PWP-GT, Efron and Breslow ties).

## 🧪 Example

The demo script simulates a randomized trial with recurrent events (250 subjects, 24 months of follow-up), in which the baseline hazard increases by 20% with each successive event and the true treatment hazard ratio is exp(−0.5) ≈ 0.61. It then fits AG, PWP-CP and PWP-GT and compares the treatment hazard ratios:

    pwp_demo

The simulated data remain in the workspace as the table dat (id, start, stop, event, enum, treat, age, gap) and can be saved to pwp_demo_data.csv by setting saveCSV = true in the script.

Expected behaviour: PWP-CP and PWP-GT recover a treatment hazard ratio close to 0.61, while the PWP baseline cumulative hazards are higher and steeper for later event numbers, the pattern built into the simulation and the one AG cannot capture.

The simulation follows the MetricGate tutorial:  
https://metricgate.com/blogs/prentice-williams-peterson-recurrent/

## 🧾 Citation
If you use pwp in research, analysis, or publications, please cite:

Cardillo G. (2026). PWP-models: Prentice–Williams–Peterson Cox models for recurrent events in MATLAB.  
Available at: https://github.com/dnafinder/PWP-models

References:
- Prentice RL, Williams BJ, Peterson AV. On the regression analysis of multivariate failure time data. Biometrika 1981;68:373–379.
- Andersen PK, Gill RD. Cox's regression model for counting processes: a large sample study. Ann Stat 1982;10:1100–1120.
- Lin DY, Wei LJ. The robust inference for the Cox proportional hazards model. J Am Stat Assoc 1989;84:1074–1078.

## 👤 Author
Giuseppe Cardillo  
Email: giuseppe.cardillo.75@gmail.com  
GitHub: https://github.com/dnafinder

## 📄 License
The code is provided as-is, without any explicit warranty.  
Please refer to the repository for licensing details if a LICENSE file is present.
