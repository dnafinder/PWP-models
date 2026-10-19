%% PWP_DEMO Prentice-Williams-Peterson models on simulated recurrent events
% This script simulates a randomized trial with recurrent events, in which
% the baseline hazard increases with each successive event (state
% dependence) and the true treatment hazard ratio is exp(-0.5) ~ 0.61.
% It then fits the Andersen-Gill model, the PWP counting-process (total
% time) model and the PWP gap-time model with PWP, and compares the
% treatment hazard ratios side by side.
% The simulation follows the MetricGate tutorial on PWP models:
% https://metricgate.com/blogs/prentice-williams-peterson-recurrent/
%
% The simulated data are left in the workspace as the table DAT
% (one row per at-risk interval: id, start, stop, event, enum, treat, age,
% gap) and, if saveCSV is true, written to pwp_demo_data.csv.
%
% See also pwp
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo.75@gmail.com
%           GitHub: https://github.com/dnafinder/PWP-models

rng(2026)

%% 1. Simulate recurrent events for n subjects
n = 250;                 % subjects
adminCensor = 24;        % months of administrative follow-up
trueBetaTreat = -0.50;   % treatment hazard ratio = exp(-0.5) ~ 0.61
betaAge = 0.01;          % log-hazard per year of age (centred at 60)
inflation = 0.20;        % baseline hazard increases ~20% per event
baseRate = 0.06;         % baseline rate of the first event (per month)
maxEvents = 20;          % safety cap
saveCSV = false;         % true: write the data to pwp_demo_data.csv

treat = double(rand(n,1) < 0.5);    % randomized treatment
age = round(60 + 8*randn(n,1));     % age at baseline

nmax = n*(maxEvents+1);
id = zeros(nmax,1); start = zeros(nmax,1); stop = zeros(nmax,1);
event = zeros(nmax,1); enum = zeros(nmax,1);
r = 0;
for i = 1:n
    t = 0;
    lp = trueBetaTreat*treat(i) + betaAge*(age(i) - 60);
    for ev = 1:maxEvents+1
        % gap-time hazard scales with event number (state dependence)
        rate = baseRate*(1 + (ev - 1)*inflation)*exp(lp);
        tnew = t - log(rand)/rate;      % exponential gap
        r = r + 1;
        id(r) = i; start(r) = t; enum(r) = ev;
        if tnew >= adminCensor || ev > maxEvents
            stop(r) = adminCensor; event(r) = 0;
            break
        end
        stop(r) = tnew; event(r) = 1;
        t = tnew;
    end
end
id = id(1:r); start = start(1:r); stop = stop(1:r);
event = event(1:r); enum = enum(1:r);
dat = table(id,start,stop,event,enum,treat(id),age(id), ...
    'VariableNames',{'id','start','stop','event','enum','treat','age'});
dat.gap = dat.stop - dat.start;
clear id start stop event enum r t tnew rate lp i ev nmax

fprintf('=== Simulated recurrent-event data ===\n')
fprintf('Subjects: %d   Rows: %d   Total events: %d\n',n,height(dat),sum(dat.event))
evPerId = accumarray(dat.id,dat.event);
fprintf('Events per subject (truncated at >=5):\n')
fprintf('%8s','0','1','2','3','4','>=5'); fprintf('\n')
fprintf('%8d',histcounts(min(evPerId,5),-0.5:1:5.5)); fprintf('\n')
if saveCSV
    writetable(dat,'pwp_demo_data.csv')
end

covs = {'treat','age'};

%% 2. Andersen-Gill: common baseline, no stratification
fitAG = pwp(dat,'Covariates',covs,'Model','ag');

%% 3. PWP-CP (total time): strata = event number, collapsed at 5
fitCP = pwp(dat,'Covariates',covs,'Model','cp','MaxStrata',5,'Plot',true);

%% 4. PWP-GT (gap time): the clock resets after each event
fitGT = pwp(dat,'Covariates',covs,'Model','gt','MaxStrata',5,'Plot',true);

%% 5. Side-by-side comparison of the treatment effect
fits = {fitAG,fitCP,fitGT};
Model = {'Andersen-Gill';'PWP-CP (total time)';'PWP-GT (gap time)'};
HR = cellfun(@(f) f.coef{'treat','HR'},fits)';
SE_robust = cellfun(@(f) f.coef{'treat','robust_se'},fits)';
z = cellfun(@(f) f.coef{'treat','z'},fits)';
p = cellfun(@(f) f.coef{'treat','p'},fits)';
cmp = table(Model,HR,SE_robust,z,p);
fprintf('=== Treatment hazard ratio across models ===\n')
fprintf('True simulated HR = %.3f\n\n',exp(trueBetaTreat))
disp(cmp)
fprintf(['PWP separates the baseline hazard by event number (see the\n' ...
    'figures: later events have higher, steeper hazards), while AG\n' ...
    'forces them to be equal.\n\n'])

%% 6. Optional: treatment effect allowed to differ by event number
% fitBy = pwp(dat,'Covariates',covs,'Model','cp','MaxStrata',5,'ByStratum','treat');

%% 7. Optional: let PWP ask the questions and choose the model
% fitAsk = pwp(dat,'Covariates',covs);
