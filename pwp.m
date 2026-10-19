function out = pwp(varargin)
% PWP Prentice-Williams-Peterson (and Andersen-Gill) Cox models for recurrent events
% The Prentice-Williams-Peterson (PWP) model is a conditional approach to
% recurrent-event survival analysis. The Cox partial likelihood is
% stratified by event number, so that each recurrence has its own baseline
% hazard, and a subject is at risk for the k-th event only after having
% experienced the (k-1)-th event. Two clocks are available:
%   - PWP-CP (counting process / total time): time since study entry; the
%     k-th event is modelled on the interval (t(k-1), t(k)].
%   - PWP-GT (gap time): the clock is reset to zero after every event.
% The Andersen-Gill (AG) model is also available for comparison: same
% counting-process data, one common baseline hazard (no stratification).
% Robust (sandwich) standard errors clustered on subject are computed, as
% coxph(... + cluster(id)) does in R.
%
% If the model is not specified, PWP shows a summary of the data and asks
% a few questions to choose the appropriate model (clock, collapsing of
% high-order strata, stratum-specific covariate effects).
%
% Syntax:
%   out = pwp(tbl,'Covariates',{'x1','x2',...})
%   out = pwp(tbl,'Covariates',covs,'Model',model, Name,Value, ...)
%   out = pwp(id,tstart,tstop,event,X, Name,Value, ...)
%
% Inputs:
%   tbl    - table in counting-process format: one row per at-risk
%            interval per subject (as in R: id, start, stop, event).
%   id     - subject identifier (numeric, string, cellstr or categorical)
%   tstart - start of the interval (time since study entry)
%   tstop  - end of the interval (time since study entry), tstop > tstart
%   event  - 1 if the interval ends with an event, 0 if censored
%   X      - n-by-p numeric matrix of covariates
%   Rows of the same subject must not overlap. Event numbers (strata) and
%   gap times are computed internally from the event history, so subjects
%   may have several rows between two events (time-varying covariates).
%
% Name-Value pairs:
%   Model      - 'ask' (default): interactive choice of the model
%                'cp' (or 'tt','total'): PWP counting process / total time
%                'gt' (or 'gap'): PWP gap time
%                'ag': Andersen-Gill
%   Covariates - (table input only) names of the covariate variables.
%                Categorical, string and cellstr variables are dummy coded
%                with the first level as reference.
%   ID, Start, Stop, Event - (table input only) names of the variables
%                (defaults: 'id','start','stop','event')
%   VarNames   - (matrix input only) names of the columns of X
%   MaxStrata  - event numbers >= MaxStrata are collapsed into one stratum
%                (R: pmin(enum, MaxStrata)). Default Inf (no collapsing).
%   ByStratum  - covariates whose effect is allowed to differ across
%                strata (R: x:strata(enum)). Names (or column indices for
%                matrix input). Default none.
%   Ties       - 'efron' (default, as in R) or 'breslow'
%   Robust     - true (default): robust SE clustered on id are used for
%                z tests and confidence intervals; false: model-based SE
%   Alpha      - significance level for confidence intervals (default 0.05)
%   Display    - true (default): print an R-like summary
%   Plot       - true: plot the stratum-specific baseline cumulative
%                hazards (default false)
%   MaxIter    - maximum Newton-Raphson iterations (default 30)
%   Tol        - relative convergence tolerance on log-likelihood
%                (default 1e-9)
%
% Output:
%   out.model       - model description
%   out.call        - equivalent PWP call (useful after the interactive mode)
%   out.Rcall       - equivalent R survival::coxph call
%   out.coef        - table: coef, HR, se, robust_se, z, p, HR lower/upper
%   out.beta, out.V, out.Vrobust - estimates and covariance matrices
%   out.loglik      - [null model, fitted model] partial log-likelihood
%   out.tests       - likelihood ratio, Wald and score tests [stat df p]
%   out.strata      - table of rows and events per stratum
%   out.basehaz     - baseline cumulative hazard per stratum, at the
%                     covariate means (as R basehaz/survfit default)
%   out.data        - analysed intervals (id, enum, stratum, start, stop,
%                     event) on the time scale of the model
%   out.n, out.nsubjects, out.nevents, out.nmissing, out.iterations,
%   out.converged
%
% Example:
%   % see pwp_demo.m for a complete simulated example
%   out = pwp(dat,'Covariates',{'treat','age'},'Model','cp','MaxStrata',5);
%   out = pwp(dat,'Covariates',{'treat','age'});   % interactive choice
%
% References:
%   Prentice RL, Williams BJ, Peterson AV. On the regression analysis of
%   multivariate failure time data. Biometrika 1981;68:373-379.
%   Andersen PK, Gill RD. Cox's regression model for counting processes:
%   a large sample study. Ann Stat 1982;10:1100-1120.
%   Lin DY, Wei LJ. The robust inference for the Cox proportional hazards
%   model. J Am Stat Assoc 1989;84:1074-1078.
%
% See also pwp_demo
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo.75@gmail.com
%           GitHub: https://github.com/dnafinder/PWP-models

% -------------------- Input parsing --------------------
narginchk(1,Inf)
if istable(varargin{1})
    isTbl = true;
    tbl = varargin{1};
    nv = varargin(2:end);
else
    if nargin < 5
        error('pwp:input', ...
            'Matrix syntax requires id, tstart, tstop, event and X.')
    end
    isTbl = false;
    [id,tstart,tstop,event,X] = deal(varargin{1:5});
    nv = varargin(6:end);
end

isTxt = @(x) ischar(x) || (isstring(x) && isscalar(x));
isNames = @(x) ischar(x) || isstring(x) || iscellstr(x) || isnumeric(x);
isFlag = @(x) isscalar(x) && (islogical(x) || isnumeric(x));

p = inputParser;
p.FunctionName = 'pwp';
addParameter(p,'Model','ask',isTxt);
addParameter(p,'Covariates',{},isNames);
addParameter(p,'ID','id',isTxt);
addParameter(p,'Start','start',isTxt);
addParameter(p,'Stop','stop',isTxt);
addParameter(p,'Event','event',isTxt);
addParameter(p,'VarNames',{},isNames);
addParameter(p,'MaxStrata',Inf,@(x) isnumeric(x) && isscalar(x) && x >= 1 && ...
    (isinf(x) || x == fix(x)));
addParameter(p,'ByStratum',{},isNames);
addParameter(p,'Ties','efron',@(x) any(strcmpi(x,{'efron','breslow'})));
addParameter(p,'Robust',true,isFlag);
addParameter(p,'Alpha',0.05,@(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(p,'Display',true,isFlag);
addParameter(p,'Plot',false,isFlag);
addParameter(p,'MaxIter',30,@(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p,'Tol',1e-9,@(x) isnumeric(x) && isscalar(x) && x > 0);
parse(p,nv{:});
opt = p.Results;
clear p
for f = {'Model','ID','Start','Stop','Event','Ties'}
    opt.(f{1}) = char(opt.(f{1}));
end

% -------------------- Data extraction --------------------
if isTbl
    vn = tbl.Properties.VariableNames;
    need = {opt.ID,opt.Start,opt.Stop,opt.Event};
    miss = need(~ismember(need,vn));
    if ~isempty(miss)
        error('pwp:input','Variable(s) not found in table: %s.',strjoin(miss,', '))
    end
    covs = cellstr(string(opt.Covariates));
    covs = covs(:)';
    if isempty(covs)
        error('pwp:input','Specify the covariates: pwp(tbl,''Covariates'',{...}).')
    end
    miss = covs(~ismember(covs,vn));
    if ~isempty(miss)
        error('pwp:input','Covariate(s) not found in table: %s.',strjoin(miss,', '))
    end
    id = tbl.(opt.ID);
    tstart = tbl.(opt.Start);
    tstop = tbl.(opt.Stop);
    event = tbl.(opt.Event);
    [X,names,colVar,naX] = buildDesign(tbl,covs);
    labels = struct('id',opt.ID,'start',opt.Start,'stop',opt.Stop,'event',opt.Event);
else
    if ~isnumeric(X) && ~islogical(X)
        error('pwp:input','X must be a numeric matrix.')
    end
    X = double(X);
    if isvector(X) && numel(X) == numel(tstart)
        X = X(:);
    end
    if isempty(opt.VarNames)
        names = arrayfun(@(j) sprintf('x%d',j),1:size(X,2),'UniformOutput',false);
    else
        names = cellstr(string(opt.VarNames));
        if numel(names) ~= size(X,2)
            error('pwp:input','VarNames must have one name per column of X.')
        end
    end
    names = names(:)';
    covs = names;
    colVar = 1:size(X,2);
    naX = any(isnan(X),2);
    labels = struct('id','id','start','start','stop','stop','event','event');
end

tstart = double(tstart(:));
tstop = double(tstop(:));
event = double(event(:));
id = id(:);
n0 = numel(tstart);
if numel(tstop) ~= n0 || numel(event) ~= n0 || numel(id) ~= n0 || size(X,1) ~= n0
    error('pwp:input','id, start, stop, event and covariates must have the same number of rows.')
end
if isempty(X)
    error('pwp:input','At least one covariate is required.')
end

% missing values (R: na.omit)
naId = ismissing(id);
keep = ~(isnan(tstart) | isnan(tstop) | isnan(event) | naX | naId);
nmissing = sum(~keep);
rowIdx = find(keep);    % original row numbers
id = id(keep); tstart = tstart(keep); tstop = tstop(keep);
event = event(keep); X = X(keep,:);
n = numel(tstart);
if n == 0
    error('pwp:input','No complete rows left.')
end
if any(~ismember(event,[0 1]))
    error('pwp:input','Event must be coded 0 (censored) or 1 (event).')
end
if any(tstop <= tstart)
    error('pwp:input','Each interval must have stop > start.')
end
if any(tstart < 0)
    error('pwp:input','Times must be non-negative (time since study entry).')
end

% -------------------- Event history: event number and gap time --------------------
[~,~,cid] = unique(id);
[~,ord] = sortrows([cid tstop]);
cid = cid(ord); id = id(ord); tstart = tstart(ord); tstop = tstop(ord);
event = event(ord); X = X(ord,:);

enum = zeros(n,1);   % event number for which the row is at risk
ref = zeros(n,1);    % time of the previous event (gap-time origin)
for i = 1:n
    if i == 1 || cid(i) ~= cid(i-1)
        k = 1; r = 0;
    elseif tstart(i) < tstop(i-1) - 1e-12*max(1,abs(tstop(i-1)))
        error('pwp:input','Overlapping intervals for subject at row %d.',rowIdx(ord(i)))
    end
    enum(i) = k;
    ref(i) = r;
    if event(i) == 1
        k = k + 1;
        r = tstop(i);
    end
end
nsubjects = max(cid);

% -------------------- Model choice --------------------
model = lower(char(opt.Model));
maxStrata = opt.MaxStrata;
byStratum = opt.ByStratum;
asked = strcmp(model,'ask');
if asked
    [model,maxStrata,byStratum] = askModel(enum,event,covs);
end
switch model
    case {'cp','tt','total','pwp-cp','pwp-tt'}
        model = 'cp';
    case {'gt','gap','pwp-gt'}
        model = 'gt';
    case {'ag','andersen-gill'}
        model = 'ag';
    otherwise
        error('pwp:input','Unknown model ''%s''. Use ''cp'', ''gt'', ''ag'' or ''ask''.',model)
end

% strata and time scale
if strcmp(model,'ag')
    strat = ones(n,1);
    maxStrata = Inf;
    a = tstart; b = tstop;
else
    strat = min(enum,maxStrata);
    if strcmp(model,'gt')
        a = tstart - ref; b = tstop - ref;
    else
        a = tstart; b = tstop;
    end
end
us = unique(strat);
nstrata = numel(us);

% stratum-specific covariate effects (R: x:strata(enum))
if ischar(byStratum) || isstring(byStratum)
    byStratum = cellstr(byStratum);
end
byStratum = byStratum(:)';
if isnumeric(byStratum)
    bsIdx = byStratum(:)';
    if any(bsIdx < 1 | bsIdx > numel(covs) | bsIdx ~= fix(bsIdx))
        error('pwp:input','ByStratum indices out of range.')
    end
else
    [tf,bsIdx] = ismember(byStratum,covs);
    if any(~tf)
        error('pwp:input','ByStratum covariate(s) not in model: %s.',strjoin(byStratum(~tf),', '))
    end
end
bsIdx = unique(bsIdx);
if ~isempty(bsIdx)
    if strcmp(model,'ag')
        error('pwp:input','ByStratum requires a stratified (PWP) model.')
    end
    Xn = zeros(n,0); nn = {};
    for j = 1:size(X,2)
        if ismember(colVar(j),bsIdx)
            for k = 1:nstrata
                Xn(:,end+1) = X(:,j).*(strat == us(k)); %#ok<AGROW>
                nn{end+1} = sprintf('%s:%s',names{j},stratumLabel(us(k),maxStrata)); %#ok<AGROW>
            end
        else
            Xn(:,end+1) = X(:,j); %#ok<AGROW>
            nn{end+1} = names{j}; %#ok<AGROW>
        end
    end
    X = Xn; names = nn;
end
pp = size(X,2);
cst = all(abs(X - X(1,:)) < eps(max(1,abs(X(1,:)))),1);
if any(cst)
    error('pwp:input','Constant covariate column(s): %s.',strjoin(names(cst),', '))
end

% -------------------- Risk sets --------------------
% for each stratum and each distinct event time t: rows with a < t <= b
nev = 0;
for k = 1:nstrata
    rk = strat == us(k);
    nev = nev + numel(unique(b(rk & event == 1)));
end
RS = cell(nev,1); DS = cell(nev,1); tj = zeros(nev,1); sj = zeros(nev,1);
c = 0;
for k = 1:nstrata
    rows = find(strat == us(k));
    tk = unique(b(rows(event(rows) == 1)));
    for t = tk'
        c = c + 1;
        RS{c} = rows(a(rows) < t & b(rows) >= t);
        DS{c} = rows(b(rows) == t & event(rows) == 1);
        tj(c) = t; sj(c) = us(k);
    end
end
if nev == 0
    error('pwp:input','No events in the data.')
end
isEfron = strcmpi(opt.Ties,'efron');

% -------------------- Newton-Raphson --------------------
xm = mean(X,1);
Xc = X - xm;
beta = zeros(pp,1);
[ll0,g0,H0] = coxLik(beta,Xc,RS,DS,isEfron);
ll = ll0; g = g0; H = H0;
converged = false;
for iter = 1:opt.MaxIter
    step = H\g;
    bnew = beta + step;
    [lln,gn,Hn] = coxLik(bnew,Xc,RS,DS,isEfron);
    h = 0;
    while (~isfinite(lln) || lln < ll - opt.Tol*abs(ll)) && h < 30
        step = step/2;
        bnew = beta + step;
        [lln,gn,Hn] = coxLik(bnew,Xc,RS,DS,isEfron);
        h = h + 1;
    end
    done = abs(lln - ll) <= opt.Tol*max(1,abs(lln));
    beta = bnew; ll = lln; g = gn; H = Hn;
    if done
        converged = true;
        break
    end
end
if ~converged
    warning('pwp:convergence','Newton-Raphson did not converge in %d iterations.',opt.MaxIter)
end
if rcond(H) < 1e-12
    warning('pwp:singular','Information matrix is close to singular: check covariates.')
end

% -------------------- Variances --------------------
V = inv(H);
[U,dH] = coxScore(beta,Xc,RS,DS,isEfron);
G = sparse(cid,(1:n)',1,nsubjects,n);
Dc = G*(U*V);          % dfbeta summed by subject
Vr = full(Dc'*Dc);
V = (V + V')/2; Vr = (Vr + Vr')/2;

se = sqrt(diag(V));
rse = sqrt(diag(Vr));
if opt.Robust
    use = rse; Vuse = Vr;
else
    use = se; Vuse = V;
end
z = beta./use;
pval = erfc(abs(z)/sqrt(2));
zc = sqrt(2)*erfinv(1 - opt.Alpha);
HR = exp(beta);
lo = exp(beta - zc*use);
hi = exp(beta + zc*use);

chi2p = @(x,df) gammainc(x/2,df/2,'upper');
LR = 2*(ll - ll0);
W = beta'*(Vuse\beta);
S = g0'*(H0\g0);

% -------------------- Output --------------------
switch model
    case 'cp', mdesc = 'PWP-CP (counting process / total time)';
    case 'gt', mdesc = 'PWP-GT (gap time)';
    case 'ag', mdesc = 'Andersen-Gill';
end

out.model = mdesc;
out.modelCode = model;
out.ties = lower(char(opt.Ties));
out.maxStrata = maxStrata;
out.n = n;
out.nsubjects = nsubjects;
out.nevents = sum(event);
out.nmissing = nmissing;
out.converged = converged;
out.iterations = iter;

ci = strrep(sprintf('%g',100*(1 - opt.Alpha)),'.','_');
out.coef = table(beta,HR,se,rse,z,pval,lo,hi,'RowNames',names(:), ...
    'VariableNames',{'coef','HR','se','robust_se','z','p', ...
    ['HR_lower_' ci],['HR_upper_' ci]});
out.beta = beta;
out.V = V;
out.Vrobust = Vr;
out.robust = logical(opt.Robust);
out.loglik = [ll0 ll];
out.tests.LR = [LR pp chi2p(LR,pp)];
out.tests.Wald = [W pp chi2p(W,pp)];
out.tests.Score = [S pp chi2p(S,pp)];

srows = zeros(nstrata,1); sev = zeros(nstrata,1); ssub = zeros(nstrata,1);
slab = cell(nstrata,1);
for k = 1:nstrata
    rk = strat == us(k);
    srows(k) = sum(rk); sev(k) = sum(event(rk));
    ssub(k) = numel(unique(cid(rk)));
    slab{k} = stratumLabel(us(k),maxStrata);
end
if strcmp(model,'ag')
    slab = {'all'};
end
out.strata = table(slab,ssub,srows,sev,'VariableNames',{'Stratum','Subjects','Rows','Events'});

bh_s = cell(nev,1); bh_c = zeros(nev,1);
for k = 1:nstrata
    jk = find(sj == us(k));
    [~,o2] = sort(tj(jk));
    jk = jk(o2);
    bh_c(jk) = cumsum(dH(jk));
    bh_s(jk) = {stratumLabel(us(k),maxStrata)};
end
if strcmp(model,'ag')
    bh_s(:) = {'all'};
end
out.basehaz = table(bh_s,sj,tj,bh_c,'VariableNames',{'Stratum','StratumCode','Time','CumHaz'});

out.data = table(id,enum,strat,a,b,event,'VariableNames', ...
    {'id','enum','stratum','start','stop','event'});

[out.call,out.Rcall] = buildCalls(model,maxStrata,opt.Ties,covs,bsIdx,labels,isTbl);

if opt.Display
    printSummary(out,asked)
end
if opt.Plot
    plotBasehaz(out)
end
end

% =====================================================================
function [ll,g,H] = coxLik(beta,Xc,RS,DS,isEfron)
% stratified Cox partial log-likelihood, gradient and information matrix
eta = Xc*beta;
w = exp(eta);
p = size(Xc,2);
ll = 0; g = zeros(p,1); H = zeros(p);
for j = 1:numel(RS)
    ri = RS{j}; di = DS{j}; dd = numel(di);
    wr = w(ri); Xr = Xc(ri,:);
    S0 = sum(wr); S1 = (wr'*Xr)'; S2 = Xr'*(Xr.*wr);
    ll = ll + sum(eta(di));
    g = g + sum(Xc(di,:),1)';
    if isEfron && dd > 1
        wd = w(di); Xd = Xc(di,:);
        D0 = sum(wd); D1 = (wd'*Xd)'; D2 = Xd'*(Xd.*wd);
        for r = 0:dd-1
            f = r/dd;
            den = S0 - f*D0;
            m = (S1 - f*D1)/den;
            ll = ll - log(den);
            g = g - m;
            H = H + (S2 - f*D2)/den - m*m';
        end
    else
        m = S1/S0;
        ll = ll - dd*log(S0);
        g = g - dd*m;
        H = H + dd*(S2/S0 - m*m');
    end
end
end

% =====================================================================
function [U,dH] = coxScore(beta,Xc,RS,DS,isEfron)
% score residuals (rows) and baseline hazard increments at covariate means
[n,p] = size(Xc);
w = exp(Xc*beta);
U = zeros(n,p);
dH = zeros(numel(RS),1);
for j = 1:numel(RS)
    ri = RS{j}; di = DS{j}; dd = numel(di);
    wr = w(ri); Xr = Xc(ri,:);
    S0 = sum(wr); S1 = wr'*Xr;
    if isEfron && dd > 1
        wd = w(di); Xd = Xc(di,:);
        D0 = sum(wd); D1 = wd'*Xd;
        A0 = 0; B0 = zeros(1,p); A1 = 0; B1 = zeros(1,p); mbar = zeros(1,p);
        for r = 0:dd-1
            f = r/dd;
            den = S0 - f*D0;
            m = (S1 - f*D1)/den;
            A0 = A0 + 1/den;        B0 = B0 + m/den;
            A1 = A1 + (1-f)/den;    B1 = B1 + (1-f)*m/den;
            mbar = mbar + m/dd;
        end
    else
        m = S1/S0;
        A0 = dd/S0; B0 = dd*m/S0; A1 = A0; B1 = B0; mbar = m;
    end
    dH(j) = A0;
    U(di,:) = U(di,:) + Xc(di,:) - mbar;
    isD = ismember(ri,di);
    A = repmat(A0,numel(ri),1); A(isD) = A1;
    B = repmat(B0,numel(ri),1); B(isD,:) = repmat(B1,sum(isD),1);
    U(ri,:) = U(ri,:) - wr.*(A.*Xr - B);
end
end

% =====================================================================
function [X,names,colVar,naX] = buildDesign(tbl,covs)
% numeric design matrix; categorical variables dummy coded (first level = reference)
n = height(tbl);
X = zeros(n,0); names = {}; colVar = []; naX = false(n,1);
for v = 1:numel(covs)
    x = tbl.(covs{v});
    if isnumeric(x) || islogical(x)
        x = double(x);
        if size(x,2) ~= 1
            error('pwp:input','Covariate %s must be a single column.',covs{v})
        end
        naX = naX | isnan(x);
        X(:,end+1) = x; %#ok<AGROW>
        names{end+1} = covs{v}; %#ok<AGROW>
        colVar(end+1) = v; %#ok<AGROW>
    elseif iscategorical(x) || isstring(x) || iscellstr(x) || ischar(x)
        if ~iscategorical(x)
            if ischar(x)
                x = cellstr(x);
            end
            x = categorical(x);
        end
        x = removecats(x(:));
        naX = naX | isundefined(x);
        lv = categories(x);
        if numel(lv) < 2
            error('pwp:input','Covariate %s has fewer than two levels.',covs{v})
        end
        for L = 2:numel(lv)
            X(:,end+1) = double(x == lv{L}); %#ok<AGROW>
            names{end+1} = [covs{v} lv{L}]; %#ok<AGROW>
            colVar(end+1) = v; %#ok<AGROW>
        end
    else
        error('pwp:input','Unsupported type for covariate %s.',covs{v})
    end
end
end

% =====================================================================
function [model,maxStrata,byStratum] = askModel(enum,event,covs)
% interactive choice of the model
K = max(enum);
fprintf('\n=== PWP model selection ===\n')
fprintf('Event number   Rows   Events\n')
evk = zeros(K,1);
for k = 1:K
    evk(k) = sum(event(enum == k));
    fprintf('%12d %6d %8d\n',k,sum(enum == k),evk(k))
end
fprintf('\n')

fprintf(['Q1. Is the baseline risk expected to change with event order?\n' ...
    '    (e.g. a 2nd hospitalization tends to come faster than the 1st)\n' ...
    '    y = yes, n = no (events are exchangeable, Poisson-like), u = not sure\n'])
a1 = askChoice('    Answer [y/n/u] (default u):',{'y','n','u'},'u');
if strcmp(a1,'n')
    model = 'ag';
    maxStrata = Inf;
    byStratum = {};
    fprintf('\n-> Andersen-Gill model: one common baseline hazard.\n\n')
    return
end
if strcmp(a1,'u')
    fprintf(['    PWP is chosen: it does not force a common baseline.\n' ...
        '    Tip: fit ''ag'' too and compare the results.\n'])
end

fprintf(['\nQ2. What drives the risk of the next event?\n' ...
    '    1 = time since study entry (ageing, disease progression,\n' ...
    '        calendar time, cumulative exposure)  -> PWP-CP (total time)\n' ...
    '    2 = time since the previous event (rest-and-relapse,\n' ...
    '        renewal-like process)                 -> PWP-GT (gap time)\n' ...
    '    3 = not sure\n'])
a2 = askChoice('    Answer [1/2/3] (default 3):',{'1','2','3'},'3');
switch a2
    case '1'
        model = 'cp';
    case '2'
        model = 'gt';
    case '3'
        model = 'cp';
        fprintf(['    PWP-CP is chosen. Tip: fit ''gt'' too; if the hazard ratios\n' ...
            '    disagree noticeably, that is itself a finding.\n'])
end

sug = find(evk >= 10,1,'last');
if isempty(sug)
    sug = 1;
end
fprintf(['\nQ3. High event numbers have few subjects and unstable strata.\n' ...
    '    Event numbers >= K can be collapsed into one stratum.\n' ...
    '    Suggested K = %d (last event number with at least 10 events).\n'],sug)
while true
    s = strtrim(lower(input('    K [Enter = suggested, inf = no collapsing]: ','s')));
    if isempty(s)
        maxStrata = sug; break
    elseif strcmp(s,'inf')
        maxStrata = Inf; break
    end
    v = str2double(s);
    if isfinite(v) && v >= 1 && v == fix(v)
        maxStrata = v; break
    end
    fprintf('    Please type a positive integer or inf.\n')
end

fprintf(['\nQ4. Can covariate effects differ across event numbers?\n' ...
    '    (e.g. treatment losing efficacy on later events)\n'])
a4 = askChoice('    Answer [y/n] (default n):',{'y','n'},'n');
byStratum = {};
if strcmp(a4,'y')
    fprintf('    Available covariates: %s\n',strjoin(covs,', '))
    while true
        s = input('    Covariates (comma separated): ','s');
        byStratum = strtrim(strsplit(s,','));
        byStratum = byStratum(~cellfun(@isempty,byStratum));
        bad = byStratum(~ismember(byStratum,covs));
        if isempty(bad)
            break
        end
        fprintf('    Unknown: %s\n',strjoin(bad,', '))
    end
end
fprintf('\n')
end

% =====================================================================
function s = askChoice(prompt,valid,default)
while true
    s = strtrim(lower(input([prompt ' '],'s')));
    if isempty(s)
        s = default;
    end
    if any(strcmp(s,valid))
        return
    end
    fprintf('    Please answer one of: %s\n',strjoin(valid,'/'))
end
end

% =====================================================================
function s = stratumLabel(k,maxStrata)
if isfinite(maxStrata) && k == maxStrata
    s = sprintf('enum>=%d',k);
else
    s = sprintf('enum=%d',k);
end
end

% =====================================================================
function [mcall,rcall] = buildCalls(model,maxStrata,ties,covs,bsIdx,L,isTbl)
q = @(c) strjoin(cellfun(@(x) ['''' x ''''],c,'UniformOutput',false),',');
if isTbl
    mcall = sprintf('pwp(tbl,''Covariates'',{%s},''Model'',''%s''',q(covs),model);
else
    mcall = sprintf('pwp(id,tstart,tstop,event,X,''Model'',''%s''',model);
end
if isfinite(maxStrata) && ~strcmp(model,'ag')
    mcall = sprintf('%s,''MaxStrata'',%d',mcall,maxStrata);
end
if ~isempty(bsIdx)
    mcall = sprintf('%s,''ByStratum'',{%s}',mcall,q(covs(bsIdx)));
end
if ~strcmpi(ties,'efron')
    mcall = sprintf('%s,''Ties'',''%s''',mcall,lower(ties));
end
mcall = [mcall ')'];

if isfinite(maxStrata)
    sv = sprintf('pmin(enum, %d)',maxStrata);
else
    sv = 'enum';
end
terms = covs;
for j = bsIdx(:)'
    terms{j} = sprintf('%s:strata(%s)',covs{j},sv);
end
rhs = strjoin(terms,' + ');
switch model
    case 'ag'
        rcall = sprintf('coxph(Surv(%s, %s, %s) ~ %s + cluster(%s)', ...
            L.start,L.stop,L.event,rhs,L.id);
    case 'cp'
        rcall = sprintf('coxph(Surv(%s, %s, %s) ~ %s + strata(%s) + cluster(%s)', ...
            L.start,L.stop,L.event,rhs,sv,L.id);
    case 'gt'
        rcall = sprintf(['coxph(Surv(gap_start, gap_stop, %s) ~ %s + ' ...
            'strata(%s) + cluster(%s)'],L.event,rhs,sv,L.id);
end
rcall = sprintf('%s, ties = "%s")',rcall,lower(ties));
end

% =====================================================================
function printSummary(out,asked)
fprintf('\n=== %s ===\n',out.model)
if asked
    fprintf('Equivalent call: %s\n',out.call)
end
fprintf('R equivalent:    %s\n',out.Rcall)
if strcmp(out.modelCode,'gt')
    fprintf('                 (gap_start, gap_stop = out.data.start, out.data.stop)\n')
end
fprintf('n = %d rows, %d subjects, %d events',out.n,out.nsubjects,out.nevents)
if out.nmissing > 0
    fprintf(' (%d rows deleted due to missingness)',out.nmissing)
end
fprintf('\nTies: %s; %d iterations',out.ties,out.iterations)
if ~out.converged
    fprintf(' (NOT converged)')
end
fprintf('\n\n')

C = out.coef;
nm = C.Properties.RowNames;
w = max(12,max(cellfun(@length,nm)));
fprintf(['%-' num2str(w) 's %10s %10s %10s %10s %8s %10s\n'],'', ...
    'coef','exp(coef)','se(coef)','robust se','z','Pr(>|z|)')
for j = 1:height(C)
    fprintf(['%-' num2str(w) 's %10.5f %10.5f %10.5f %10.5f %8.3f %10.4g\n'],nm{j}, ...
        C.coef(j),C.HR(j),C.se(j),C.robust_se(j),C.z(j),C.p(j))
end
vn = C.Properties.VariableNames;
fprintf(['\n%-' num2str(w) 's %10s %10s %10s %10s\n'],'','exp(coef)','exp(-coef)', ...
    strrep(vn{7},'HR_',''),strrep(vn{8},'HR_',''))
for j = 1:height(C)
    fprintf(['%-' num2str(w) 's %10.5f %10.5f %10.5f %10.5f\n'],nm{j}, ...
        C.HR(j),1/C.HR(j),C.(vn{7})(j),C.(vn{8})(j))
end
if out.robust
    fprintf('(z tests and confidence intervals use robust SE clustered on id)\n')
end

fprintf('\nStrata:\n')
S = out.strata;
fprintf('  %-10s %9s %7s %7s\n','Stratum','Subjects','Rows','Events')
for k = 1:height(S)
    fprintf('  %-10s %9d %7d %7d\n',S.Stratum{k},S.Subjects(k),S.Rows(k),S.Events(k))
end

T = out.tests;
fprintf('\nLikelihood ratio test = %.2f on %d df, p = %.4g',T.LR)
if out.robust
    fprintf('  (not valid with clustered data)')
end
if out.robust
    fprintf('\nWald test (robust)    = %.2f on %d df, p = %.4g',T.Wald)
else
    fprintf('\nWald test             = %.2f on %d df, p = %.4g',T.Wald)
end
fprintf('\nScore (logrank) test  = %.2f on %d df, p = %.4g\n\n',T.Score)
end

% =====================================================================
function plotBasehaz(out)
B = out.basehaz;
us = unique(B.StratumCode);
figure('Color','w');
hold on
cols = lines(numel(us));
lab = cell(numel(us),1);
for k = 1:numel(us)
    rk = B.StratumCode == us(k);
    t = B.Time(rk); h = B.CumHaz(rk);
    stairs([0; t],[0; h],'LineWidth',2,'Color',cols(k,:));
    lab{k} = regexprep(B.Stratum{find(rk,1)},'^enum=?','Event ');
end
hold off
grid on
if strcmp(out.modelCode,'gt')
    xlabel('Time since previous event')
else
    xlabel('Time since study entry')
end
ylabel('Baseline cumulative hazard')
title([out.model ': baseline cumulative hazards'])
if numel(us) > 1
    legend(lab,'Location','northwest','Box','off')
end
end
