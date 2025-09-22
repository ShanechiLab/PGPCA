function [matPi, matGa, probQForZ, resColl] = EStepAll(sampY, sampZ, paraStru, funStru, varargin)
%{
10-03-2023: (A) Add fields "logLL & logLLZS" in "resColl" to check the log likelihood (w/ & w/o z^*).
            (B) Update "z^*" estimation process by including "funValTrans & funValYZOpt".

10-20-2023: (A) Add option "probQForZType" to control the way of estimating q_j(z).
            (B) Add option "probZInteWei" for discretizing integration with measure "probZ" generally.
            (C) Add output "probQForZ" to include the estimated probability "q_j(z_k)" in E-step.
            (D) Add & remove some fields in "resColl" to match the generalization.

10-26-2023: (A) Add sanity check to guarantee "logLL ~ logLLProbQ (each y_j, not just the mean)".
            (B) Update option "saveType" to be more parsimonious in fields of "resColl" for saving memory.
-----------------------------------------------------------------------------------------------------------

This function computes the PGPCA E-step for general manifold model. Note that it only computes single
  E-step. The iteration is handled by the main U/I function. Clearly, the correct E-step uses posterior
  distribution, but we keep the old wrong method here for diversity and comparison.

Model:  y_t = phi(z_t) + K(z_t)*C*x_t + r_t (COV: R = (sigma^2)*I_n)

CASE 1: posterior distribution.

Distribution q_j(z) = p(z | y_j), the posterior distribution. This is the optimal EM distribution.

Step:
1.  For each y_j & z_k, compute q_j(z_k) = funProbYlZ(y_j, z_k) * funProbZ(z_k).
2.  Normalize q_j(z_k) by making sure that sum( q_j(z_k=1:numZ) ) = 1, so they are probability.
3.  Compute following symmetric matrices:

    (a) S_j(z_k)  =   [y_j - phi(z_k)] * [y_j - phi(z_k)]';
    (b) Pi(q)     =   1/numY * sum_j{sum_k{ S_j(z_k) * q_j(z_k) }};
    (c) Gamma(q)  =   1/numY * sum_j{sum_k{ K(z_k)' * S_j(z_k) * K(z_k) * q_j(z_k) }};
-----------------------------------------------------------------------------------------------------------

CASE 2: intermediate value theorem (IVT).

Distribution q_j(z) is a delta function with shifting z_j^* computed by IVT. This q_j(z) performs bad!

Step:
1.  Compute the integral of "valY(sampY) = funValYZ(sampY, sampZ) * funProbZ(sampZ)".
2.  For each y_j in sampY, find z_j^* in "funValYZ(y_j, z_j^*) = valY(y_j)" (from IVT) by "fminunc".
3.  Compute following symmetric matrices:

    (a) S_j(z_j^*)  =   [y_j - phi(z_j^*)] * [y_j - phi(z_j^*)]';
    (b) Pi(z^*)     =   1/numY * sum[ S_j(z_j^*) ];
    (c) Gamma(z^*)  =   1/numY * sum[ K(z_j^*)' * S_j(z_j^*) * K(z_j^*) ];
-----------------------------------------------------------------------------------------------------------

NOTE: To unify the "matPi & matGa" computation, all types are summarized under 3 variables:

1.  matZGp:         a 2D array [numZGp, dimZ]. Each [k, :] is a selected sampling "z" for matrices.
2.  matZIndSet:     a 2D array [numY, numZ]. Each [j, k] is a selected pair (y_j, z_k).
3.  matZProb:       a 2D array [numY, numZ]. Each [j, k] is a weight (= probability) for pair (y_j, z_k).
-----------------------------------------------------------------------------------------------------------

To speed-up E-step, all functions must accept "array-wise I/O". We list them below for reference.

Type 1: Learn z^* from IVT.

    MATLAB fminunc solves: funValYZOpt(z_j^*) = funValTrans( funValYZ(sampY, sampZ) * funProbZ(sampZ) );

    funProbZ:       [numZ, 1]    = funProbZ(sampZ)              p(z) of "sampZ".
    funValYZ:       [numY, numZ] = funValYZ(sampY, sampZ)       Element [i,j] = valYZ(y_i|z_j) for IVT.
    funValYZOpt:    [numY, numZ] = funValYZOpt(sampY, sampZ)    A version of "funValYZ" for "fminunc".
    funValTrans:    [scalar]     = funValTrans(x)               Transfer "funValYZ" output for "fminunc".

Type 2: Learn matrix "matPi & matGa" for M-step.

    funMean:        [numMatZGp, dimY] = funMean(matZGp)         Each row [j,:] = phi( zGp(j,:) ).
    funMatK:        [dimY, dimL, numMatZGp] = funMatK(matZGp)   Each slice [:,:,j] = K( zGp(j,:) ).

Type 3: Compute the final log-LL as the sanity check.

    funProbYlZ:     [numY, numZ] = funProbYlZ(y, z, varStru)    Element [i,j] = p(y_i|z_j) for log-LL.

NOTE:
1.  All output's 1st dimension is "time/sample" except function "funMatK", which is the 3rd. This is good
      for matrix-slicing in computing "matGa".
2.  "funValYZ" is not necessary to be a conditional probability p(y_i|z_j). This is good for functional
      flexibility and saving time (e.g., ignore the normalization term in probability.)
3.  K(z) must be an "orthonormal matrix" for any "z" to use the above algorithm.
-----------------------------------------------------------------------------------------------------------

NOTE:
1.  "Functional separation" is the main problem here because "a general E-step" is contradict to "a handy
      E-step". We compromise them by including "sampZ" and "multiple" functions so different models can be 
      handled by the same function handles. The price is that we need another function to take care them.
2.  For I/O flexibility, we collect all functions in input "funStru".
3.  "z_j^* initialization" in MATLAB method "fminunc" is critical, so we let user-defined as an option.
4.  M-step only needs matrices "Pi & Gamma". Other results are collected in output "resColl".
5.  Fundamentally, matrices "Pi & Gamma" from different cases are from the same integration with different
      distribution (posterior, delta function, etc.). Their exact computation above looks different due to
      discretization for numerical integration.
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
sampY:              a 2D array [numY, dimY]. Every row is an observation y_j.
sampZ:              a 2D array [numZ, dimZ]. Every row is an integral point z_k.
paraStru:           a struct. It contains all system parameters below.
    matC:               a 2D array [dimL, dimC] / empty []. The dimensionality reduction matrix.
    mainVar:            a scalar >= 0. The average noise VAR of the "modeled space".
    sideVar:            a scalar >= 0 / empty []. The average noise VAR of the "unmodeled space".
funStru:            a struct. It collects all function handles in following fields:
    ------------------- (LEARN Z^* BY IVT)
    funProbZ:           a func. handle. The probability p(z).
    funValYZ:           a func. handle. The conditional value valYZ(y|z) (may not be probability!).
    funValYZOpt:        a func. handle. It's for optimization done by "fminunc".
    funValTrans:        a func. handle. It transfers "funValYZ" for "fminunc".
    ------------------- (LEARN matPi & matGa)
    funMean:            a func. handle. The mean vector "phi(z)" for each "z" sample.
    funMatK:            a func. handle. The "orthonormal" coordinate matrix "K(z)" for each "z" sample.
    ------------------- (COMPUTE log-LL)
    funProbYlZ:         a func. handle. The probability p(y|z).
    ------------------- (FINISH)
------------------- Optional ('tag' + 'value' with (number))
"probQForZType":    (1) (probQForZType) a string "posterior / delta_func" (def: "posterior")
                        This option decides the type of q_j(z), the main point of E-step. Options are:
    "posterior":            The posterior distribution q_j(z_k) = p(z_k | y_j).
    "delta_func":           The delta function at z_j^* for each y_j by IVT.
"probZInteWei":     (1) (probZInteWei) a column [numZ, 1]. (def: ones(numZ, 1))
                        The weights for discretizing integration with measurement "probZ".
"iniZS":            (1) (iniZS) a 2D array [numY, dimZ]. (def: closest points in funValYZ(sampY, sampZ).)
                        The initialization of "fminunc" to find all z_{1:numY}^*.
"saveType":         (1) (saveType) a string "all/core" (def: "all").
                        It decides the recorded variables in "resColl" for saving memory.
    "all":                  all variables (value, initialZ, sampZS, matrix, etc.) are saved.
    "core":                 only the compressed information (e.g., logLL) and reasonable size variables 
                              are saved. Other fields are tagged "op" in the output list. They are mainly: 
    -> array dim > 2:           Mostly, tensor arrays are huge. 
    -> [numY * /big dim/]:      Since "numY" can be 10^4 ~ 10^5, those arrays are usually huge.

Outputs:
------------------- Essential
matPi:              a symmetric matrix [dimY]. The covariance of Y in "all dimensions".
matGa:              a symmetric matrix [dimL]. The covariance of Y in "K matrix subspace".
probQForZ:          a 2D array [numY, numZ] / empty []. Element [j, k] is q_j(z_k). Empty means N/A.
resColl:            a struct. It collects all related results in following fields:
    saveType:           a string. The selected option for saving variables.
    logLL:              a column [numY, 1]. The log-LL from "sampY & sampZ" (no E-step learning).
    logLLProbQ:         a column [numY, 1]. The log-LL ELBO based on q_j(z) from E-step learning.
    probZ:              a column [numZ, 1]. The output of "funProbZ(sampZ)".
op  probYlZ:            a 2D array [numY, numZ]. The output of "funProbYlZ(sampY, sampZ, paraStru)".
    probQForZType:      a string. It's equal to the optional input "probQForZType".
    probZInteWei:       a column [numZ, 1]. It's equal to the optional input "probZInteWei".
    iniZS:              a 2D array [numY, dimZ]. Every row is an initialization of "fminunc" for IVT.
    matZGp:             a 2D array [numMatZGp, dimZ]. Each [j,:] is a sample in group "z".
op  matZIndSet:         a 2D array [numY, numMatZInd]. Each [j,k] is the index of "matZGp" for q_j(z_k).
op  matZProb:           a 2D array [numY, numMatZInd]. Each [j,k] is the probability q_j(z_k).
    sampZSErrVal:       a column [numY, 1]. The residual error from optimization "fminunc".
    sampMeanGp:         a 2D array [numY, dimY]. The output of "funMean(matZGp)".
op  sampMatKGp:         a 3D array [dimY, dimL, numY]. The output of "funMatK(matZGp)".
%}

%% Assign inputs.
% Constant.
[numY, dimY] = size(sampY);
[numZ, dimZ] = size(sampZ);
errRatio = PGPCA.errRatio;

% Parser for optional inputs.
optionP = inputParser;
optionP.addParameter('probQForZType',   "posterior",    @(x) any(x == ["posterior", "delta_func"]));
optionP.addParameter('probZInteWei',    ones(numZ, 1),  @(x) iscolumn(x) && all(x > 0));
optionP.addParameter('iniZS',           [],             @(x) isempty(x) || all(size(x) == [numY, dimZ]));
optionP.addParameter('saveType',        "all",          @(x) any(x == ["all", "core"]));

% Parse the inputs.
optionP.parse(varargin{:});
probQForZType = optionP.Results.probQForZType;
probZInteWei  = optionP.Results.probZInteWei;
iniZS         = optionP.Results.iniZS;
saveType      = optionP.Results.saveType;

%% Initialize output "resColl".
% NOTE: Since not all variables exist in all "probQForZType" cases, I initialize it here and collect the 
%         avaiable ones in follow-up sessions.
resColl = struct('saveType',        saveType, ...
                 'logLL',           [], ...
                 'logLLProbQ',      [], ...
                 'probZ',           [], ...
                 'probYlZ',         [], ...
                 'probQForZType',   probQForZType, ...
                 'probZInteWei',    probZInteWei, ...
                 'iniZS',           [], ...
                 'matZGp',          [], ...
                 'matZIndSet',      [], ...
                 'matZProb',        [], ...
                 'sampZSErrVal',    [], ...
                 'sampMeanGp',      [], ...
                 'sampMatKGp',      []);

%% Prepare p(z) & p(y|z).
probZ   = funStru.funProbZ(sampZ);
probYlZ = funStru.funProbYlZ(sampY, sampZ, paraStru);

% Assert: their dimensions must be correct.
CheckMultiDimIndex(size(probZ),   [numZ, 1]);
CheckMultiDimIndex(size(probYlZ), [numY, numZ]);

%% Compute: estimated distribution q_j(z).
switch probQForZType
    case "posterior"
        % CASE: q_j(z_k) = p(z_k | y_j), the posterior distribution with "weight".
        probQForZ = probYlZ .* (probZ .* probZInteWei)';
        % Normalization: sum( q_j(z_k=1:numZ) ) = 1.
        probQForZ = probQForZ ./ sum(probQForZ, 2);
        
        % Prepare summary variables.
        matZGp      = sampZ;
        matZIndSet  = repmat([1:1:numZ], numY, 1);
        matZProb    = probQForZ;
        
    case "delta_func"
        % CASE: q_j(z) = \delta(z - z_j^*) where z_j^* from IVT.
        [sampZS, resColl.sampZSErrVal, iniZS] = IVTEstimateSampZS(sampY, sampZ, iniZS, probZ, funStru);
        
        % Prepare summary variables.
        matZGp      = sampZS;
        matZIndSet  = [1:1:numY]';
        matZProb    = ones(numY, 1);
        
    otherwise
        error('The given "probQForZType" (%s) is invalid.', probQForZType);
end

%% Compute: All pairs of "y_j - phi(z_k)", which is "sampDiff".
% All "phi(matZGp)" for further manipulation.
sampMeanGp = funStru.funMean(matZGp);

% Assert: the dimension of each "phi(matZGp) == dimY".
dimSampMean = size(sampMeanGp, 2);
assert( dimSampMean == dimY, ...
        'The dimension of "sampMean" (%d) must = "dimY" (%d).', dimSampMean, dimY );

% Compute: 3D matrix "sampDiff = [dimY, numMatZInd, numY]". "numY" is the last index for slicing!
numMatZInd = size(matZIndSet, 2);
sampDiff = zeros(dimY, numMatZInd, numY);
for(j = 1:1:numY)
    sampDiff(:,:,j) = (sampY(j,:) - sampMeanGp(matZIndSet(j,:), :))';
end

%% Compute: matPi.
% Compute "sum_j{sum_k{ S_j(z_k) * q_j(z_k) }}" then take average over "j".
matPi = zeros(dimY, dimY);
for(j = 1:1:numY)
    matPi = matPi + (sampDiff(:,:,j) .* matZProb(j,:)) * sampDiff(:,:,j)';
end
matPi = matPi / numY;

% Make sure it's symmetric.
matPi = (matPi + matPi') / 2;

% Assert: make sure that "matPi" is PSD.
PDTrans(matPi, "PSD", "none", "negaTol", false, [], []);

%% Compute: matGa. 
% All "K(matZGp)" for further manipulation.
sampMatKGp = funStru.funMatK(matZGp);

% Sanity check: the dimension (except the 2nd index) must be correct.
numMatZGp = size(matZGp, 1);
CheckDimAndTimeIndex( sampMatKGp, [1,3], [dimY, numMatZGp], [], [] );

% Compute "sum_j{sum_k{ K(z_k)' * S_j(z_k) * K(z_k) * q_j(z_k) }}" and take average over "j".
dimL = size(sampMatKGp, 2);
matGa = zeros(dimL, dimL);
for(j = 1:1:numY)
    % The selected index of matZ-related.
    matZInd = matZIndSet(j,:);
    
    % First, compute "sampMatKGp(:,:,matZInd)' * sampDiff(:,:,j)" (not K'*S*K) for saving memory.
    sampKDiff = pagemtimes(sampMatKGp(:,:,matZInd),             'transpose', ...
                           permute(sampDiff(:,:,j), [1,3,2]),   'none');
    sampKDiff = permute(sampKDiff, [1,3,2]);
    
    % Add "sum_k{ K(z_k)' * S_j(z_k) * K(z_k) * q_j(z_k) }" to "matGa".
    matGa = matGa + (sampKDiff .* matZProb(j,:)) * sampKDiff';
end
matGa = matGa / numY;

% Make sure it's symmetric.
matGa = (matGa + matGa') / 2;

% Assert: make sure that "matGa" is PSD.
PDTrans(matGa, "PSD", "none", "negaTol", false, [], []);

%% Compute: logLL & logLLZS.
% logLL: it's based on "sampY & sampZ".
logLL = log( probYlZ * probZ );
assert( all(IsNumericProper(logLL)), '"logLL" must be properly numeric!' );

% logLLProbQ: the ELBO of "logLL".
inteZStru = struct('zGp',       matZGp, ...
                   'zIndSet',   matZIndSet, ...
                   'zProb',     matZProb);
logLLProbQ = PGPCA.ELBOForEM(probQForZType, sampY, inteZStru, paraStru, funStru);

% Assert: logLL ~ logLLProbQ (sample-wise). We use "~" not ">=" due to numerical discretizing error.
if( any( probQForZType == ["posterior", "delta_func"] ) )
    logLLDiffAbs = abs(logLL - logLLProbQ);
    errTol = abs(mean(logLL)) * errRatio;
    assert( all(logLLDiffAbs <= errTol), ...
            'All "|logLL - logLLProbQ|" (max: %2.2e) must = 0 (errTol: %2.2e).', ...
            max(logLLDiffAbs), errTol );
end

%% Prepare outputs.
% Collect the variables.
resColl.logLL       = logLL;
resColl.logLLProbQ  = logLLProbQ;
resColl.probZ       = probZ;
resColl.iniZS       = iniZS;
resColl.matZGp      = matZGp;
resColl.sampMeanGp  = sampMeanGp;

% Save "large" variable if the user asks.
if(saveType == "all")
    resColl.probYlZ    = probYlZ;
    resColl.matZIndSet = matZIndSet;
    resColl.matZProb   = matZProb;
    resColl.sampMatKGp = sampMatKGp;
end

end

function [sampZS, sampZSErrVal, iniZS] = IVTEstimateSampZS(sampY, sampZ, iniZS, probZ, funStru)
%{
This auxiliary function computes z_j^* for each y_j by IVT numerically. It includes following steps:

1.  Compute "valYZ & valY" for IVT reference.
2.  Set initial value "iniZS" for numerical computation.
3.  Compute z_j^* by IVT.

Inputs:
------------------- Essential
sampY:              a 2D array [numY, dimY]. Every row is an observation y_j.
sampZ:              a 2D array [numZ, dimZ]. Every row is an integral point z_k.
iniZS:              a 2D array [numY, dimZ] / empty []. The initial value of "sampZS" for IVT.
probZ:              a column [numZ, 1]. The probability of "sampZ".
funStru:            a struct. It collects all function handles in following fields:
    funValYZ:           a func. handle. The conditional value valYZ(y|z) (may not be probability!).
    funValYZOpt:        a func. handle. It's for optimization done by "fminunc".
    funValTrans:        a func. handle. It transfers "funValYZ" for "fminunc".

Outputs:
------------------- Essential
sampZS:             a 2D array [numY, dimZ]. The collection of all z_{1:numY}^* based on IVT.
sampZSErrVal:       a column [numY, 1]. The residual error from optimization "fminunc".
iniZS:              a 2D array [numY, dimZ]. Every row is an initialization of "fminunc".
%}

%% Constants.
[numY, ~   ] = size(sampY);
[numZ, dimZ] = size(sampZ);

%% Compute "valYZ & valY".
valYZ = funStru.funValYZ(sampY, sampZ);
valY = valYZ * probZ;

% Sanity check: the dimensions of above variables must be correct.
CheckMultiDimIndex(size(valYZ), [numY, numZ]);

% Since the values are transferred by "funValTrans" in "fminunc", the initial values should follow it.
valYZ = funStru.funValTrans(valYZ);
valY  = funStru.funValTrans(valY);

%% Initialize "iniZS".
if( isempty(iniZS) )
    % Find the closest point in "valYZ" w.r.t "valY".
    [~, minInd] = min(abs(valYZ - valY), [], 2);
    % Select "iniZS" from "sampZ".
    iniZS = sampZ(minInd, :);
end

% Sanity check: the dimension should be matched.
CheckMultiDimIndex(size(iniZS), [numY, dimZ]);

%% Compute: sampZS and corresponding error.
% Initialize variables.
sampZS = zeros(numY, dimZ);
sampZSErrVal = zeros(numY, 1);

% Solve the optimization one by one.
for(j = 1:1:numY)
    % Print learning progress in MATLAB command window.
    if( mod(j, 100) == 0 )
        fprintf('E-step computes "sampZS" by IVT for-loop index: %d\n', j);
    end
    
    % Select parameters.
    seleY = sampY(j,:);
    seleValY = valY(j);
    
    % Define the function for minimization.
    targFun = @(x) abs(funStru.funValYZOpt(seleY, x) - seleValY);
    
    % The options for "fminunc" (I use "quasi-newton" since the derivative is very complex generally!).
    optionSet = optimoptions(@fminunc, 'Algorithm',    'quasi-newton', ...
                                       'GradObj',      'off', ...
                                       'Display',      'notify', ...
                                       'MaxFunEvals',  10000, ...
                                       'MaxIter',      10000);
                                 
    % find the optimal "t".
    [sampZS(j,:), sampZSErrVal(j)] = fminunc(targFun, iniZS(j,:), optionSet);
end

end

