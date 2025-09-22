function [paraStru, funStru, dimC, resColl] = ...
    MStepAll(matPi, matGa, probQForZ, dimC, varStru, funStru, varargin)
%{
10-03-2023: (A) Add field "logLL" in "resColl" and input "funStru" to check the log likelihood.
            (B) Collects "matC, mainVar, & sideVar" by "paraStru" to make I/O in PGPCA consistent.

10-05-2023: (A) Add field "logLLZS" in "resColl" to check the log likelihood (the real cost!).
            (B) Change input "varStru" including "sampY, sampZ, & sampZS" to generalize the input type!

10-23-2023: (A) Add option "fitFunProbZSW" to fit new "probZ" when it's "true".
            (B) Add input "probQForZ" to include the distribution q_j(z_k).
            (C) Add option "probQForZType" to determine the way of computing "logLL & logLLProbQ".
            (D) Add output "funStru" to take care updating "probZ" in function "funProbZ".
            (E) Add & remove some fields in "resColl" to match the generalization.
    
10-26-2023: (A) Add sanity check to guarantee "logLL >= logLLProbQ (each y_j, not just the mean)".
-----------------------------------------------------------------------------------------------------------

This function computes the PGPCA M-step for general manifold model. Note that it only computes single
  M-step. The iteration is handled by the main U/I function. Comparing to standard M-step for linear 
  dynamic model (LDM), this function estimates parameters "C & sigma^2" under "various model settings".

PGPCA derivation shows that when rank(K(z_t)) = dimL < dimY, set R = (sigma^2)*I_n may cause problem since 
  sigma^2 may be too large for catching up noise in the direction unmodeled by K(z_t). One way to solve it
  is modeling those directions separately. We include these possibilities in this function.
    
An optional fitting variable is "probZ (-> from "funProbZ")" when "probQForZ ~= []". The idea is updating 
  p(z) based on the posterior distribution. This part is basically the same as the mixture PPCA. To extend 
  the application scope, I include option "fitFunProbZSW" to fit "probZ" or not.
-----------------------------------------------------------------------------------------------------------

Model:  y_t = phi(z_t) + K(z_t)*C*x_t + r_t (COV: R -> 2 forms)
    
Universal part: update "funProbZ" -> newProbZ(k) = mean(probQForZ(:,k));
    
1.  This is only valid when "probQForZ ~= []".
2.  To guarantee probability normalization "sum(probZ) = 1", we make the output of new "funStru.funProbZ" 
      be a constant column vector, the updated probability of all samples in "varStru.sampZ".

Case-dependent part: C & sigma^2.

1.  rank(K) = dimY:     No restricted direction. R = (sigma^2)*I_n.
2.  rank(K) < dimY:     Some options for the restricted directions:
    
    (A) Model them as before:       R = (sigma^2)*I_n. The modeling effect may be poor.
    (B) Model them separately:      R = blkdiag(s1^2, s2^2). The modeled/unmodeled spaces are splitted.
-----------------------------------------------------------------------------------------------------------

Standard progress:

1.  The leaking (residual) noise variance "leakVar" in the unmodeled space
    (A) If dimL = dimY, leakVar must be 0.
    (B) If dimL < dimY, leakVar = trace(matPi) - trace(matGa) >= 0.

2.  Compute the eigenvalues (descendent) with eigenvectors of "matGa".

3.  "mainVar" (= sigma^2) depends on following cases:
    (A) If dimL = dimY, then:       mainVar = [leakVar (= 0) + sum( eig(dimC+1:dimL) )] / (dimY - dimC);
                                    sideVar = [] (Because we don't model it, so we don't use "0" here!)
    (B) If dimL < dimY, then:
        [1] leakSepaSW = false:     mainVar = [leakVar + sum( eig(dimC+1:dimL) )] / (dimY - dimC);
                                    sideVar = [] (Because we still don't model it!)
        [2] leakSepaSW = true:      mainVar = sum( eig(dimC+1:dimL) / (dimL - dimC);
                                    sideVar = leakVar / (dimY - dimL);

4.  "matC" = UD where:  U includes first "dimC" eigenvectors from "matGa".
                        D = diag( sqrt[eig(1:dimC) - mainVar] ).

NOTE:
1.  If "dimL < dimY & leakSepaSW = false", then it's possible that "eig(dimC) < mainVar"! In this case, 
      we reduce dimC till the inequality is valid, and turn the warning flag on so the user knows that
      the model dimension has been changed!
-----------------------------------------------------------------------------------------------------------

Special Condition:

1.  dimC = 0:       matC    = [] (since it's not modeled.)
                    mainVar = the same as usual.
                    sideVar = the same as usual.

2.  dimC = dimL:    matC    = the same as usual.
                    mainVar = "0" except "dimL < dimY & leakSepaSW = false" -> then the same as above.
                    sideVar = the save as usual.

NOTE:
1.  This condition 0 <= dimC <= dimL always holds. Therefore, we only handle the extreme cases, which
      correspond to "isotropic model (no matC)" and "fully expanded model (no modeled noise)".
2.  If dimL < dimY, the leaking variance is always there regardless of dimC value!
-----------------------------------------------------------------------------------------------------------

NOTE:
1.  "Covering all models" is the main concern of a "general M-step", including all extreme cases above.
2.  Functions for E-step only need "paraStru". Other results are collected in "resColl".
3.  Inputs "sampY, sampZ, & funStru" are for computing "logLL" only.
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
matPi:              a symmetric matrix [dimY]. The covariance of Y in "all dimensions".
matGa:              a symmetric matrix [dimL]. The covariance of Y in "K matrix subspace".
probQForZ:          a 2D array [numY, numZ]. Each [j,k] is the probability q_j(z_k).
dimC:               an integer >= 0. The dimension of "matC", the latent state dimension.
varStru:            a struct. It collects all variable samples for log-lokelihood below:
    sampY:              a 2D array [numY, dimY]. Every row is an observation y_j.
    sampZ:              a 2D array [numZ, dimZ]. Every row is an integral point in computing p(y|z).
op  sampZS:             a 2D array [numY, dimZ]. Every row is a z_j^* sample for corresponding y_j.
funStru:            a struct. It collects all function handles in following fields:
    funProbZ:           a func. handle. The probability p(z).
    funProbYlZ:         a func. handle. The probability p(y|z).
------------------- Optional ('tag' + 'value' with (number))
"probQForZType":    (1) (probQForZType) a string "posterior / delta_func". (def: "posterior")
                        This option decides the type of q_j(z), the main point of E-step. Options are:
    "posterior":            The posterior distribution q_j(z_k) = p(z_k | y_j).
    "delta_func":           The delta function at z_j^* for each y_j by IVT and "probQForZ = []".
"fitFunProbZSW":    (1) (fitFunProbZSW) a boolean. (def: false)
                        "true" means "funStru.funProbZ" is updated based on "probQForZ". Detail is above.
"leakSepaSW":       (1) (leakSepaSW) a boolean. (def: "true" if dimL < dimY, otherwise "false".)
                        "true" means the leaking noise in the unmodeled space is handled separately.
"failMsgType":      (1) (failMsgType) a string "error / warning / none". (def: "none")
                        It decides how this function responds to user-unexpected modification, e.g., 
                          dimC is reduced. The options are:
    "error":                The command window shows an error and the code stop.
    "warning":              The command window shows a warning message.
    "none":                 Nothing.

Outputs:
------------------- Essential
paraStru:           a struct. It contains all system parameters below.
    matC:               a 2D array [dimL, dimC] / empty []. The dimensionality reduction matrix.
    mainVar:            a scalar >= 0. The average noise VAR of the "modeled space".
    sideVar:            a scalar >= 0 / empty []. The average noise VAR of the "unmodeled space".
funStru:            a struct. Besides original fields from input, the new/revised fields are:
    funProbZ:           a func. handle. It's only revised when "fitFunProbZSW = true". Detail is above.
dimC:               an integer >= 0. The "updated" dimension of "matC".
resColl:            a struct. It collects all related results in following fields:
    probQForZType:      a string. It's equal to the optional input "probQForZType".
    fitFunProbZSW:      a boolean. It's equal to the optional input "fitFunProbZSW".
    leakSepaSW:         a boolean. "true" means the leaking unmodeled noise is handled by "sideVar".
    failMsgType:        a string. It decides how this function responds to user-unexpected modification.
    logLL:              a column [numY, 1]. The log-LL from "sampY & sampZ" with new learned "paraStru".
    logLLProbQ:         a column [numY, 1]. The log-LL ELBO based on parameters from M-step learning.
    leakVar:            a scalar >= 0. The total noise variance in the unmodeled space.
    eigValGa:           a row [1, dimL]. The descendent eigenvalues of "matGa".
    eigVecGa:           a square [dimL]. Every column eigVecGa(:,j) is an eigenvector of "matGa".
    dimCOri:            an integer >= 0. The "original" dimC given by user.
    dimCEnd:            an integer >= 0. The "ending" dimC after necessary modification.
    mainVarCandi:       a row [1, dimL]. The "mainVar" under different "dimC = [0:1:dimL]".
    matCEigCandi:       a row [1, dimL]. The minimal eigenvalues criterion given "dimC = [0:1:dimL]".
    dimCUpperBd:        an integer >= 0. The upper bound of "dimC" (larger ones violate the condition.)
    matU:               a 2D array [dimL, dimC]. The matrix "U" of "matC = U*D".
    vecD:               a row [1, dimC]. The diagonal of matrix "D" of "matC = U*D".
    newProbZ:           a column [numZ, 1] / empty []. The updated "probZ" (if asked).
    warnFlag:           a struct. It collects all warning flags to highlight hidden problems:
        dimC:               "true" means "dimCOri ~= dimCEnd". The dimension is modified.
%}

%% Assign inputs.
% Extract variables for convenience.
sampY  = varStru.sampY;
sampZ  = varStru.sampZ;

% Constant.
dimY = size(matPi, 1);
dimL = size(matGa, 1);
numY = size(sampY, 1);
numZ = size(sampZ, 1);
errRatio = PGPCA.errRatio;

% Initialization.
warnFlag.dimC = false;

% Parser for optional inputs.
optionP = inputParser;
optionP.addParameter('probQForZType',   "posterior",    @(x) any(x == ["posterior", "delta_func"]));
optionP.addParameter('fitFunProbZSW',   false,          @islogical);
optionP.addParameter('leakSepaSW',      [],             @(x) isempty(x) || islogical(x));
optionP.addParameter('failMsgType',     "none",         @(x) ismember(x, ["error", "warning", "none"]));

% Parse the inputs.
optionP.parse(varargin{:});
probQForZType   = optionP.Results.probQForZType;
fitFunProbZSW   = optionP.Results.fitFunProbZSW;
leakSepaSW      = optionP.Results.leakSepaSW;
failMsgType     = optionP.Results.failMsgType;

%% Initialization & sanity check.
% Initialize: if "leakSepaSW" is not given.
if( isempty(leakSepaSW) )
    leakSepaSW = (dimL < dimY);
end
    
% Assert: matPi & matGa are symmetric matrix.
assert( issymmetric(matPi), '"matPi" must be symmetric.' );
assert( issymmetric(matGa), '"matGa" must be symmetric.' );

% Assert: 0 <= dimC <= dimL <= dimY.
assert( 0 <= dimC && dimC <= dimL && dimL <= dimY, ...
        'The order must hold: 0 <= dimC (%d) <= dimL (%d) <= dimY (%d).', dimC, dimL, dimY );

%% Prepare: matrix traces & eigen-decomposition.
% Assert: both traces (= noise variance) are >= 0.
matPiTra = trace(matPi);
matGaTra = trace(matGa);
assert( matPiTra >= 0, 'trace(matPi) (%2.2e) must be >= 0.', matPiTra );
assert( matGaTra >= 0, 'trace(matGa) (%2.2e) must be >= 0.', matGaTra );

% Compute "leakVar".
leakVar = matPiTra - matGaTra;
% The tolerance of "leakVar" w.r.t zero.
errTol = matPiTra * errRatio;
% Assert and calibrate "leakVar".
if( dimL == dimY )
    % CASE: it's full dimension, so leakVar must be 0.
    assert( abs(leakVar) < errTol, '"leakVar" (%2.2e) must be = 0 (errTol: %2.2e).', leakVar, errTol );
    leakVar = 0;
else
    % CASE: it's not full-D, so leakVar must be >= 0.
    assert( leakVar > -errTol, '"leakVar" (%2.2e) must be >= 0 (errTol: %2.2e).', leakVar, errTol );
    leakVar = max([0, leakVar]);
end

% Compute eigenvalues & eigenvectors. Then organize their forms.
[eigVecGa, eigValGa] = eig(matGa, "vector");
eigValGa = fliplr(eigValGa');
eigVecGa = fliplr(eigVecGa);

% Assert: eigenvalues of "matGa" must be >= 0.
assert( all(eigValGa >= 0), 'All values of "eigValGa" (min: %2.2e) should be >= 0.', min(eigValGa) );

%% Compute: all "mainVar" under all "dimC = [0, dimL]" (NOTE: they may contradict to "matC" condition.) 
% The numerator of "mainVar" for dimC = [0, dimL].
mainVarNume = [cumsum(eigValGa, 'reverse'), 0];

% Add "leakVar" if needed (NOTE: "leakVar" has been calibrated above).
if(~leakSepaSW)
    mainVarNume = mainVarNume + leakVar;
end

% The denominator of "mainVar". Take care when dimC = dimL = dimY (avoid divided by zero).
mainVarDeno = dimY - [0:1:dimL];
mainVarDeno(end) = max([1, mainVarDeno(end)]);

% Compare the "mainVar" candidates with corresponding eig(dimC), and find the upper bound of "dimC".
mainVarCandi = mainVarNume ./ mainVarDeno;

%% Calibrate: dimC (so we don't compute "matC & mainVar" again and again!)
matCEigCandi = [inf, eigValGa];
dimCCondiList = matCEigCandi >= mainVarCandi;
dimCUpperBd = find(dimCCondiList, 1, 'last') - 1;

% Assert: make sure all indexes lower then "dimCUpperBd" are "true". Note the shift "+1".
assert( all(dimCCondiList(1:dimCUpperBd + 1)), ...
        'All conditions lower than "dimCUpperBd" (%d) should be true!', dimCUpperBd );
    
% Finally, we decide "dimC".
dimCOri = dimC;

% Fix "dimC" if it's larger than the upper bound.
if( dimC > dimCUpperBd )
    dimC = dimCUpperBd;
    warnFlag.dimC = true;
    
    % Show messages following user requirement.
    switch failMsgType
        case "none"
            % Do nothing.            
        case "error"
            error('User-given "dimC" (%d) must be <= dimC upper bound (%d).', dimC, dimCUpperBd);
        case "warning"
            warning('User-given "dimC" (%d) must be <= dimC upper bound (%d).', dimC, dimCUpperBd);
        otherwise
            error('The given failed message type "%s" is invalid.', failMsgType);
    end
end

% Record the result.
dimCEnd = dimC;

%% Compute: mainVar, sideVar, and matC.
% Extract "mainVar" from its candidates. (NOTE: the index shifting of "dimC".)
mainVar = mainVarCandi(dimC + 1);

% Compute: sideVar.
if( dimL < dimY && leakSepaSW )
    % CASE: the unmodeled space EXISTS & the user want to model it separately.
    sideVar = leakVar / (dimY - dimL);
else
    % CASE: the unmodeled space (if it exists) is included in "mainVar" already.
    sideVar = [];
end

% Compute: matC.
if( dimC == 0 )
    % CASE: There is no "matC". All randomness is captured by "mainVar & sideVar".
    matC = [];
    matU = [];
    vecD = [];
else
    % CASE: The model includes "matC".
    matU = eigVecGa(:, 1:dimC);
    
    % Assert: all "eigValGa(1:dimC) - mainVar >= 0". Just to play safe!
    vecD = eigValGa(1:dimC) - mainVar;
    assert( all(vecD >= 0), 'All values of "matD" (min: %2.2e) must be >= 0.', min(vecD) );
    vecD = sqrt(vecD);
    
    % Finalize "matC".
    matC = matU * diag(vecD);
end

%% (If asked) Update "funStru.funProbZ".
if( fitFunProbZSW )
    % Assert: probQForZ must exist with dimension [numY, numZ].
    assert( ~isempty(probQForZ), 'To update "probZ", "probQForZ" cannot be empty!' );
    CheckMultiDimIndex(size(probQForZ),   [numY, numZ]);

    % Compute new "probZ" and update "funStru.funProbZ".
    newProbZ = mean(probQForZ, 1)';
    funStru.funProbZ = @(x) newProbZ;
else
    % This empty variable is for output field integrity.
    newProbZ = [];
end

%% Compute: logLL.
% Collect all system parameters. 
paraStru = struct('matC',       matC, ...
                  'mainVar',    mainVar, ...
                  'sideVar',    sideVar);

% logLL: it's based on "sampY & sampZ".
probYlZ = funStru.funProbYlZ(sampY, sampZ, paraStru);
probZ   = funStru.funProbZ(sampZ);
logLL   = log( probYlZ * probZ );
assert( all(IsNumericProper(logLL)), '"logLL" must be properly numeric!' );

%% Compute: logLLProbQ
% Prepare "inteZStru" depending on "probQForZType".
switch probQForZType
    case "posterior"
        % CASE: q_j(z_k) = p(z_k | y_j), the posterior distribution with "weight".
        inteZStru.zGp      = sampZ;
        inteZStru.zIndSet  = repmat([1:1:numZ], numY, 1);
        inteZStru.zProb    = probQForZ;

    case "delta_func"
        % CASE: q_j(z) = \delta(z - z_j^*) where z_j^* from IVT.
        inteZStru.zGp      = varStru.sampZS;
        inteZStru.zIndSet  = [1:1:numY]';
        inteZStru.zProb    = ones(numY, 1);
        
    otherwise
        error('The given "probQForZType" (%s) is invalid.', probQForZType);
end

% logLLProbQ: the ELBO of "logLL".
logLLProbQ = PGPCA.ELBOForEM(probQForZType, sampY, inteZStru, paraStru, funStru);

% (Only "posterior") assert: logLL >= logLLProbQ (sample-wise).
if( probQForZType == "posterior" )
    % NOTE: We put "negative tolerance" for the case that "logLL ~ logLLProbQ", which is possible from 
    %         Jensen's inequality.
    logLLDiff = logLL - logLLProbQ;
    errTol = abs(mean(logLL)) * errRatio;
    assert( all(logLLDiff >= -errTol), ...
            'Under "posterior", all "logLL - logLLProbQ" (min: %2.2e) should >= 0 (errTol: %2.2e).', ...
            min(logLLDiff), -errTol );
end

%% Prepare outputs.
resColl = struct('probQForZType',       probQForZType, ...
                 'fitFunProbZSW',       fitFunProbZSW, ...
                 'leakSepaSW',          leakSepaSW, ...
                 'failMsgType',         failMsgType, ...
                 'logLL',               logLL, ...
                 'logLLProbQ',          logLLProbQ, ...
                 'leakVar',             leakVar, ...
                 'eigValGa',            eigValGa, ...
                 'eigVecGa',            eigVecGa, ...
                 'dimCOri',             dimCOri, ...
                 'dimCEnd',             dimCEnd, ...
                 'mainVarCandi',        mainVarCandi, ...
                 'matCEigCandi',        matCEigCandi, ...
                 'dimCUpperBd',         dimCUpperBd, ...
                 'matU',                matU, ...
                 'vecD',                vecD, ...
                 'newProbZ',            newProbZ, ...
                 'warnFlag',            warnFlag);

end

