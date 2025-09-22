function [dimC, paraStru, funStru, resEStep, resMStep, resLogLL] = ...
    EMAlgo(sampY, sampZ, dimC, paraStru, funStru, varargin)
%{
09-25-2023: IMPORTANT: THIS IS JUST A TEMPORARY FUNCTION FOR HAVING THE PRETEST RESULT FASTER. MANY OTHER
              GENERALIZATIONS AND OPTIMIZATIONS ARE NEEDED!

10-04-2023: (A) Update coding and I/O following new "EStepAll & MStepAll" framework.
            (B) Include "PGPCA.ParaFun2Model" for new func. handles in "EStepAll & MStepAll".
    
10-23-2023: (A) Add session "Print PGPCA Setting" to avoid accidental mistake.
            (B) Add options "probQForZType, probZInteWei, & fitFunProbZSW" based on "EStepAll & MStepAll".
            (C) Add sanity check to make sure "probQForZType & fitFunProbZSW" are consistent.
            (D) Add & remove some fields in "resEStep, resMStep, & resLogLL" due to "EStepAll & MStepAll".
    
10-27-2023: (A) Add sanity check to guarantee "logLLProbQ" from "posterior E-step & M-step" are in order.
            (B) Update "mainVar" null initialization to be more precise by considering "phi(z_t)" effect.
-----------------------------------------------------------------------------------------------------------
 
This function is the complete EM algorithm of PGPCA model. The generalizability is the key, so model-based
  inputs are function handles for dealing different object types automatically.

Model:      y_t = phi(z_t) + K(z_t)*C*x_t + r_t [COV: R -> Ku(z_t) * diag(mainVar, sideVar) * Ku(z_t)']

Parameter:  (1) matC:       The critical subspace in K(z_t) for maximizing LL.
            (2) mainVar:    The main noise variance (those diagonal terms) of r_t.
            (3) sideVar:    (optional) The complementary noise variance for the unmodeled space.

Optional:   (1) funProbZ:   It can be updated when "fitFunProbZSW = true" & "probQForZType = posterior".

NOTE:
1.  Some "fields = []" in input "paraStru" means default initialization. I put them in essential rather
      optional inputs so the I/O interface is more concise.
-----------------------------------------------------------------------------------------------------------

(follow "EStepAll") all functions must accept "array-wise I/O". We list them below for reference.

funProbZ:   [numZ, 1]    = funProbZ(sampZ).         The output is p(z).
funMean:    [numZ, dimY] = funMean(sampZ).          Each output row = phi(z_j) [the mean in model].
funMatK:    [dimY, dimL, numZ] = funMatK(sampZ).    Each slice [:,:,j] is a matrix K(z_j) in model.
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
sampY:              a 2D array [numY, dimY]. Every row is an observation y_j.
sampZ:              a 2D array [numZ, dimZ]. Every row is an integral point in computing "valY(sampY)".
dimC:               an integer >= 0. The dimension of "matC", the latent state dimension.
paraStru:           a struct. It contains all (learned) system parameters below.
    matC:               a 2D array [dimL, dimC] / empty []. The dimensionality reduction matrix.
    mainVar:            a scalar >= 0 / empty []. The average noise VAR of the "modeled space".
    sideVar:            a scalar >= 0 / empty []. The average noise VAR of the "unmodeled space".
funStru:            a struct. It collects all function handles in following fields:
    funProbZ:           a func. handle. The probability p(z).
    funMean:            a func. handle. The mean vector "phi(z)" for each "z" sample.
    funMatK:            a func. handle. The "orthonormal" coordinate matrix "K(z)" for each "z" sample.
------------------- Optional ('tag' + 'value' with (number))
"iterLimit":        (1) (iterLimit) a +integer (def: 100).
                        The maximal #iterations for IVT. It's only used for the "errRatio" case.
"probQForZType":    (1) (probQForZType) a string "posterior / delta_func" (def: "posterior")
                        This option decides the type of q_j(z), the main point of E-step. Options are:
    "posterior":            The posterior distribution q_j(z_k) = p(z_k | y_j).
    "delta_func":           The delta function at z_j^* for each y_j by IVT.
"probZInteWei":     (1) (probZInteWei) a column [numZ, 1]. (def: ones(numZ, 1))
                        The weights for discretizing integration with measurement "probZ".
"fitFunProbZSW":    (1) (fitFunProbZSW) a boolean. (def: false)
                        "true" means "funStru.funProbZ" is updated based on "probQForZ". Detail is above.
"iniZS":            (1) (iniZS) a 2D array [numY, dimZ]. (def: [])
                        The initialization of "fminunc" to find all z_{1:numY}^* in E-step.
"leakSepaSW":       (1) (leakSepaSW) a boolean. (def: []).
                        "true" means the leaking noise in the unmodeled space is handled separately.
"failMsgType":      (1) (failMsgType) a string "error / warning / none". (def: "none")
                        It decides how this function responds to user-unexpected modification, e.g., 
                          dimC is reduced. The options are:
    "error":                The command window shows an error and the code stop.
    "warning":              The command window shows a warning message.
    "none":                 Nothing.
"saveType":         (1) (saveType) a string "all/core" (def: "all").
                        It decides the recorded variables in "resColl" for saving memory.
    "all":                  all variables (value, initialZ, sampZS, matrix, etc.) are saved.
    "core":                 only the compressed information (e.g., logLL) and reasonable size variables 
                              are saved. Other fields are tagged "op" in the output list. They are mainly: 
    -> array dim > 2:           Mostly, tensor arrays are huge. 
    -> [numY * /big dim/]:      Since "numY" can be 10^4 ~ 10^5, those arrays are usually huge.
-----------------------------------------------------------------------------------------------------------

Outputs:
------------------- Essential
dimC:               an integer >= 0. The "updated" dimension of "matC".
paraStru:           a struct. It contains all parameters after learning below.
    matC:               a 2D array [dimL, dimC] / empty []. The dimensionality reduction matrix.
    mainVar:            a scalar >= 0. The average noise VAR of the "modeled space".
    sideVar:            a scalar >= 0 / empty []. The average noise VAR of the "unmodeled space".
funStru:            a struct. BESIDES ORIGINAL FIELDS, we update/add following fields:
    funValYZ            a func. handle. The function derived from p(y|z) for E-step z^* estimation.
    funValYZOpt         a func. handle. The function for MATLAB "fminunc" to learn z^*.
    funValTrans         a func. handle. The function transfering "funValYZ" to speed up "fminunc".
    funProbYlZ          a func. handle. The conditional probability p(y|z).
    funProbZ            a func. handle. The probability p(z).
------------------- (E-step info.)
resEStep:           a struct row [1, numIter]. Each element records the E-step values in fields below:
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
    matPi:              a symmetric matrix [dimY]. The covariance of Y in "all dimensions".
    matGa:              a symmetric matrix [dimL]. The covariance of Y in "K matrix subspace".
op  probQForZ:          a 2D array [numY, numZ] / empty []. Element [j, k] is q_j(z_k). Empty means N/A.
------------------- (M-step info.)
resMStep:           a struct row [1, numIter]. Each element records the M-step values in fields below:
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
    paraStru:           a struct. It records the optimally learned system parameters:
        matC:               a 2D array [dimL, dimC] / empty []. The optimal D-reduction matrix.
        mainVar:            a scalar >= 0. The average noise VAR of the "modeled space".
        sideVar:            a scalar >= 0 / empty []. The average noise VAR of the "unmodeled space".
------------------- (logLL info.)
resLogLL:           a struct row [1, numIter]. Each element records final "logLL" after an EM iteration:
    probMode:           a string. The main mode of the probability computation.
    probValForm:        a string. The sub-mode of the probability computation.
    probVal:            a 2D array [numY, numZ] / a column [numY, 1]. It's "logLL" in this function.
    mvnMean:            a 2D array [dimY, numZ]. Each column [:,j] is "phi(z)" at matching "sampZ".
    mvnCov:             a 3D array [dimY, dimY, numZ]. Each slice [:,:,j] is COV of "y_t" at "sampZ".
op  covResColl:         a struct. The info. in computing "mvnCov" by "PGPCA.ModelCov".
op  mvnResColl:         a struct. The info. in computing MVNPDF-related values.
%}

%% Assign inputs.
% Constant - from inputs.
[numY, dimY] = size(sampY);
[numZ, dimZ] = size(sampZ);
errRatio = PGPCA.errRatio;

% Parser for optional inputs.
optionP = inputParser;
optionP.addParameter('iterLimit',       100,            @(x) x > 0 && mod(x,1) == 0)
optionP.addParameter('probQForZType',   "posterior",    @(x) any(x == ["posterior", "delta_func"]));
optionP.addParameter('probZInteWei',    ones(numZ, 1),  @(x) iscolumn(x) && all(x > 0));
optionP.addParameter('fitFunProbZSW',   false,          @islogical);
optionP.addParameter('iniZS',           [],             @(x) isempty(x) || all(size(x) == [numY, dimZ]));
optionP.addParameter('leakSepaSW',      [],             @(x) isempty(x) || islogical(x));
optionP.addParameter('failMsgType',     "none",         @(x) ismember(x, ["error", "warning", "none"]));
optionP.addParameter('saveType',        "all",          @(x) any(x == ["all", "core"]));

% Parse the inputs.
optionP.parse(varargin{:});
iterLimit       = optionP.Results.iterLimit;
probQForZType   = optionP.Results.probQForZType;
probZInteWei    = optionP.Results.probZInteWei;
fitFunProbZSW   = optionP.Results.fitFunProbZSW;
iniZS           = optionP.Results.iniZS;
leakSepaSW      = optionP.Results.leakSepaSW;
failMsgType     = optionP.Results.failMsgType;
saveType        = optionP.Results.saveType;

%% Initialization: matC.
% Assert: "dimC" is an integer >= 0.
assert( dimC >= 0 && mod(dimC, 1) == 0, '"dimC" (%2.2e) must be an integer >= 0.', dimC );

% Constant: "dimL"
dimL = size(funStru.funMatK(sampZ(1,:)), 2);

% Extract system parameters for coding conciseness.
matC = paraStru.matC;

% Initialize: "matC".
if( isempty(matC) )
    % NOTE: this statement includes case "dimC = 0".
    matC = zeros(dimL, dimC);
end
% Assert: make sure the dimension is correct.
dimMatC = size(matC);
assert( all(dimMatC == [dimL, dimC]), ...
        '"matC" Dim. [%d, %d] must be [dimL, dimC] = [%d, %d].', dimMatC(1), dimMatC(2), dimL, dimC );

% Initialize: if "leakSepaSW" is not given.
if( isempty(leakSepaSW) )
    leakSepaSW = (dimL < dimY);
end  
    
%% Initialization: estimate a reasonable initial value of "mainVar & sideVar".
%{
Intuitively, "mainVar" captures the COV of y_t "minus phi(z_t)", so we need to estimate each term. Since 
  "paraStru" is not complete at the beginning, we solve this problem in following method.

1.  Create a temporary paraStru, "tempParaStru".
2.  Compute p(z) & p(y|z) by "tempParaStru & tempFunStru".
3.  Estimate z_t = maximal p(z|y_t).
4.  Estimate "mainVar & sideVar" in each dimension separately.
%}  

% Estimate "mainVar" first time and create "tempParaStru".
tempParaStru = struct('matC',       [], ...
                      'mainVar',    trace(sampY' * sampY) / (numY * dimY), ...
                      'sideVar',    []);

% Create "tempFunStru" and then compute p(z) & p(y|z).
[~, tempFunStru] = PGPCA.ParaFun2Model(dimY, dimL, tempParaStru, funStru, "probYlZType", "mvnpdf");
probZ   = tempFunStru.funProbZ(sampZ);
probYlZ = tempFunStru.funProbYlZ(sampY, sampZ, tempParaStru);

% Estimate "z_t = maximal p(z|y_t)".
probZlY = probYlZ .* (probZ .* probZInteWei)';
probZlY = probZlY ./ sum(probZlY, 2);
[~, indZlY] = max(probZlY, [], 2);
sampZlYPhi = tempFunStru.funMean(sampZ(indZlY, :));
   
% Estimate a reasonable initializing value for "mainVar & sideVar".
iniVarVal = sum((sampY - sampZlYPhi).^2, "all") / (numY * dimY);

%% Initialization: mainVar & sideVar.
% Extract system parameters for coding conciseness.
mainVar = paraStru.mainVar;
sideVar = paraStru.sideVar;

% Initialize: "mainVar".
if( isempty(mainVar) )
    % For simplicity, we use the average variance of "sampY" with "iniVarRatio" to control magnitude.
    mainVar = iniVarVal;
end
% Assert: "mainVar" is a nonnegative scalar.
assert( mainVar >= 0, '"mainVar" (%2.2e) must be >= 0.', mainVar );

% Initialize: "sideVar".
if( ~leakSepaSW )
    % CASE: There is no space partition, so "sideVar = []".
    sideVar = [];
else
    % CASE: Space partition exists.
    if( isempty(sideVar) )
        sideVar = iniVarVal;
    end
    % Assert: "sideVar" is a nonnegative scalar.
    assert( sideVar >= 0, '"sideVar" (%2.2e) must be >= 0.', sideVar );
end
    
%% Prepare function "funValYZ, funValYZOpt, funValTrans, & funProbYlZ".
% Collect variables.
paraStru.matC    = matC;
paraStru.mainVar = mainVar;
paraStru.sideVar = sideVar;

% Create all func. handles.
[~, funStru] = PGPCA.ParaFun2Model(dimY, dimL, paraStru, funStru, "probYlZType", "mvnpdf");

%% Assert setting consistency & print it in command window.
% Check: "fitFunProbZSW = true" only if "probQForZType = posterior".
assert( fitFunProbZSW == false || probQForZType == "posterior", ...
        '"fitFunProbZSW" (%d) is "true" only if "probQForZType" (%s) is "posterior".', ...
        fitFunProbZSW, probQForZType );

% Print the current setting.
dispSettingStru = struct('numY',                numY, ...
                         'dimY',                dimY, ...
                         'numZ',                numZ, ...
                         'dimZ',                dimZ, ...
                         'dimL',                dimL, ...
                         'dimC',                dimC, ...
                         'matC',                matC, ...
                         'mainVar',             mainVar, ...
                         'sideVar',             sideVar, ...
                         'iterLimit',           iterLimit, ...
                         'probQForZType',       probQForZType, ...
                         'probZInteWei',        probZInteWei, ... 
                         'fitFunProbZSW',       fitFunProbZSW, ...
                         'iniZS',               iniZS, ...
                         'leakSepaSW',          leakSepaSW, ...
                         'failMsgType',         failMsgType, ...
                         'saveType',            saveType);

fprintf('\n\nHere is the current PGPCA EM setting:\n\n');
disp(dispSettingStru);

%% Run EM iteration.
% Decide #iteration.
numIter = iterLimit;

% Collect variables for M-step I/O.
varStru = struct('sampY', sampY, 'sampZ', sampZ, 'sampZS', []);

% Run each iteration one by one.
fprintf('\n\nStart running PGPCA EMAlgo:\n\n');
for(j = 1:1:numIter)
    %% Compute one iteration.
    fprintf('PGPCA EMAlgo iteration: %3d\n', j);
    % Run E-step.
    [matPi, matGa, probQForZ, resEStepEach] = ...
        PGPCA.EStepAll(sampY, sampZ, paraStru, funStru, ...
                       'probQForZType',     probQForZType, ...
                       'probZInteWei',      probZInteWei, ...
                       'iniZS',             iniZS, ...
                       'saveType',          saveType);
    % Add fields to "resEStepEach".
    resEStepEach.matPi = matPi;
    resEStepEach.matGa = matGa;
    
    if( saveType == "all" )
        resEStepEach.probQForZ = probQForZ;
    else
        resEStepEach.probQForZ = [];
    end
    % (If available) Collect "sampZS" from E-step.
    if( probQForZType == "delta_func" )
        varStru.sampZS = resEStepEach.matZGp;
    end
    % Run M-step.
    [paraStru, funStru, dimC, resMStepEach] = ...
        PGPCA.MStepAll(matPi, matGa, probQForZ, dimC, varStru, funStru, ...
                       'probQForZType',     probQForZType, ...
                       'fitFunProbZSW',     fitFunProbZSW, ...
                       'leakSepaSW',        leakSepaSW, ... 
                       'failMsgType',       failMsgType);       
    % Add fields to "resMStepEach".
    resMStepEach.paraStru = paraStru;
    
    % Update all function handles.
    [~, funStru] = PGPCA.ParaFun2Model(dimY, dimL, paraStru, funStru, "probYlZType", "mvnpdf");

    % Compute "logLL".
    [~, resLogLLEach] = ...
        PGPCA.ModelPara2Prob("Probability", sampY, sampZ, paraStru, funStru, ...
                             "probValForm",     "ComplLL", ...
                             "saveType",        saveType);
    
    %% Sanity check: the logLL-related values must follow some orders.
    % Assert: M-step ELBO >= E-step ELBO.
    aveELBO_E = mean(resEStepEach.logLLProbQ);
    aveELBO_M = mean(resMStepEach.logLLProbQ);
    errTol = abs(mean([aveELBO_E, aveELBO_M])) * errRatio;
    assert( aveELBO_M - aveELBO_E >= -errTol, ...
            'In iteration %2d, M-step ELBO - E-step ELBO (%2.2e) should >= 0 (errTol: %2.2e).', ...
            j, aveELBO_M - aveELBO_E, -errTol );
    
    % (Only "posterior") assert: E-step ELBO (j) >= M-step ELBO (j-1).
    if(probQForZType == "posterior" && j >= 2)
        aveELBO_MPre = mean(resMStep(j-1).logLLProbQ);
        assert( aveELBO_E - aveELBO_MPre >= -errTol, ...
                'E-step ELBO [iter %d] - M-step ELBO [iter %d] (%2.2e) should >= 0 (errTol: %2.2e).', ...
                j, j-1, aveELBO_E - aveELBO_MPre, -errTol );
    end
    
    %% Collect results in this iteration.
    % Group "resColl" from E-step, M-step, and logLL.
    resEStep(j) = resEStepEach;
    resMStep(j) = resMStepEach;
    resLogLL(j) = resLogLLEach;
    
    % If j == 1, allocate the saving space for further iterations.
    if(j == 1)
        resEStep(numIter + 1) = resEStep(1);
        resMStep(numIter + 1) = resMStep(1);
        resLogLL(numIter + 1) = resLogLL(1);
        % Remove the last index to avoid misconfiguration.
        resEStep(numIter + 1) = [];
        resMStep(numIter + 1) = [];
        resLogLL(numIter + 1) = [];
    end
    
end

end

