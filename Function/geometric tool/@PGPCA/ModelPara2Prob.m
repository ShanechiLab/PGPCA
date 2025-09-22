function [probVal, resColl] = ModelPara2Prob(probMode, sampY, sampZ, paraStru, funStru, varargin)
%{
09-25-2023: IMPORTANT: THIS IS JUST A TEMPORARY FUNCTION FOR HAVING THE PRETEST RESULT FASTER. MANY OTHER
              GENERALIZATIONS AND OPTIMIZATIONS ARE NEEDED!

10-02-2023: (A) Add option "ComplExpo & SquNorm" so the output "probVal" is more flexible in E-step.
            (B) The noise COV at each "sampZ" is computed by "PGPCA.ModelCov" now for coding integrity.

10-12-2023: (A) Add option "JensenLowBd" in probMode "Probability" for the Jensen's lower bound.

10-26-2023: (A) Add option "saveType" to restrict the output fields in "resColl" for many EM iterations.
-----------------------------------------------------------------------------------------------------------

This function computes the probability-related values by model's parameters. The model is:

Model: y_t = phi(z_t) + K(z_t)*C*x_t + r_t [COV: R -> Ku(z_t) * diag(mainVar, sideVar) * Ku(z_t)']

To be useful in PGPCA EM algorithm, we have "two main modes" plus follow-up sub-modes:

1.  probMode = "ExpoPart":      The exponential components from function "MVNPDFExpoPart".

    (A) "ComplExpo":    a 2D array [numY, numZ]. The complete exponential part "exp(-1/2 * [squNorm])".
    (B) "SquNorm":      a 2D array [numY, numZ]. Only the squared norm part (no exp & -1/2)!

2.  probMode = "Probability":   The probability of each sample from function "MVNPDFBatch".

    (A) "MixProb":      a 2D array [numY, numZ]. The probability "p(y_i|z_j)".
    (B) "ComplProb":    a column [numY, 1]. The complete probability "p(y_i)".
    (C) "ComplLL":      a column [numY, 1]. The complete log-likelihood "ln(p(y_i))".
    (D) "JensenLowBd":  a column [numY, 1]. The Jensen's lower bound "ln(p(y_i|z_j)) * p(z_j)".

NOTE:
1.  Option "SquNorm" is useful in E-step state estimation by IVT.
2.  To compute "ComplProb & ComplLL", an additional function "funProbZ" must be provided.
3.  "JensenLowBd" is a lower bound of "ComplLL". We include this option for convenient analysis.
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
probMode:           a string "ExpoPart/Probability". It decides the type of output "probVal".
sampY:              a 2D array [numY, dimY]. Every row is an observation y_j.
sampZ:              a 2D array [numZ, dimZ]. Every row is an integral point in computing "valY(sampY)".
paraStru:           a struct. It contains all (learned) system parameters below.
    matC:               a 2D array [dimL, dimC] / empty []. The dimensionality reduction matrix.
    mainVar:            a scalar >= 0. The average noise VAR of the "modeled space".
    sideVar:            a scalar >= 0 / empty []. The average noise VAR of the "unmodeled space".
funStru:            a struct. It collects all function handles in following fields:
op  funProbZ:           a func. handle. The probability p(z).
    funMean:            a func. handle. The mean vector "phi(z)" for each "z" sample.
    funMatK:            a func. handle. The "orthonormal" coordinate matrix "K(z)" for each "z" sample.
------------------- Optional ('tag' + 'value' with (number))
"probValForm":      (1) (probValForm) a string depending on input "probMode". (def: see below)
                        It decides the output form of "probVal". The details are above. Basically,
    w/ "ExpoPart":          Options are "ComplExpo / SquNorm". (def: "ComplExpo")
    w/ "Probability":       Options are "MixProb / ComplProb / ComplLL / JensenLowBd". (def: "MixProb")
"saveType":         (1) (saveType) a string "all/core" (def: "all").
                        It decides the recorded variables in "resColl" for saving memory.
    "all":                  all variables are saved.
    "core":                 only the main variables with reasonable size are saved. Those fields tagged 
                              "op" in the output list are ignored with field = empty [].

Outputs:
------------------- Essential
probVal:            a 2D array [numY, numZ] / a column [numY, 1]. It depends on "probMode & probValForm".
resColl:            a struct. It records intermediate values during computation in fields below:
    probMode:           a string. The main mode of the probability computation.
    probValForm:        a string. The sub-mode of the probability computation.
    probVal:            a 2D array [numY, numZ] / a column [numY, 1]. It just copies the 1st output.
    mvnMean:            a 2D array [dimY, numZ]. Each column [:,j] is "phi(z)" at matching "sampZ".
    mvnCov:             a 3D array [dimY, dimY, numZ]. Each slice [:,:,j] is COV of "y_t" at "sampZ".
op  covResColl:         a struct. The info. in computing "mvnCov" by "PGPCA.ModelCov".
op  mvnResColl:         a struct. The info. in computing MVNPDF-related values.
%}

%% Assign inputs.
% Default - "probValForm" depends on "probMode".
switch probMode
    case "ExpoPart"
        probValFormDef   = "ComplExpo";
        probValFormCandi = ["ComplExpo", "SquNorm"];
    case "Probability"
        probValFormDef   = "MixProb";
        probValFormCandi = ["MixProb", "ComplProb", "ComplLL", "JensenLowBd"];
    otherwise
        error('The given "probMode" (%s) is invalid.', probMode);
end

% Parser for optional inputs.
optionP = inputParser;
optionP.addParameter('probValForm',     probValFormDef,     @(x) ismember(x, probValFormCandi));
optionP.addParameter('saveType',        "all",              @(x) any(x == ["all", "core"]));

% Parse the inputs.
optionP.parse(varargin{:});
probValForm = optionP.Results.probValForm;
saveType    = optionP.Results.saveType;

%% Prepare "mean & COV".
% Compute mean.
mvnMean = funStru.funMean(sampZ);
% Formalize its form.
mvnMean = mvnMean';

% Compute COV.
[mvnCov, covResColl] = PGPCA.ModelCov("covY", sampZ, paraStru, funStru);

%% Compute "probVal".
% Check "probMode".
switch probMode
    case "ExpoPart"
        %% CASE: Only the "exponential part" of normal distribution.
        % Decide the output is the "complete exponential part" or not.
        switch probValForm
            case "ComplExpo"
                pureNormSW = false;
            case "SquNorm"
                pureNormSW = true;
            otherwise
                error('The given "probValForm" (%s) is invalid.', probValForm);
        end
        
        % Compute: the mvnpdf "exponential part".
        [probVal, mvnResColl] = MVNPDFExpoPart(sampY, mvnMean, mvnCov, 'pureNormSW', pureNormSW);

    case "Probability"
        %% CASE: The probability of normal distribution.
        [probVal, mvnResColl] = MVNPDFBatch(sampY, mvnMean, mvnCov);
        
        % Prepare "probZ" if needed.
        if( probValForm ~= "MixProb" )
            probZ = funStru.funProbZ(sampZ);
        end
        
        % Further process "probVal" following "probValForm".
        switch probValForm
            case "MixProb"
                % CASE: the conditional probability p(y_i|z_j).
                % Do nothing.
                
            case "ComplProb"
                % CASE: the complete probability.
                probVal = probVal * probZ;
                
            case "ComplLL"
                % CASE: the complete log-likelihood.
                probVal = log(probVal * probZ);
                
            case "JensenLowBd"
                % CASE: the Jensen's lower bound.
                probVal = log(probVal) * probZ;
                
            otherwise
                error('The given "probValForm" (%s) is invalid.', probValForm);
        end
        
    otherwise
        error('The given "probMode" (%s) is invalid.', probMode);
end

% Finally, make sure "probVal" is properly numeric.
assert( all(IsNumericProper(probVal), 'all'), '"probVal" must be properly numeric!' );

%% Prepare outputs.
resColl = struct('probMode',        probMode, ...
                 'probValForm',     probValForm, ...
                 'probVal',         probVal, ...
                 'mvnMean',         mvnMean, ...
                 'mvnCov',          mvnCov, ...
                 'covResColl',      [], ...
                 'mvnResColl',      []);
             
% (If asked) save all variables.
if( saveType == "all" )
    resColl.covResColl = covResColl;
    resColl.mvnResColl = mvnResColl;
end
    
end

