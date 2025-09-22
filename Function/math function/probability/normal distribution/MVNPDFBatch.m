function [probMat, resColl] = MVNPDFBatch(sampY, mvnMean, mvnCov, varargin)
%{
This function computes "multiple samples" in "multiple normal distributions" in BATCH for speed-up. It's
  different from MATLAB "mvnpdf" because sampling probabilities from single model have shared components.
  Not computing redundant/unnecessay values is critical when the #sample is large!

Goal:   Compute the exponential function of "each sample" under "each model" efficiently.

Func:   1/[ (2*pi)^{n/2} * |Cov|^{1/2} ] * exp{ -1/2 * (Cov-based squared norm) }

Step:   Compute "sample-model" pair in batch to avoid any unnecessary computation.

1.  Compute "exponential part" by function "MVNPDFExpoPart".
2.  Compute "determinant" for each model (not each sample!).
3.  Compute the final probability.
-----------------------------------------------------------------------------------------------------------

NOTE:
1.  It's based on MATLAB function "mvnpdf", but more general and efficient.
2.  "mvnMean & mvnCov" are with model index at the end considering the slicing efficiency.
3.  Since det(A^{-1}) = [det(A)]^{-1}, given that "mvnCov is inverse" won't slow the computation.
4.  Some fields in "resColl" are inherent from "MVNPDFExpoPart", so they're not options here.
5.  When "pairModeSW = true", It reduces to MATLAB "mvnpdf" case such that each sample map to one system 
      only, so output "expMat" is a vector.
6.  In some cases, given the pre-computed determinant can speed-up the process a lot! Note that we don't
      double-check its consistency with "mvnCov". Otherwise, given "mvnCovDet" is meaningless.
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
sampY:              a 2D array [numSamp, dimY]. Each row is a sample.
mvnMean:            a 2D array [dimY, 1/numModel]. Each row is the mean of the corresponding sample.
mvnCov:             a 3D array [dimY, dimY, 1/numModel]. Each slice [:,:,j] is a covariance matrix.
------------------- Optional ('tag' + 'value' with (number))
"invCovSW":         (1) (invCovSW) a boolean (def: false)
                        "true" means "mvnCOV" is the inverse matrix already. This can save time!
"pairModeSW":       (1) (pairModeSW) a boolean (def: false)
                        "true" means "sampY" and "models" are 1-1 mapping. Details are above.
"mvnCovDet":        (1) (mvnCovDet) a vector [1, numModel] (def: [])
                        The user-given determinant of COV (not inverse COV!).

Outputs:
------------------- Essential
probMat:            a 2D array [numSamp, 1/numModel]. Each element is the experiment result.
resColl:            a struct. It collects all results in following fields:
    invCovSW:           a boolean. "true" means "mvnCov" is the inversed matrix already.
    pureNormSW:         a boolean. "true" means output "expMat = S (above)", no "exp() & -1/2".
    pairModeSW:         a boolean. "true" means "sampY & models" are 1-1 paired.
    dimY:               a +integer. The dimension of "sampY".
    numSamp:            a +integer. The number of samples.
    numModel:           a +integer. The number of probabilistic models.
    expMat:             a 2D array [numSamp, 1/numModel]. The exponential part from "MVNPDFExpoPart".
    mvnCovDet:          a vector [1, numModel]. The COV determinant of each model.
%}

%% Assign inputs.
% Constant.
dimY = size(sampY, 2);
pureNormSW = false;

% Default.
invCovSW = false;
pairModeSW = false;
mvnCovDet = [];

% Assign options (if exist.)
numVarargin = length(varargin);
if(numVarargin > 0)
    tagInd = 1;
    while(tagInd <= numVarargin)
        % Select the option.
        tagName = varargin{tagInd};
        switch tagName
            case "invCovSW"
                tagShift = 1;
                invCovSW = varargin{tagInd + 1};
            case "pairModeSW"
                tagShift = 1;
                pairModeSW = varargin{tagInd + 1};
            case "mvnCovDet"
                tagShift = 1;
                mvnCovDet = varargin{tagInd + 1};
            otherwise
                error('The given option "%s" is invalid.', tagName);
        end
        % Shift the index.
        tagInd = tagInd + tagShift + 1;
    end
end

%% Compute: the exponential part.
% NOTE: the sanity check is also done by "MVNPDFExpoPart". This improves coding integrity.
[expMat, resColl] = MVNPDFExpoPart(sampY, mvnMean, mvnCov, ...
                                   'invCovSW', invCovSW, ...
                                   'pairModeSW', pairModeSW, ...
                                   'pureNormSW', pureNormSW);

% Extract some values.
numModel = resColl.numModel;

%% Compute: the determinant.
% Assert: "invCovSW" is boolean (just to play safe).
invCovSWClass = class(invCovSW);
assert( invCovSWClass == "logical", ...
        'The class of "invCovSW" (%s) must be logical.', invCovSWClass );

% (If not given) compute determinant.
if( isempty(mvnCovDet) )
    mvnCovDet = zeros(1, numModel);
    % Compute determinant one by one.
    for(j = 1:1:numModel)
        mvnCovDet(j) = det(mvnCov(:,:,j));
    end
    % Inverse the result if the given Cov are inverse already.
    if( invCovSW )
        mvnCovDet = 1 ./ mvnCovDet;
    end
end

% Assert: the spec of "mvnCovDet" is correct. This is a sanity check.
mvnCovDetNum = length(mvnCovDet);
assert( isvector(mvnCovDet), ...
        '"mvnCovDet" should be a vector.' );
assert( mvnCovDetNum == numModel, ...
        '"mvnCovDet" length (%d) must = "numModel" (%d).', mvnCovDetNum, numModel );
assert( all(mvnCovDet > 0), ...
        '"mvnCovDet" components must be > 0.' );

%% Compute: the final probability.
% Check the output spec.
if( pairModeSW )
    % CASE: samples are models are paired, so expMat = [numSamp = numModel, 1].
    probMat = expMat ./ sqrt(mvnCovDet)';
else
    % CASE: expMat = [numSamp, numModel].
    probMat = expMat ./ sqrt(mvnCovDet); 
end

% Finalize "probMat" by the constant part.
probMat = probMat / (2*pi)^(dimY/2);

%% Prepare outputs.
resColl.expMat = expMat;
resColl.mvnCovDet = mvnCovDet;

end

