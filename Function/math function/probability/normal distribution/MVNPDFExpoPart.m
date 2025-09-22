function [expMat, resColl] = MVNPDFExpoPart(sampY, mvnMean, mvnCov, varargin)
%{
This function is a component of my multi-variate normal distribution function "MVNPDFBatch". I separate it 
  because it can be independently applied to, for example, parameter optimization sometimes.

Goal:   Compute the exponential function of "each sample" under "each model" efficiently.

Step:   Compute "sample-model" pair in batch to avoid any unnecessary computation.

1.  Transpose "sampY = sampY'" to align its dimension format with "mvnMean & mvnCov".
2.  expMat(i,j) = exp{ -1/2 * S }
    where S = (sampY(:,i) - mvnMean(:,j))' * inv(mvnCov(:,:,j)) * (sampY(:,i) - mvnMean(:,j))
-----------------------------------------------------------------------------------------------------------

Option: "pureNormSW" decides if the output is "S" only, no "exp() & -1/2". 

Option: "pairModeSW" decides the mapping between "sampY" and "mvnMean & mvnCov".

(A) pairModeSW = false:

    All samples map to all models, so the output is a 2D array [numSamp, numModel].

(B) pairModeSW = true:

    Each sample maps to its corresponding model, so the output is a vector [numSamp, 1]. Therefore,
      Condition "numSamp == numModel" must be hold.
-----------------------------------------------------------------------------------------------------------

NOTE:
1.  It's based on MATLAB function "mvnpdf", but more general and efficient.
2.  For complete normal pdf, please refer to function "MVNPDFBatch".
3.  "mvnMean & mvnCov" are with model index at the end considering the slicing efficiency.


4.  When "pairModeSW = true", It reduces to MATLAB "mvnpdf" case such that each sample map to one system 
      only, so output "expMat" is a vector.
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
sampY:              a 2D array [numSamp, dimY]. Each row is a sample.
mvnMean:            a 2D array [dimY, 1/numModel]. Each row is the mean of the corresponding sample.
mvnCov:             a 3D array [dimY, dimY, 1/numModel]. Each slice [:,:,j] is a covariance matrix.
------------------- Optional ('tag' + 'value' with (number))
"invCovSW":         (1) (invCovSW) a boolean (def: false)
                        "true" means "mvnCOV" is the inverse matrix already. This can save time!
"pureNormSW":       (1) (pureNormSW) a boolean (def: false)
                        "true" means output "expMat = S (above)", no "exp() & -1/2".
"pairModeSW":       (1) (pairModeSW) a boolean (def: false)
                        "true" means "sampY" and "models" are 1-1 mapping. Details are above.

Outputs:
------------------- Essential
expMat:             a 2D array [numSamp, 1/numModel]. Each element is the experiment result.
resColl:            a struct. It collects all results in following fields:
    invCovSW:           a boolean. "true" means "mvnCov" is the inversed matrix already.
    pureNormSW:         a boolean. "true" means output "expMat = S (above)", no "exp() & -1/2".
    pairModeSW:         a boolean. "true" means "sampY & models" are 1-1 paired.
    dimY:               a +integer. The dimension of "sampY".
    numSamp:            a +integer. The number of samples.
    numModel:           a +integer. The number of probabilistic models.
%}

%% Assign inputs.
% Constant.
[numSamp, dimY] = size(sampY);

% Default.
invCovSW = false;
pureNormSW = false;
pairModeSW = false;

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
            case "pureNormSW"
                tagShift = 1;
                pureNormSW = varargin{tagInd + 1};
            case "pairModeSW"
                tagShift = 1;
                pairModeSW = varargin{tagInd + 1};
            otherwise
                error('The given option "%s" is invalid.', tagName);
        end
        % Shift the index.
        tagInd = tagInd + tagShift + 1;
    end
end

%% Prepare parameters & sanity check.
% Extract sizes (NOTE: THESE VALUES WON'T BE USED AFTER THIS SECTION!)
[dimMean, numMean] = size(mvnMean);
[dimCov1, dimCov2, numCov] = size(mvnCov);

% Compute "numModel".
numModel = CopyDataCheck([numMean, numCov]);
% Create the index set for "mean & cov".
[~, mvnMeanIndSet] = CheckMultiDimIndex(numModel, numMean);
[~, mvnCovIndSet] = CheckMultiDimIndex(numModel, numCov);

%% Sanity check.
% Assert: "dimY == dimMean".
assert( dimY == dimMean, ...
        'The dimension of "sampY" (%d) and "mvnMean" (%d) must be the same.', dimY, dimMean );

% Assert: mvnCov(:,:,j) is squared and its dimension is correct.
assert( dimCov1 == dimCov2, ...
        'The slice of "mvnCov" [%d, %d, j] must be a squared.', dimCov1, dimCov2 );
assert( dimCov1 == dimY, ...
        'The dimension of "sampY" (%d) and "mvnCov" (%d) must be the same.', dimY, dimCov1' );

% Assert: "invCovSW" is boolean.
invCovSWClass = class(invCovSW);
assert( invCovSWClass == "logical", ...
        'The class of "invCovSW" (%s) must be logical.', invCovSWClass );
    
% Assert: "pureNormSW" is boolean.
pureNormSWClass = class(pureNormSW);
assert( pureNormSWClass == "logical", ...
        'The class of "pureNormSW" (%s) must be logical.', pureNormSWClass );

% Assert: "pairModeSW" is boolean.    
pairModeSWClass = class(pairModeSW);
assert( pairModeSWClass == "logical", ...
        'The class of "pairModeSW" (%s) must be logical.', pairModeSWClass );

% Assert: if "pairModeSW = true", then "numSamp = numModel".
assert( (pairModeSW == false) || (numSamp == numModel), ...
        'When "pairModeSW = true", "numSamp" (%d) must = "numModel" (%d).', numSamp, numModel );
 
%% Compute the exponential.
% Customize the output form.
if(pairModeSW); expMat = zeros(numSamp, 1); else; expMat = zeros(numSamp, numModel); end

% Since "inv(mvnCov)" is the bottleneck, we compute all samples for each model in batch.
for(j = 1:1:numModel)
    % The selected difference and Cov.
    if(pairModeSW)
        mvnDiff = sampY(j,:) - mvnMean( :, mvnMeanIndSet(j) )';
    else
        mvnDiff = sampY - mvnMean( :, mvnMeanIndSet(j) )';
    end
    seleCov = mvnCov( :, :, mvnCovIndSet(j) );
    
    % Decide the inversion or not.
    if(invCovSW)
        mvnDiffInv = mvnDiff * seleCov;
    else
        mvnDiffInv = mvnDiff / seleCov;
    end
    
    % Only compute the diagonal term in Y*D^{-1}*Y for saving time.
    if(pairModeSW)
        expMat(j) = sum(mvnDiffInv .* mvnDiff, 2);
    else
        expMat(:,j) = sum(mvnDiffInv .* mvnDiff, 2);
    end
end

% Exponential it.
if( ~pureNormSW )
    expMat = exp( (-1/2)*expMat );
end

%% Prepare output.
resColl = struct('invCovSW',    invCovSW, ...
                 'pureNormSW',  pureNormSW, ...
                 'pairModeSW',  pairModeSW, ...
                 'dimY',        dimY, ...
                 'numSamp',     numSamp, ...
                 'numModel',    numModel);

end

