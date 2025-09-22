function [AHat, transInfo] = PDTrans(A, propType, mdType, varargin)
%{
02/07/2022: 
1.  Renew the function so it can handle vector/matrix & PD/PSD various kinds of cases.

07/25/2022:
1.  Add option "negaTol" to avoid/detect very negative matrix (like eig = -1e3).

08/19/2022:
1.  Renew option "negaTol" so it can handle "scale/absolute" cases. Remove symmetric calibration.

09/07/2022:
1.  Renew output "transInfo" so #layer in struct is only 1. This is critical for reducing the memory.
-----------------------------------------------------------------------------------------------------------

This function makes a symmetric matrix positive-definite (PD) or positive semi-definite (PSD). PD/PSD of 
  the matrix "A" is determined by "eig" & "det".

(A) The formats of matrix A:
    1.  vector:         a column/row. The diagonal of the matrix.
    2.  2D array:       a squared matrix.

(B) The properties of matrix A:
    1.  PD:             positive definite.
    2.  PSD:            positive semi-definite.

(C) The calibration methods:
    1.  "constBias":    add a constant bias (NOTE: this output may still not be PD!).
    2.  "eig":          make the smallest eigenvalue be positive. (+ default constant bias).
    3.  "nearSPD":      find the nearest symmetric PD matrix.
    4.  "none":         no calibration. In this case, the function just test the property.

(D) The awaring/failing messages (when the final type doesn't match):
    1.  "error":        the program will be stopped w/ an error message.
    2.  "warning":      the program keeps running w/ a warning message.

NOTE:
1.  To avoid mis-operation, we show an error when "A" is not symmetric (even just tiny mismatch)! The user
      should guarantee this. This way is less error-prone.
-----------------------------------------------------------------------------------------------------------

The tolerance of the initial negative eigenvalues in "A":

By default, input "A" cannot be too negative. This function is used to fix numerial error of negative 
  "A" when it should be positive mathematically, so very negative eigs usually means some problems. The
  question is: what is too negative. This relates to "rcond" of A. Therefore, there are two methods:

1.  negaTolMd   = "scale".
    negaTolVal  = a scaling scalar in [0,1).

->  In this case, the tolerance "negaTolBd = MaxEig(A) * negaTolVal * (-1)".
->  It depends on the range of elements in "A". MaxEig(A) must be positive or an error shows.

2.  negaTolMd   = "absolute":
    negaTolVal  = a non-positive scalar.

->  In this case, the tolerance "negaTolBd = negaTolVal"
->  It's a fixed value no matter the range of "A". This is useful when some conditions are very strict.
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
A:                  a 2D square array [dimA] / a vector [1, dimA] or [dimA, 1]. The targeted matrix.
propType:           a string "PD / PSD". The required property of the matrix "A".
mdType:             a string "constBias / eig / nearSPD / none". The method for the calibration. 
------------------- Optional ('tag' + 'value' with (number))
"failType":         (1) (failType) a string "error / warning" (def: "error")
                        When final "AHat" still doesn't satisfy the property, "error" stops the function, 
                          and "warning" just show a warning message (the program can keep running).
"formAHat":         (1) (formAHat) a string "row/column/2D square" / empty []. (def: [])
                        The matrix will be transformed to this format at the output. [] means the same as 
                          the input (i.e., no transformation).
"negaTol":          (3) (negaTolSW)  a boolean. (def: "true")
                        (negaTolMd)  a string "scale / absolute" / empty []. (def: [])
                        (negaTolVal) a scalar / empty []. (def: [])
                        When "negaTolSW" is true, the most negative eigenvalue of "A" must be larger than
                          "negaTolBd" whose computation is listed above. Empty [] means default settings.

Outputs:
------------------- Essential
AHat:               a 2D square matrix [dimA, dimA]. The transformed A.
transInfo:          a struct. It records information about this transformation in the following fields:
    ------------------- (BASIC SETTING)
    propType:           a string. The required property.
    mdType:             a string. The calibrating method.
    ------------------- (BASIC PD/PSD EXAMINATION)
    fixFlag:            a boolean. "true" means the "Ahat" is modified.
    badFlag:            a boolean. "true" means the "Ahat" is still not PD. 
    oldMinEig:          a scalar. The minimal eigenvalue BEFORE the tranformation.
    minEig:             a scalar. The minimal eigenvalue AFTER the tranformation.
    detA:               a scalar. The determinant "det(AHat)".
    ------------------- (SPECIAL: THE TOLERANCE OF NEGATIVE EIGENVALUES)
    negaTolSW:          a boolean. "true" means the initial most negative Eig cannot exceed "negaTolBd".
    negaTolMd:          a string. The method to generate "negaTolBd".
    negaTolVal:         a scalar. The parameter to compute "negaTolBd".
    negaTolBd:          a -scalar. The negative lower bound of the initial most negative Eig.
    negaBadFlag:        a boolean. "true" means the initial most negative Eig < "negaTolBd". This is bad.
%}      

%% Setting.
% Constant - Eigenvalue-related settings.
defShiftBias        =  1e-10;
defNegaTolMd        =  "scale";
defNegaTolScale     =  1e-10;
defNegaTolBd        = -1e-10;
% Constant - variable status/structure.
fixFlag = false;
transInfo = struct('propType', propType, 'mdType', mdType, 'fixFlag', [], ...
                   'badFlag',     [], 'oldMinEig', [], 'minEig',     [], 'detA',      [], ...
                   'negaTolSW',   [], 'negaTolMd', [], 'negaTolVal', [], 'negaTolBd', [], ...
                   'negaBadFlag', []);

% Default - user-defined.
failType = "error";
formAHat = [];
negaTolSW = true;
negaTolMd = [];
negaTolVal = [];

%% Assign options (if they exist).
numVarargin = length(varargin);
if( numVarargin > 0 )
    tagInd = 1;
    while( tagInd <= numVarargin )
        % Select the option.
        tagName = varargin{tagInd};
        switch tagName
            case 'failType'
                tagShift = 1;
                failType = varargin{ tagInd + 1 };
            case 'formAHat'
                tagShift = 1;
                formAHat = varargin{ tagInd + 1 };
            case 'negaTol'
                tagShift = 3;
                negaTolSW = varargin{ tagInd + 1 };
                negaTolMd = varargin{ tagInd + 2 };
                negaTolVal = varargin{ tagInd + 3 };
            otherwise
                error('The given option "%s" is invalid.', tagName);
        end
        % Shift the index.
        tagInd = tagInd + tagShift + 1;
    end
end

%% Decide the form of "A".
[dimR, dimC] = size(A);
if( dimR == dimC )
    % CASE: 2D square.
    formA = "2D square";
    dimA = dimR;
    % Assert: A is symmetric.
    assert( issymmetric(A), 'Input "A" is not symmetric!' );
    
elseif( dimR == 1 )
    % CASE: row.
    formA = "row";
    dimA = dimC;
    
elseif( dimC == 1 )
    % CASE: column.
    formA = "column";
    dimA = dimR;
    
else
    % Not accepted format.
    error('The dimension of A [%d, %d] must be a vector/2D square.', dimR, dimC)
end

%% Sanity check the negative bound.
% Assert: "negaTolSW" must be a boolean.
assert( islogical(negaTolSW), 'The "negaTolSW" class (%s) must be boolean.', class(negaTolSW) );
% Initialize "negaTolMd".
if( isempty(negaTolMd) ); negaTolMd = defNegaTolMd; end

% Compute the negative tolerance bound.
switch negaTolMd
    case "scale"
        % CASE: tolerance "negaTolBd = MaxEig(A) * negaTolVal * (-1)".
        if( isempty(negaTolVal) ); negaTolVal = defNegaTolScale; end
        % Assert: it's in [0, 1).
        assert( 0 <= negaTolVal && negaTolVal < 1, '"negaTolVal" (%2.2e) must be in [0,1).', negaTolVal );
        % Get the maximal eigenvalue & make sure it's non-negative.
        maxEigA = eig(A);
        assert( maxEigA(end) >= 0, 'The max eig of "A" (%2.2e) must be >= 0.', maxEigA(end) );
        % Avoid maxEigA is not too small.
        maxEigA = max( maxEigA(end), defShiftBias );
        % Set the tolerance boundary.
        negaTolBd = maxEigA * negaTolVal * (-1);
        
    case "absolute"
        % CASE: tolerance "negaTolBd = negaTolVal".
        if( isempty(negaTolVal) ); negaTolVal = defNegaTolBd; end
        % Set the tolerance boundary.
        negaTolBd = negaTolVal;
    
    otherwise
        error('The given negative tolerance method "%s" is invalid.', negaTolMd);
end

% Assert: the eigenvalue lower bound must be negative.
assert( negaTolBd <= 0, 'The eigenvalue lower bound (%2.2e) must be non-positive.', negaTolBd );
% Collect the result.
transInfo.negaTolSW     = negaTolSW;
transInfo.negaTolMd     = negaTolMd;
transInfo.negaTolVal    = negaTolVal;
transInfo.negaTolBd     = negaTolBd;

%% Check if "A" is good already.
[badFlag, checkInfo] = CheckMatrixProp(A, formA, propType, negaTolSW, negaTolBd);
oldMinEig = checkInfo.minEig;
if( badFlag )
    fixFlag = true;
else
    % Prepare the output & following the format.
    AHat = TransMatrixForm(A, formAHat);
    % Basic fields.
    transInfo.fixFlag       = fixFlag;
    transInfo.badFlag       = badFlag;
    transInfo.oldMinEig     = oldMinEig;
    % PD/PSD fields.
    transInfo.minEig        = checkInfo.minEig;
    transInfo.detA          = checkInfo.detA;
    transInfo.negaBadFlag   = checkInfo.negaBadFlag;    
    return;
end

%% Prepare the variable.
% The shifting array.
switch formA
    case "row"
        shiftArray = ones(1, dimA);
    case "column"
        shiftArray = ones(dimA, 1);
    case "2D square"
        shiftArray = eye(dimA);
    otherwise
        error('The matrix "A" format (%s) is invalid.', formA);
end

% The shifting bias.
switch propType
    case "PD"
        shiftBias = defShiftBias;
    case "PSD"
        % We use "eps" to avoid numerical error!
        shiftBias = eps;
    otherwise
        error('The given property "%s" is invalid.', propType);
end

%% Transform "A" to "AHat".
switch mdType
    case 'constBias'
        % CASE: add the constant bias. Note that the new A may still not be PD!
        AHat = A + shiftBias * shiftArray;
        
    case 'eig'
        % CASE: find the most negative eigenvalues, and then make it be "addBias".
        AHat = A + (shiftBias - checkInfo.minEig) * shiftArray;
        
    case 'nearSPD'
        % CASE: find the nearest SPD matrix w.r.t A.
        % WARNING: if "A" is a vector, this may cause unexpected behavior!
        if( any(formA == ["row", "column"]) )
            warning('Format "%s" using calibration "%s" may induces problems!', formA, mdType);
        end
        if( propType == "PSD" )
            warning('Property "%s" using calibration "%s" is actually "PD"!', propType, mdType);
        end
        % Transform "A" to 2D square, calibrate it, and then change it back.
        transA = TransMatrixForm(A, '2D square');
        [AHat, ~] = NearestSPD(transA);
        AHat = TransMatrixForm(AHat, formA);
        
    case 'none'
        % CASE: no calibration.
        AHat = A;
        % NOTE: this is a special case that "A" is not fixed in any way, so "fixFlag = false".
        fixFlag = false;
        
    otherwise
        error('The given option is "%s". Please fix it.', mdType);
end

%% Prepare the output.
[badFlag, checkInfo] = CheckMatrixProp(AHat, formA, propType, negaTolSW, negaTolBd);
if( badFlag )
    switch failType
        case "error"
            error('After method "%s", matrix "A" is still not PD! (minEig: %2.2e, det: %2.2e)', ...
                  mdType, checkInfo.minEig, checkInfo.detA)
              
        case "warning"
            warning('After method "%s", matrix "A" is still not PD! (minEig: %2.2e, det: %2.2e)', ...
                    mdType, checkInfo.minEig, checkInfo.detA);
            
        otherwise
            error('The criterion option "%s" is invalid.', failType);
    end
end

% Prepare the output.
AHat = TransMatrixForm(AHat, formAHat);
% Basic fields.
transInfo.fixFlag       = fixFlag;
transInfo.badFlag       = badFlag;
transInfo.oldMinEig     = oldMinEig;
% PD/PSD fields.
transInfo.minEig        = checkInfo.minEig;
transInfo.detA          = checkInfo.detA;
transInfo.negaBadFlag   = checkInfo.negaBadFlag;    

end


function [badFlag, checkInfo] = CheckMatrixProp(A, formA, propType, negaTolSW, negaTolBd)
%{
This function checks if the symmetric matrix "A" satisfies the given property. The cases are:
1.  "PD":       min(eig) > 0  & det(A) > 0.
2.  "PSD":      min(eig) >= 0 & det(A) >= 0.

Also, if asked, the minimal eigenvalue must be larger than the lower bound, or the error is shown!

Inputs:
------------------- Essential
A:                  a 2D square [dimA] / a vector. A symmetric matrix.
formA:              a string "row/column/2D square". The format of "A".
propType:           a string "PD/PSD". The required property of "A".
negaTolSW:          a boolean. "true" means the minimal eigenvalue of "A" must satisfy the bound.
negaTolBd:          a -scalar. The lower bound of the negative eigenvalues.

Outputs:
------------------- Essential
badFlag:            a boolean. "true" means "A" doesn't satisfy the property "PD/PSD".
checkInfo:          a struct. It records the checking information including:
    minEig:             a scalar. The minimal eigenvalue of "A".
    detA:               a scalar. The determinant of "A".
    negaTolSW:          a boolean. "true" means the most negative eigenvalue of "A" is bounded.
    negaTolBd:          a negative scalar. The lower bound of the negative eigenvalues.
    negaBadFlag:        a boolean / empty []. "true" means the "minEig" fails the eigenvalue lower bound.
%}

%% Function block.
checkInfo = struct('minEig', [], 'detA', [], 'negaTolSW', [], 'negaTolBd', [], 'negaBadFlag', []);
% Compute measurements.
switch formA
    case "2D square"
        checkInfo.minEig = min(eig(A));
        checkInfo.detA = det(A);
    case {"row", "column"}
        checkInfo.minEig = min(A);
        checkInfo.detA = prod(A);
    otherwise
        error('The format of "A" (%s) is invalid.', formA);
end

% If asked, check the eigenvalue bound.
if( negaTolSW )
    if( checkInfo.minEig >= negaTolBd )
        % CASE: the minEig >= the lower bound. This is good.
        checkInfo.negaBadFlag = false;
    else
        % CASE: the minEig < the lower bound. This is bad.
        checkInfo.negaBadFlag = true;
        error( 'The "minEig" (%2.2e) must be >= the lower bound (%2.2e).', checkInfo.minEig, negaTolBd );
    end
end

% Check if "A" follows the property.
switch propType
    case "PD"
        passFlag = (checkInfo.minEig > 0 && checkInfo.detA > 0);
    case "PSD"
        passFlag = (checkInfo.minEig >= 0 && checkInfo.detA >= 0);
    otherwise
        error('The property of "A" (%s) is invalid.', propType);
end

% Create the output.
badFlag = ~passFlag;

% Collect the result.
checkInfo.negaTolSW = negaTolSW;
checkInfo.negaTolBd = negaTolBd;

end


function transA = TransMatrixForm(A, formTrans)
%{
This function transform the matrix following the given format. The output format can be:
1.  "row":              a row vector.
2.  "column":           a column vector.
3.  "2D square":        a 2D squared array.
4.  []:                 nothing. Following the original format.

The format of "A":
1.  vector:             a column/row. The diagonal of the matrix.
2.  2D array:           a squared matrix.

Inputs:
--------------- Essential
A:              a 2D square [dimA] / a vector. A symmetric matrix.
formTrans:      a string "row/column/2D square" / empty []. The transformed format of "A".

Outputs:
--------------- Essential
transA:         a 2D square [dimA] / a vector. The transformed matrix "A".
%}

%% Function block.
if( isempty(formTrans) )
    % Do nothing.
    transA = A;
else
    % Decide the transformation.
    switch formTrans
        case "row"
            % CASE: a row.
            if( isrow(A) );             transA = A;
            elseif( iscolumn(A) );      transA = A';
            else;                       transA = diag(A)';
            end
            
        case "column"
            % CASE: a column.
            if( isrow(A) );             transA = A';
            elseif( iscolumn(A) );      transA = A;
            else;                       transA = diag(A);
            end
            
        case "2D square"
            % CASE: a 2D squared array.
            if( isvector(A) );          transA = diag(A);
            else;                       transA = A;
            end
            
        otherwise
            error('The transforming format "%s" is invalid.', formTrans);
    end
end

end

