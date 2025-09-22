function [givenFlag, indRef, multiDInfo] = CheckMultiDimIndex(dimVar, dimRef, varargin)
%{
This function checks if the variable "dimVar" can match the reference "dimRef". The cases are:

1.  dimRef = a +integer:    It's the numel(dimVar).
2.  dimRef = a vector:      It's the dimVar.

The followings are examples: assume dimVar = [3,4,5], then

Example 1:  dimRef = 1          ->  A constant.         -> indRef = ones(1,60)
Example 2:  dimRef = 2          ->  2 ~= 60.            -> Error.
Example 3:  dimRef = 60         ->  The numel matches.  -> indRef = [1:1:60]
Example 4:  dimRef = [1,1,1]    ->  [3,4,5] ~= [1,1,1]. -> Error.
Example 5:  dimRef = [3,4,5]    ->  The size matches.   -> indRef = [1:1:60]

There are options can treat special cases:

1.  If "prodOnly" is true:      
->  We only care prod(dimRef) even it's a vector.

2.  If "prodIfVector" if true:
->  We only care prod(dimRef) if dimRef represents a vector.
->  For example: dimRef = [4,1], [1,3], [1,1], but not [2,3,4] or [1,1,1]).
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
dimVar:             a row [1, numDim]. The multi-D dimension of an external variable outside the function. 
dimRef:             a row [1, numDim] / a scalar. The reference that "dimVar" should match.
------------------- Optional ('tag' + 'value' with (number))
"prodOnly":         (1) (prodOnly) a boolean (def: false)
                        "true" means only prod(dimRef) is important (see above).
"prodIfVector":     (1) (prodIfVector) a boolean (def: false)
                        "true" takes effect when "dimRef" represents a vector (see above).

Outputs:
------------------- Essential
givenFlag:          a boolean. "true" means "dimVar" does NOT represent an empty variable (e.g., [0,0]).
indRef:             a row [1, numVar]. It depends on if "dimRef" represents a scalar or a vector.
multiDInfo:         a struct. It records info. about this function in following fields:
    prodOnly:           a boolean. It's the same as the option "prodOnly".
    prodIfVector:       a boolean. It's the same as the option "prodIfVector".
    dimRefType:         a string "scalar/vector/multiD". The type of reference represented by "dimRef".
%}

%% Assign inputs.
% Constant.
numVar = prod(dimVar);
numRef = prod(dimRef);

% Default.
prodOnly = false;
prodIfVector = false;

% Assign options (if asked).
numVarargin = length(varargin);
if( numVarargin > 0 )
    tagInd = 1;
    while( tagInd <= numVarargin )
        % Select the option.
        tagName = varargin{tagInd};
        switch tagName
            case 'prodOnly'
                tagShift = 1;
                prodOnly = varargin{tagInd + 1};
                
            case 'prodIfVector'
                tagShift = 1;
                prodIfVector = varargin{tagInd + 1};
                
            otherwise
                error('The given option "%s" is invalid.', tagName);
        end
        % Shift the index.
        tagInd = tagInd + tagShift + 1;
    end
end

%% Initialize the output & Return if "dimVar" represents an empty variable.
multiDInfo = struct('prodOnly',         prodOnly, ...
                    'prodIfVector',     prodIfVector, ...
                    'dimRefType',       []);

% Check if "dimVar" represents an empty variable.
givenFlag = (~isempty(dimVar) && numVar ~= 0);
if( ~givenFlag )
    indRef = ones(1, numRef);
    return;
end

%% Sanity check.
% Assert: "dimVar & dimRef" are rows.
assert( isrow(dimVar), '"dimVar" must be a row.' );
assert( isrow(dimRef), '"dimRef" must be a row.' );

% Assert: "dimVar & dimRef" must be non-negative since they represent dimensions.
assert( all(dimVar >= 0), '"dimVar" [%s] must be non-negative.', num2str(dimVar, '%d ') );
assert( all(dimRef >= 0), '"dimRef" [%s] must be non-negative.', num2str(dimRef, '%d ') );

%% Decide the represented type of "dimRef".
dimRefL = length(dimRef);
assert( dimRefL > 0, 'The length of "dimRef" (%d) must be > 0.', dimRefL );

% Check the type.
switch dimRefL
    case 1
        % CASE: It's a scalar.
        dimRefType = "scalar";
        
    case 2
        % CASE: It's a vector or multi-D. Check the values.
        if( any(dimRef == 1) )
            dimRefType = "vector";
        else
            dimRefType = "multi-D";
        end
        
    otherwise
        % CASE: It's a multi-D.
        dimRefType = "multi-D";
end

% Record the value.
multiDInfo.dimRefType = dimRefType;

%% Standardize "dimRef" following the options.
% PS. We use "dimFinal" for variable separation.
dimFinal = dimRef;

% CASE: Consider "numRef" only if "dimRef" represents a vector.
if( prodIfVector && dimRefType == "vector" )
    dimFinal = numRef;
end

% CASE: Consider "numRef" only.
if( prodOnly )
    dimFinal = numRef;
end

%% Assert & generate the multi-D index.
if( length(dimFinal) == 1 )
    % CASE: scalar. Distinguish the "constant & variant" cases.
    if( dimFinal == 1 )
        % CASE: the reference is a constant.
        indRef = ones(1, numVar);
    else
        % CASE: the reference is variant (i.e., element-wise)
        assert( numVar == dimFinal, '#dimVar (%d) must be the same as #dimRef (%d).', numVar, dimFinal );
        indRef = [1:1:numVar];
    end

else
    % CASE: vector. 
    assert( all(dimVar == dimFinal), ...
            'dimVar [%s] & dimRef [%s] should be the same.', ...
            num2str(dimVar, '%d '), num2str(dimFinal, '%d ') );
    indRef = [1:1:numRef];
    
end

end

