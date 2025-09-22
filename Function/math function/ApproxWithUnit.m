function valArray = ApproxWithUnit(valArray, varargin)
%{
This function approximates numbers with the given unit (e.g., 20, 0.5, etc.), which is not 1 anymore!

The approximating methods:
1.  round
2.  ceil
3.  floor

Steps:
1.  mod(valArray, unit) to find the closest value given the unit.
2.  Decide the approximated value based on the method.

NOTE:
1.  numVal = dimX1*dimX2*...*dimXn.
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
valArray:           a multi-D array [dimX1, dimX2, ..., dimXn]. The numeric array of values.
------------------- Optional ('tag' + 'value' with (number))
"method":           (1) (approxMd) a scalar/vector/array of cell/string/char [1, numVal]. (def: "round") 
                        The method to approximate values. Vector/array version are element-specific.
"unit":             (1) (approxUnit) a +scalar/vector/array. (def: 1)
                        The unit of module. Vector/array version are element-specific.

Outputs:
------------------- Essential
valArray:           a multi-D array [dimX1, dimX2, ..., dimXn]. The module version from the input.
%}

%% Assign.
% Constant.
dimVal = size(valArray);
numVal = numel(valArray);

% Default.
approxMd = "round";
approxUnit = 1;

% Assign options (if asked).
numVarargin = length(varargin);
if( numVarargin > 0 )
    tagInd = 1;
    while( tagInd <= numVarargin )
        % Select the option.
        tagName = varargin{tagInd};
        switch tagName
            case 'method'
                tagShift = 1;
                approxMd = varargin{tagInd + 1};
                
            case 'unit'
                tagShift = 1;
                approxUnit = varargin{tagInd + 1};
                
            otherwise
                error('The given option "%s" is invalid.', tagName);
        end
        % Shift the index.
        tagInd = tagInd + tagShift + 1;
    end
end

%% Sanity & value normalization.
% Assert: the unit must positive.
assert( all(approxUnit > 0), 'All elements of "approxUnit" must be positive.' );

% Standardize "approxMd" into a scalar/vector/array of "string" (not char or cell!) 
approxMd = string(approxMd);

% Create the indexing set for "approxMd & approxUnit" (constant / element-wise)
[~, indRefMd, ~]    = CheckMultiDimIndex(dimVal, size(approxMd), 'prodIfVector', true);
[~, indRefUnit, ~]  = CheckMultiDimIndex(dimVal, size(approxUnit), 'prodIfVector', true);

%% Module the values.
for( j = 1:1:numVal )
    % Get the value, method, & unit.
    seleVal     = valArray(j);
    seleMd      = approxMd(indRefMd(j));
    seleUnit    = approxUnit(indRefUnit(j));
    
    % Compute the "scale" based on the method.
    valRatio = seleVal/seleUnit;
    switch seleMd
        case "round"
            valRatio = round(valRatio);
        case "ceil"
            valRatio = ceil(valRatio);
        case "floor"
            valRatio = floor(valRatio);
        otherwise
            error('The given approximating method "%s" (element %d) is invalid.', seleMd, j);
    end
    
    % Multiply the unit back.
    valArray(j) = valRatio * seleUnit;
end

end

