function [givenFlag, tIndRef] = CheckDimAndTimeIndex( varX, dimInd, dimRef, tInd, numStep )
%{
This function examines the multi-dimensional variable (double/cell/struct/etc.) in two steps:
1.  Check if the indexed dimensions (dimInd) are the same as the reference (dimRef).
2.  Determine the variable is time-variant (numStep) or invariant (1). 
    (A) if time-variant:    tIndRef = [1:numStep].
    (B) if time-invariant:  tIndRef = ones(1,numStep).

NOTE:
1.  This function is very common while handling time-variant/invariant cases, so I package this function.
2.  We also deal with the "empty" case because it's a common "special" mark in many functions/situations.
3.  In this function, "empty" implies tIndRef = ones(1,numStep).

Inputs:
------------------- Essential
varX:               a multi-dimensional variable / empty. The variable whos dimensions need to be checked.
dimInd:             a row [1, numDimInd] / empty. The indexes of dimensions of "varX" to be checked.
dimRef:             a row [1, numDimInd] / a cell row {1, numCase} / empty. The reference/answer of checked dimensions.
                      When it's a row, every element is [1, numDimInd], a possible dimension.
tInd:               a +integer. The index of time.
numStep:            a +integer. The number of step. This is essential for time indexing reference.

Outputs:
------------------- Essential
givenFlag:          a boolean. "true" means "varX" is given (i.e., non-empty).
tIndRef:            a row [1, numSamp]. It depends on if "varX" is time-variant or invariant.
%}

%% Computation.
givenFlag = ~isempty(varX);

if( ~givenFlag )
    % No further check is necessary. However, we set tIndRef = ones(1, numStep) because it's useful in many cases.
    tIndRef = ones(1, numStep);
    return;
else
    % Check the dimension.
    if( ~isempty(dimInd) )
        dimVarX = size(varX, dimInd);
        % Make sure the dimension reference is given.
        assert( ~isempty(dimRef), 'The dimension index is given, so its reference cannot be empty.' );
        % Scan all dimension cases.
        dimRefClass = class(dimRef);
        switch dimRefClass
            case 'double'
                % CASE: only one reference.
                dimCheckRes = all(dimVarX == dimRef);
            case 'cell'
                % CASE: Multiple references.
                numCase = length(dimRef);
                dimCheckRes = false(1, numCase);
                for( j = 1:1:numCase )
                    dimCheckRes(j) = all(dimVarX == dimRef{j});
                end
            otherwise
                error('The class of "dimRef" (%s) must be either "double" or "cell".', dimRefClass);
        end
        % Check if one of the dimension matches.
        assert( any(dimCheckRes), 'The size of variable (%s) is invalid.', num2str(dimVarX, '%d, ') );
    end
    % Decide the time index reference.
    if( ~isempty(tInd) )
        tVarX = size(varX, tInd);
        switch tVarX
            case 1
                % CASE: time-invariant.
                tIndRef = ones(1, numStep);
            case numStep
                % CASE: time-variant.
                tIndRef = [1:1:numStep];
            otherwise
                error('The time dimension (%d) should be either 1 or %d.', tVarX, numStep);
        end
    else
        % In this case, the user doesn't want to know time info., so we set the output empty.
        tIndRef = [];
    end
end

end

