function ceVal = CrossEntropy(ceMode, probP, probQ, varargin)
%{
This function computes the crossed entropy "-p(i)*log(q(i))" while taking care the extreme case 0*log(0).
  To be general, there are two modes:

1.  ceMode = "element":     Only compute the element-wise value "-p(i)*log(q(i))". No summation.

    (A) "probP & probQ" arrays have the same size for element-wise computation. 
    (B) Their values must be properly numeric >= 0 (i.e., no NaN or Inf)!
    (C) We define 0*log(0) = 0. Note that if q(i) = 0, then p(i) must be 0, too (but not conversely).

2.  ceMode = "complete":    The cross entropy H(p,q) = sum_{i=1:numSamp}[ -p(i)*log(q(i)) ].

    (A) All conditions and element-wise computations in "element" form are needed.
    (B) For summation, "probP & probQ" must be 2D arrays with "numSamp = size(probP, 2)". 
    (C) sum(probP, 2) & sum(probQ, 2) must be ones, so they are probability.

NOTE:
1.  We include the "element" mode since we only need "-p(i)*log(q(i))" sometimes (like EM log-likelihood).
2.  MATLAB's function "crossentropy" is for neural network cost function, so I write my own function.
3.  Due to numerical precision limitation, sum(probP) & sum(probQ) may not be exactly one. Therefore, we
      have a tolerance range to compensate it.
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
ceMode:             a string "element / complete". The mode of computing the crossed entropy.
probP:              a n-D array / 2D array [numCase, numSamp]. Each element is a probability p([i,j,...]).
probQ:              a n-D array / 2D array [numCase, numSamp]. Each element is a probability q([i,j,...]).
------------------- Optional ('tag' + 'value' with (number))
"negaSW":           (1) (negaSW) a boolean. (def: false)
                        "true" means "ceVal" is multiplied by (-1) (= p(i)*log(q(i))). This is useful in
                          some computations (e.g., EM log-likelihood).
"sumCheckSW":       (1) (sumCheckSW) a boolean row [1,2]. (def: [true, true])
                        "sumCheckSW = [checkProbP, checkProbQ]" controls if the conditions sum(probP) = 1
                          and sum(probQ) = 1 is applied.

Outputs:
------------------- Essential
ceVal:              a n-D array / a column [numCase, 1]. The element-wise/complete cross-entropy.
%}

%% Assign inputs.
% Constant.
dimProbP = size(probP);
dimProbQ = size(probQ);
% Constant - the numerical tolerance.
errTol = 1e-10;

% Parser for optional inputs.
optionP = inputParser;
optionP.addParameter('negaSW',      false,      @islogical);
optionP.addParameter('sumCheckSW',  true(1,2),  @(x) islogical(x) && isvector(x) && length(x) == 2);

% Parse the inputs.
optionP.parse(varargin{:});
negaSW      = optionP.Results.negaSW;
sumCheckSW  = optionP.Results.sumCheckSW;

%% Sanity check: probP & probQ.
% Assert: probP & probQ have the same size.
assert( all(dimProbP == dimProbQ), ...
        'Sizes of "probP" [%s] & "probQ" [%s] must be the same.', num2str(dimProbP), num2str(dimProbQ) );

% Assert: probP & probQ are properly numerical.
assert( all(IsNumericProper(probP), 'all'), '"probP" must be properly numerical!' );
assert( all(IsNumericProper(probQ), 'all'), '"probQ" must be properly numerical!' );

% Assert: probP & probQ are >= 0.
assert( all(probP >= 0, 'all'), '"probP" must be >=0.' );
assert( all(probQ >= 0, 'all'), '"probQ" must be >=0.' );

%% Compute: element-wise cross-entropy.
% Assert: if q(i) = 0, p(i) = 0. Otherwise, cross-entropy is not well-defined.
zeroIndP = find(probP == 0);
zeroIndQ = find(probQ == 0);
assert( all(ismember(zeroIndQ, zeroIndP)), ...
        'Some probQ(i) = 0 (Index: %s) but corresponding probP(i) ~= 0.', ...
        num2str(setdiff(zeroIndQ, zeroIndP)) );

% Compute element-wise cross-entropy & substitute the "0*log(0) = 0" elements.
ceVal = probP .* log(probQ) * (-1);
ceVal(zeroIndQ) = 0;

% Assert: all values are proper (a sanity check just in case).
assert( all(IsNumericProper(ceVal), 'all'), '"ceVal" must be properly numerical!' );

%% (If asked) Compute: complete cross-entropy.
switch ceMode
    case "element"
        % CASE: the element-wise cross-entropy.
        % Do nothing.

    case "complete"
        % CASE: the complete cross-entropy.
        % Assert: probP & probQ are 2D arrays.
        assert( length(dimProbP) == 2, ...
                'When "ceMode = complete", the #dim of "probP & probQ" (%d) must be 2.', ...
                length(dimProbP) );

        % Assert: probP & probQ are row-based "real PMF".
        if( sumCheckSW(1) )
            errProbSum = abs(sum(probP, 2) - 1);
            assert( all( errProbSum < errTol ), ...
                    'All row-based sum of "probP" must = 1 (max error: %2.2e, errTol: %2.2e).', ...
                    max(errProbSum), errTol );
        end
        
        if( sumCheckSW(2) )
            errProbSum = abs(sum(probQ, 2) - 1);
            assert( all( errProbSum < errTol ), ...
                    'All row-based sum of "probQ" must = 1 (max error: %2.2e, errTol: %2.2e).', ...
                    max(errProbSum), errTol );
        end
        
        % Compute: complete cross-entropy.
        ceVal = sum(ceVal, 2);
        
    otherwise
        error('The given "ceMode" (%s) is invalid.', ceMode);
end

%% Prepare outputs.
if( negaSW )
    ceVal = ceVal * (-1);
end

end

