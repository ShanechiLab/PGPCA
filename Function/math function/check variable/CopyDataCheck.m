function [maxX, copySW] = CopyDataCheck(x)
%{
This function decides if some indexes need to be copied.

For example:
There are 3 data with x = [1,5,1] samples, respectively. In this case, the 1st and the 3rd data should make 5 copies to
  be matched to the 2nd data.

NOTE:
1.  All indexes must be positive integers since they represent # samples.
2.  Only one "x" component can > 1 or all compoenets are equal. Otherwise, the copy is not well-defined. 
      (e.g., how to copy 2 to match 5?)

Inputs:
--------------- Essential
x:              a row [1, numX]. The # samples in each dataset.

Outputs:
--------------- Essential
maxX:           a positive scalar. The maximal number of samples.
copySW:         a boolean row [1, numX]. "true" means the corresponding data needs to be copied.
%}

%% Check input integrity.
numX = length(x);
assert( all(x > 0) && all(mod(x,1) == 0), 'All # data must be positive integers.' );

%% Determine the data needs to be copied.
maxX = max(x);
maxInd = find(x == maxX);
if( all(x == maxX) )
    % Special case: no one needs to be copied.
    copySW = false(1, numX);
else
    % Make sure the other values are ones.
    restInd = [1:1:numX];
    restInd(maxInd) = [];
    assert( all(x(restInd) == 1), 'Some small # samples are not ones. This is invalid.' );
    % Prepare the copy.
    copySW = true(1, numX);
    copySW(maxInd) = false;
end

end

