function probYlZ = PairCondiProb(sampY, sampZ, funProbYlZ)
%{
This AUXILIARY function computes the "paired" conditional probability. It's good for coding integrity.

Function:   probYlZ(j) = funProbYlZ(sampY(j,:), sampZ(j,:));  

NOTE:
1.  This auxiliary function should be very simple. It just packages the paired code in E/M-step.
2.  Since "probYlZ" is a pairwise result, we give inputs to "funProbYlZ" one by one.
3.  To pair "sampY & sampZ", they must have the same #sample (numY == numZ).
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
sampY:              a 2D array [numY, dimY]. Each row [j,:] is a sample of y.
sampZ:              a 2D array [numZ, dimZ]. Each row [j,:] is a sample of z.
funProbYlZ:         a func. handle. It computes the conditional probability p(y|z).

Outputs:
------------------- Essential
probYlZ:            a column [numY, 1]. The conditional probability p(y_j|z_j).
%}

%% Progress

% Assert: #sample must be the same.
numY = size(sampY, 1);
numZ = size(sampZ, 1);
assert( numY == numZ, '"numY" (%d) must = "numZ" (%d).', numY, numZ );

% Compute p(y|z) one by one.
probYlZ = zeros(numY, 1);
for(j = 1:1:numY)
    probYlZ(j) = funProbYlZ(sampY(j,:), sampZ(j,:));
end

end

