function [paraStruColl, resColl] = ModelParaHD2LD(paraStru, varargin)
%{
This function computes system parameter pair {C, mainVar, sideVar} from a "high-D" case to "low-D" cases.
  It's good for coding integrity by having a function handling this transformation generally.

Model:  y_t = phi(z_t) + K(z_t)*C*x_t + r_t [COV: R -> Ku(z_t) * diag(mainVar, sideVar) * Ku(z_t)']

where   (A) Ku(z_t) = [K(z_t), Kc(z_t)]. (u: unified basis, c: complement basis).
        (B) Basis "K(z_t) & Kc(z_t)" correspond to variance "mainVar & sideVar", respectively.
-----------------------------------------------------------------------------------------------------------

We handle the following cases:

1.  matC == []:     It's error because this is already the smallest case (dimC = 0). 
2.  matC ~= []:     Then the reduced dimension can only within [0:1:dimC].

NOTE:
1.  This low-D pair may not be the optimal one for the corresponding dimension. For example, We transform
      7D -> 5D, but the optimal 5D parameter pair may be something else (unless we can prove it!).
2.  The transformed dimension can only "decreasing" not "increasing". The reason is that we don't know 
      which axis in the ball(mainVar) is the most critical, which is only known from observation samples.
3.  "sideVar" is just copied in this transformation since it's unrelated to the PGPCA EM learning. 
4.  In output "paraStruColl", we don't concatenate "paraStru" directly so its struct is still the same as
      the input one. This makes extracting them in further analyses easier.
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
paraStru:           a struct. It contains all (learned) system parameters below.
    matC:               a 2D array [dimL, dimC] / empty []. The dimensionality reduction matrix.
    mainVar:            a scalar >= 0. The average noise VAR of the "modeled space".
    sideVar:            a scalar >= 0 / empty []. The average noise VAR of the "unmodeled space".
------------------- Optional ('tag' + 'value' with (number))
"lowD":             (1) (lowD) a row [1, numLowD]. (def: [0:1:dimC-1])
                        The low-D dimensions need to be computed from "paraStru".
"concaHDSW":        (1) (concaHDSW) a boolean. (def: true)
                        "true" means input high-D "paraStru" will be concatenated at the end of output
                          "paraStruAll". It's useful in data integrity in further analysis.

Outputs:
------------------- Essential
paraStruColl:       a struct row [1, numLowD (+1)]. It collects low-D parameter sets in following fields:
    dimL:               a +integer. It's rank(K), the dimension of the "main" space.
    dimC:               a +integer. It's rank(C), the dimension of the "freely modeled main" space by x_t.
    paraStru:           a struct. It contains fields "matC, mainVar, & sideVar", the same as input. It's
                          the degenerative system paramter corresponding to "dimC".
resColl:            a struct row [1, numLowD (+1)]. It collects all related results in following fields:
    lowD:               a row [1, numLowD]. The selected low-D dimension.
    concaHDSW:          a boolean. It's the same as the option "concaHDSW".
    valCSqu:            a row [1, dimC]. The squared values of columns of "matC".
    eigVecGa:           a 2D array [dimL, dimC]. Every column is an eigenvector of "matGa".
    eigValGa:           a row [1, dimC]. Each value is an eigenvalue of "matGa".
%}

%% Regularize input.
% Constant.
errRatio = PGPCA.errRatio;

% Extract "paraStru" fields for simplicity.
matC    = paraStru.matC;
mainVar = paraStru.mainVar;
sideVar = paraStru.sideVar;

% Assert: "mainVar" >= 0.
assert( mainVar >= 0, '"mainVar" (%2.2e) must be >= 0.', mainVar );

% Assert: "matC" is nonempty.
assert( ~isempty(matC), '"matC" is empty (dimC = 0). This is the lowest dimension already.' );

% Constant: dimC & dimL.
[dimL, dimC] = size(matC);
assert( dimC <= dimL, '"dimC" (%d) must be <= "dimL" (%d).', dimC, dimL );

%% Assign options.
% Parser for optional inputs.
optionP = inputParser;
optionP.addParameter( 'lowD',       [0:1:dimC-1],   @(x) isrow(x) && min(x) >= 0 && max(x) <= dimC );
optionP.addParameter( 'concaHDSW',  true,           @islogical );

% Parse the inputs.
optionP.parse(varargin{:});
lowD      = optionP.Results.lowD;
concaHDSW = optionP.Results.concaHDSW;

%% Compute: the eigenvalues & eigenvectors of "matGa".
%{
Since "matC'*matC = diag(eig(matGa) - mainVar)", we must recover "matGa" and then compute low-D cases.
%}

% Top eigenvectors of "matGa".
valCSqu = diag( matC' * matC )';
eigVecGa = matC ./ sqrt(valCSqu);

% Top eigenvalues of "matGa".
eigValGa = valCSqu + mainVar;

% Assert: "eigValGa" is decreasing.
assert( issorted(eigValGa, 'descend'), 'The eigenvalues of "matGa" must be decreasing.' );

% Assert: columns of "matC" are mutually orthogonal.
eigVecGaSqu = eigVecGa' * eigVecGa;
maxInnProd = max( abs(eigVecGaSqu - diag(diag(eigVecGaSqu))), [], "all" );
assert( maxInnProd <= errRatio, ...
        'Inner product of eigenvectors of "matGa" (max: %2.2e) must be zero (errTol: %2.2e)', ...
        maxInnProd, errRatio );
    
% The "total" variance in the "dimL - dimC" subspace.
mainVarTotal = mainVar * (dimL - dimC);

%% Compute: high-D -> low-D one by one.
numLowD = length(lowD);

% Initialize the collecting variable.
paraStruColl = struct('dimL', cell(1, numLowD), 'dimC', [], 'paraStru', []);

% Compute low-D cases one by one.
for(j = 1:1:numLowD)
    seleLowD = lowD(j);
    
    % Compute "lowMainVar".
    if(seleLowD == dimL)
        % This is a singular case. There is no unmodeled subspace, and the main variance is zero.
        lowMainVar = 0;
    else
        lowMainVar = (sum(eigValGa(seleLowD+1:1:end)) + mainVarTotal) / (dimL - seleLowD);
    end
    
    % Compute "lowValCSqu".
    lowValCSqu = eigValGa(1:1:seleLowD) - lowMainVar;
    assert( all(lowValCSqu >= 0), ...
            'The low-D (%d) "matC" square (min: %2.2e) must >= 0.', seleLowD, min(lowValCSqu) );
        
    % Compute "lowMatC".
    lowMatC = eigVecGa(:, 1:1:seleLowD) .* sqrt(lowValCSqu);
    
    % Collect variables.
    paraStruLow = struct('matC', lowMatC, 'mainVar', lowMainVar, 'sideVar', sideVar);
    paraStruColl(j) = struct('dimL', dimL, 'dimC', seleLowD, 'paraStru', paraStruLow);
end

%% Prepare output.
% (If asked) attach the input "paraStru" at the end of "paraStruColl".
if( concaHDSW )
    paraStruColl(end + 1) = struct('dimL', dimL, 'dimC', dimC, 'paraStru', paraStru);
end

% Collect other information.
resColl = struct('lowD',        lowD, ...
                 'concaHDSW',   concaHDSW, ...
                 'valCSqu',     valCSqu, ...
                 'eigVecGa',    eigVecGa, ...
                 'eigValGa',    eigValGa);

end

