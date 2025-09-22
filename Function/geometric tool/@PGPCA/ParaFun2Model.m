function [paraStru, funStru] = ParaFun2Model(dimY, dimL, paraStru, funStru, varargin)
%{
10-03-2023: (1) Add option "probYlZType" to include various conditional probability model p(y|z).
            (2) Add output "funStru" with 3 new fields "funValYZOpt, funValTrans, & funProbYlZ".
-----------------------------------------------------------------------------------------------------------

This function transfers the model parameters (maybe from M-step) & functions into two types of modeling
  variables ready-to-use. Some variable examples are below:

1.  Model parameters:   {mainVar, sideVar} -> final noise COV "R = [dimY]".
2.  EM step functions:  "funValYZ & funProbYlZ" which require M-step's outputs.

Still, "generalizability" is the key of this function. Therefore, we don't include any specific class
  as inputs (e.g., CurveRegress) but just functional handles. This guarantees the broad application 
  including future classes (e.g., spline model with basis, other manifold, etc.)
-----------------------------------------------------------------------------------------------------------

Model:  y_t = phi(z_t) + K(z_t)*C*x_t + r_t [COV: R -> Ku(z_t) * diag(mainVar, sideVar) * Ku(z_t)']

where   (A) Ku(z_t) = [K(z_t), Kc(z_t)]. (u: unified basis, c: complement basis).
        (B) Basis "K(z_t) & Kc(z_t)" correspond to variance "mainVar & sideVar", respectively.

NOTE:
1.  Functions phi(z_t) (= funMean) and K(z_t) (= funMatK) are general. No class restriction.
2.  The term "K(z_t)*C*x_t" is not included in model if "C = [] (empty)".
3.  Two types of outputs, parameters & functions, are grouped in two struct "paraStru & funStru". This
      struct-based I/O makes this function generalizable.
-----------------------------------------------------------------------------------------------------------

Output Type 1: Model parameters.

1.  matC    = matC;
2.  vecR    = [mainVar * ones(1, dimL), sideVar * ones(1, dimY - dimL)];
3.  funMatR = @(x) matK = funMatK(x); mainVar * matK * matK' + sideVar * (I_n - matK * matK'); 

NOTE:
1.  "matC" is the same as input. We include it for completeness (including all model parameters).
2.  Since "matC" may be empty, the user must provide "dimL" to avoid glitch.
3.  "funMatR" computation is inspired by the equality In = Ku(z)*Ku(z)' =  K(z)*K(z)' + Kc(z)*Kc(z)'.
4.  The anonymous function cannot have multiple statement, so the real implementation is different.
-----------------------------------------------------------------------------------------------------------

Output Type 2: EM step functions.

1.  funValYZ    = @(y, z) exp( funValYZOpt(y, z) );
2.  funValYZOpt = @(y, z) (-1/2) * { (y - phi(z))' * inv(mvnCov(z)) * (y - phi(z)) };
3.  funValTrans = @(x) log(x);
4.  funProbYlZ  = @(y, z, varStru) mvnpdf(y, phi(z), mvnCov(z)) with parameters in "varStru";

NOTE:
1.  Regularly, function "funProbZ" is independent of model parameters & functions of inputs. Therefore,
      we don't include it here, and the user should provide it.
2.  The real implementation of these functional handles are different from above, which just demonstrates
      the idea. In implementation, efficiency/generalizability/etc. are considered.
3.  We only provide "funProbYlZ", not "funProbY = funProbYlZ * funProbZ", so further manipulation is more
      flexible. This is necessary since the E-step and M-step computes LL differently.
4.  General optimization in E-step: 

    funValYZOpt(z_j^*) = funValTrans( funValYZ(sampY, sampZ) * funProbZ(sampZ) );

    So we need 3 functions: "funValYZ, funValYZOpt, & funValTrans". This can speed up the learning process 
      since some functions can be learned easier after transformation (e.g., ln(exp(x)) = x).
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
dimY:               a +integer. The dimension of observation y_t.
dimL:               a +integer. The rank of K(z), the space of x_t.
paraStru:           a struct. It contains all (learned) system parameters below.
    matC:               a 2D array [dimL, dimC] / empty []. The dimensionality reduction matrix.
    mainVar:            a scalar >= 0. The average noise VAR of the "modeled space".
    sideVar:            a scalar >= 0 / empty []. The average noise VAR of the "unmodeled space".
funStru:            a struct. It collects all function handles in following fields:
    funProbZ:           a func. handle. The probability p(z).
    funMean:            a func. handle. The mean vector "phi(z)" for each "z" sample.
    funMatK:            a func. handle. The "orthonormal" coordinate matrix "K(z)" for each "z" sample.
------------------- Optional ('tag' + 'value' with (number))
"probYlZType":      (1) (probYlZType) a string "mvnpdf". (def: "mvnpdf")
                        It determines the type of conditional probability p(y|z). It affects functional
                          outputs, which are probability-derived.

Outputs:
------------------- Essential
paraStru:           a struct. Besides original fields, we update/add following fields:
    matC:               a 2D array [dimL, dimC] / empty []. We update the original "matC".
    vecR:               a row [1, dimY]. The eigenvalues of covariance "R(z)".
    funMatR:            a func. handle. The covariance "R(z)" given different manifold value.
funStru:            a struct. Besides original fields, we update/add following fields:
    funValYZ            a func. handle. The function derived from p(y|z) for E-step z^* estimation.
    funValYZOpt         a func. handle. The function for MATLAB "fminunc" to learn z^*.
    funValTrans         a func. handle. The function transfering "funValYZ" to speed up "fminunc".
    funProbYlZ          a func. handle. The conditional probability p(y|z).
%}

%% Assign inputs
% Constant - a hidden switch in controlling "funValYZOpt & funValTrans".
funOptType = "SquNorm";

% Parser for optional inputs.
optionP = inputParser;
optionP.addParameter('probYlZType', "mvnpdf", @(x) ismember(x, ["mvnpdf"]));

% Parse the inputs.
optionP.parse(varargin{:});
probYlZType = optionP.Results.probYlZType;

%% Compute: system parameters.
% Sanity check: dimY >= dimL > 0.
assert( dimY >= dimL && dimL > 0, ...
        'The order must follow: dimY (%d) >= dimL (%d) > 0.', dimY, dimL );

% Copy variables for convenience.
matC    = paraStru.matC;
mainVar = paraStru.mainVar;
sideVar = paraStru.sideVar;

% Prepare: vecR.
if( isempty(sideVar) )
    % CASE: the space is unified.
    vecR = mainVar * ones(1, dimY);
else
    % CASE: the space is partitioned.
    vecR = [mainVar * ones(1, dimL), sideVar * ones(1, dimY - dimL)];
end

% We create "paraStruClean & funStruClean" for avoiding unnecessary computation.
paraStruClean = struct('matC',      [], ...
                       'mainVar',   mainVar, ...
                       'sideVar',   sideVar);
funStruClean = struct('funMatK', funStru.funMatK);
% Prepare: funMatR.        
funMatR = @(z) PGPCA.ModelCov('covRPart', z, paraStruClean, funStruClean);

%% Compute: E-step functions
% Select the type of p(y|z).
switch probYlZType
    case "mvnpdf"
        % Prepare: funValYZ.
        funValYZ = @(y, z) PGPCA.ModelPara2Prob("ExpoPart", y, z, paraStru, funStru, ...
                                                "probValForm", "ComplExpo");
        % Select the optimization type.
        switch funOptType
            case "Expo"
                % Prepare: funValTrans.
                funValTrans = @(x) x;
                % Prepare: funValYZOpt.
                funValYZOpt = @(y, z) PGPCA.ModelPara2Prob("ExpoPart", y, z, paraStru, funStru, ...
                                                           "probValForm", "ComplExpo");
            case "SquNorm"                   
                % Prepare: funValTrans.
                funValTrans = @(x) log(x);
                % Prepare: funValYZOpt.
                funValYZOpt = @(y, z) PGPCA.ModelPara2Prob("ExpoPart", y, z, paraStru, funStru, ...
                                                           "probValForm", "SquNorm") * (-1/2);
            otherwise
                error('The selected "funOptType" (%s) is invalid.', funOptType);
        end
        
        % Prepare: funProbYlZ.
        funProbYlZ = @(y, z, varStru) PGPCA.ModelPara2Prob("Probability", y, z, varStru, funStru, ...
                                                           "probValForm", "MixProb");

    otherwise
        error('The given "probYlZType" (%s) is invalid.', probYlZType);
end

%% Prepare outputs.
% Collect: system parameters.
paraStru.matC       = matC;
paraStru.vecR       = vecR;
paraStru.funMatR    = funMatR;

% Collect: system functions.
funStru.funValYZ    = funValYZ;
funStru.funValYZOpt = funValYZOpt;
funStru.funValTrans = funValTrans;
funStru.funProbYlZ  = funProbYlZ;

end

