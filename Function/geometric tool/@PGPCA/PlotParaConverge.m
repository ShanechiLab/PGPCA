function varargout = PlotParaConverge(varargin)
%{
This function analyzes & plots the evolution of system parameters in learning. Since parameter dimensions
  can be high, the main output form is not parameter itself, but ||\hat{para} - true_para||. The parameter
  & modes are:

1.  matC:       error norm.
2.  mainVar:    error norm. para value.
3.  sideVar:    error norm. para value.

NOTE:
1.  If a field in "paraStruColl" is empty, we skip that parameter which is N/A.
2.  All values are plotted in the same figure. Therefore, we suggest "normalized error norm" as the mode.
3.  The normalization fails when the referenced value is zero, and show the warning message (if asked).
4.  In general, "paraStruColl" may not come from "PGPCA EM". It generalizes this function.
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
paraStruRef:        sdf
paraStruColl:       sdf



------------------- Optional ('tag' + 'value' with (number))
"metricForm":

"normSW":



"plotSW":

"plotPara":

"axisH":


"warnMsgSW":


Outputs:
------------------- Essential




%}














end

