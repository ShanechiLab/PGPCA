function xTF = IsNumericProper(xVal)
%{
This function is a simple extension of MATLAB function "isnumeric", which only checks if "xVal's class" is
  numeric, but xVal may be "NaN or Inf", which is not proper in most math operations!

Idea: xTF = true iff xVal is numeric and not NaN or Inf.
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
xVal:               a n-D array. We check its elements are "properly numeric" or not.

Outputs:
------------------- Essential
xTF:                a n-D boolean array. Each [i,j] = true means xVal(i,j) is properly numeric.
%}

xTF = isnumeric(xVal) & ~isnan(xVal) & ~isinf(xVal);

end

