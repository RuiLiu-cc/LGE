function [xs] = softThresh(x,lambda)
% compute the soft thresholding of the value or matrix  with parameter
% lambda
s = abs(x) + 1e-6;
ss = s - lambda;
ss = ss.*(ss>0);
ss = ss./s;
xs = ss.*x;

end