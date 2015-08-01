function [h]=hessian(x,c);
[nn,mm]=size(x);
h=sparse(diag(c./(x.*x)));
return
% This is the Hessian matrix of weighted barrier function
%       - \sum c_i*\log x_i