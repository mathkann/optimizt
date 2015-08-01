function [g]=gradien(x,c);
[nn,mm]=size(x);
g=-c./x;
return
% This is the gradient vector of wieghted barrier function
%       - \sum c_i\log x_i