%
%  This program solves the linearly constrained convex programming.
%
%      min   f(x)
%      s.t.  Ax = b, x >= 0.
%
%  Input 
%      x0: interior (positive) primal feasile solution
%      y0: interior Lagrange multipliers, i.e, nabla f(x0)-A'*y0 > 0
%      A: Sparse constraint matrix.
%      b: Right-hand column vector
%      toler: relative stopping tolerance: the objective value close to 
%             the optimal one in the range of tolerance. 
%             Default value: 1.e-6.
%      alpha: step size: alpha > 0. Default value: 1.
%      gamma: weight parameter: 0< gamma <=1. Default value: .7.
%
%  Matlab functions:
%      gradien(x): return the gradient vector of the objec. function at x,
%      hessian(x): return the hessian matrix of the objec. function at x.
%     
%  Output
%     x>=0  : primal solution: Ax = b,
%     y,s>=0: Lagrange (dual) solution: s = nabla f(x) - A^Ty,
%             x^Ts/n =< toler.
%  
%  How to call the solver: type splccp and return.
%

%
% Set parameters
%
 [m,n] = size(A);
 if exist('toler') ~= 1 
   toler=1.e-6; 
 end
 if exist('gamma') ~= 1 
   gamma=.5;    
 end
 if exist('alpha') ~= 1 
   alpha=1;    
 end
%
% Set initial points
%
 x=x0;
 y=y0;
 s=gradien(x,c)-A'*y;
 mu = x'*s/n;
 iter =0;
%
% Start the loop
%
 while mu >= toler,
   iter = iter + 1;
   %
   % Get the objective hessian matrix
   %
   MM = hessian(x,c);
   %
   % Form the right-hand-vector
   %
   qqq = [(gamma*mu)*ones(n,1)-x.*s;sparse(m,1)];
   %
   % Sove the system of linear equartions
   %
   XX=sparse(1:n,1:n,x);
   SS=sparse(1:n,1:n,s);
   dx=[XX*MM+SS -XX*A';A sparse(m,m)]\qqq;
   %
   % Construct primal and dual directions
   %
   dy=dx(n+1:n+m);
   dx=dx(1:n);
%
%  choose step-size
%
%   nora = min([dx./x;ds./s]);
%   nora = abs(alpha/nora);
%
% Update iterates
%
   x = x + alpha*dx;
   y = y + alpha*dy;
   s = gradien(x,c)-A'*y;
   if min([x;s])<0,
      disp('No longer feasible; gamma may be too small.');
   end;
%
% Recompute complementarity gap
%
   mu = x'*s/n
 end;
%
% Output solution
%
 iter
 return




