%
% Set parameters
%
 [m,n] = size(A);
 n = n+1;
 if exist('toler') ~= 1 
   toler=1.e-6; 
 end
 if exist('gamma') ~= 1 
   gamma=.5;    
 end
 if exist('alpha') ~= 1 
   alpha=0.7;    
 end
%
% Set initial points
%
 x = ones(n,1);
 s = ones(n,1);
 y = zeros(m,1);
 mu1= n;
 mu2= 1;
 gap = 1;
 iter =0;
 dr0  =1.e10;
%
% Start the loop
%
 while gap >= toler,
   iter = iter + 1;
   %
   % Get objective gradient and hessian
   %
   xx = x(1:n-1)/x(n);
   qq = gradien(xx,c);
   MM = hessian(xx,c);
   %
   % Form the homoginized gradient and hessian
   %
   MM = [MM qq-MM*xx;-qq'-xx'*MM xx'*MM*xx];
   qq = [x(n)*qq-A'*y-s(1:n-1);-qq'*x(1:n-1)+b'*y-s(n)];
   %
   % Solving one Newton step with the augmented system
   %
   XX=sparse(1:n,1:n,x/mu2);
   SS=sparse(1:n,1:n,s/mu2);
   gamma = min((mu2/mu1)^2,.5); % Choose centering parameter
   %
   % Check dual feasibility residual
   dr = norm(qq);
   if (dr > 2*dr0)
       %
       % If dual feasibility is not improving, set Newton step to only imprve
       % dual equation while keeping primal feasibility and dualty gap unchanged
       gamma=0;
       qqq = [ones(n,1)-x.*s/mu2-XX*qq;sparse(m,1)];
   else
       % Otherwise, set Newton step to improve all of them at same time
       qqq = [gamma*ones(n,1)-x.*s/mu2-(1-gamma)*XX*qq;-(1-gamma)*(A*x(1:n-1)-b*x(n))];
   end;
   dr0=dr;
   %
   % Linear system solve
   dx=[XX*MM+SS -XX*[A -b]';[A -b] sparse(m,m)]\qqq;
   %
   % Construct primal and dual steps
   dy=dx(n+1:n+m);
   dx=dx(1:n);
   ds=MM*dx-[A -b]'*dy+(1-gamma)*qq;
%
%  choose step-size
%
   nora = min([dx./x;ds./s]);
   nora = abs(alpha/nora);
%
% Update iterates
%
   x = x + nora*dx;
   s = s + nora*ds;
   y = y + nora*dy;
%
% Recompute duality gap
%
   mu1 = mu2;
   mu2 = x'*s/n;
   gap = mu2,
 end;
%
% Output solution or infeasibility certificate
%
 iter
 n = n-1;
 if s(n+1) < x(n+1),
   y=y/x(n+1);
   x=x(1:n)/x(n+1);
   s=gradien(x,c)-A'*y;
   s=max(sparse(n,1),s);
   disp('Find a complementarity solution');
 else
   x=x(1:n)/s(n+1);
   y=y/s(n+1);
   s=s(1:n)/s(n+1);
   if (b'*y > .5);
      disp('The problem is (near) infeasible');
   else
      disp('The problem is (near) unbounded');
   end;
 end;
 return
%
%  This program solves the linearly constrained convex programming.
%
%      min   f(x)
%      s.t.  Ax = b, x >= 0.
%
%  Input 
%      A: Sparse constraint matrix.
%      b: Right-hand column vector
%      toler: relative stopping tolerance: the objective value close to 
%             the optimal one in the range of tolerance. 
%             Default value: 1.e-6.
%      alpha: step size: 0 < alpha < 1. Default value: .9.
%      gamma: weight parameter: 0< gamma <=1. Default value: 1/sqrt(n).
%
%  Matlab functions:
%      gradien(x): return the gradient vector of the objec. function at x,
%      hessian(x): return the hessian matrix of the objec. function at x.
%     
%  Output
%     x>=0  : primal solution: Ax = b,
%     y,s>=0: dual solution: s = gradien(x) - A^Ty,
%             x^Ts =< toler.
%     OR
%     Primal infeasibility certificate:
%             y,s>=0: s=-A^Ty and b^Ty=1 
%             y is an infeasibility certificate for
%             {x: Ax=b, x>=0}
%     OR
%     Primal unbounded certificate:
%             x>=0, Ax=0, gradien(x)'*x=-1
%             x is an infeasibility certificate for
%             {(x>=0,y): gradien(x) - A^Ty>=0}
%  
% For technical details, see
%``On a homogeneous algorithm for the monotone complementarity problem,'' 
% (E. Andersen and Y. Ye), Mathematical Programming 84 (1999) 375-400.




