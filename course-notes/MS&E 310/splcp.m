%
 M = [M  q;-q' 0];
 [m,n] = size(M);
 if exist('toler') ~= 1 
   toler=1.e-6; 
 end
 if exist('gamma') ~= 1 
   gamma=n/(1+n);    
 end
 if exist('alpha') ~= 1 
   alpha=0.7;    
 end
 x = ones(n,1);
 s = ones(n,1);
 mu= 1;
 gap = 1;
 iter =0;
 while gap >= toler,
%
   iter = iter + 1;
   XX=sparse(1:n,1:n,x);
   SS=sparse(1:n,1:n,s);
   dx=[XX*M+SS]\[gamma*mu*ones(n,1)-(1-gamma)*x.*s];
   ds=M*dx;
%
%  choose step-size
%
   nora = min([dx;ds]);
   if nora/mu > -toler,
     if dx(n) > toler/n,
       s=max(ds(1:n-1),zeros(n-1,1))/abs(dx(n));
       x=max(dx(1:n-1),zeros(n-1,1))/abs(dx(n));
       M=M(1:n-1,1:n-1);
       disp('Find an approximate maximal feasible solution');
       return;
     elseif dx(1:n-1)'*ds(1:n-1) < toler,
       x=max(dx(1:n-1),zeros(n-1,1))/abs(ds(n));
       s=max(ds(1:n-1),zeros(n-1,1))/abs(ds(n));
       M=M(1:n-1,1:n-1);
       disp('The problem is infeasible');
       return;
     end;
   end;
   nora = min([dx./x;ds./s]);
   nora = abs(alpha/nora);
   x = x + nora*dx;
   s = s + nora*ds;
   mu = x'*s/n;
   gap = 1/mu,
 end;
%
%  output solution
%
 M=M(1:n-1,1:n-1);
 if x(n)/s(n) > toler,
   s=s(1:n-1)/x(n);
   x=x(1:n-1)/x(n);
   disp('Find an approximate complemetarity or maximal-feasible solution');
 else
   x=x(1:n-1)/s(n);
   s=s(1:n-1)/s(n);
   disp('The problem is infeasible');
 end;
 return
%
%  This program solves the LCP in standard form:
%
%     s = Mx + q, x^Ts=0, (x,s) >= 0.
%
%  Input 
%      M: Sparse n-dimensional monotone matrix.
%      q: n-dimensional column vector
%      toler: relative stopping tolerance: the objective value close to 
%             the local optimal one in the range of tolerance. 
%             Default value: 1.e-6.
%      alpha: step size: 0 < alpha < 1. Default value: .7.
%      gamma: weight parameter: 0< gamma <=1. Default value: n/(n+1).
%  Output
%     x,s: a maximal feasible or complementarity solution.
%
%  Homogeneous self-dual for LCP
%
%  Comment: Use eta =0 and gamma > 1.
%




