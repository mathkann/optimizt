%
 M = [M  q;-q' 0];
 [m,n] = size(M);
 if exist('toler') ~= 1 
   toler=1.e-6; 
 end
 if exist('gamma') ~= 1 
   gamma=1/sqrt(n);    
 end
 if exist('alpha') ~= 1 
   alpha=0.9;    
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
   dx=[XX*M+SS]\[gamma*mu*ones(n,1)-x.*s-(1-gamma)*XX*(M*x-s)];
   ds=M*dx+(1-gamma)*(M*x-s);
%
%  choose step-size
%
   nora = min([dx./x;ds./s]);
%
%  Check if feasibility can be reached early.
%
   if abs(nora) <= (1-gamma),
     nora = (1-gamma);
     s = s + ds/nora;
     x = x + dx/nora;
     M=M(1:n-1,1:n-1);
     if x(n) > toler^2,
       disp('Find a maximal feasible solution');
       s = s(1:n-1)/x(n);
       x = x(1:n-1)/x(n);
     else
       disp('The problem is infeasible');
       x = x(1:n-1);
       s = s(1:n-1);
     end;
     return;
   else
     nora = abs(alpha/nora);
     x = x + nora*dx;
     s = s + nora*ds;
     mu = x'*s/n;
     gap = mu,
   end;
 end;
%
% output solution
%
 n = n-1;
 M=M(1:n,1:n);
 if x(n+1) > toler^2,
   s=s(1:n)/x(n+1);
   x=x(1:n)/x(n+1);
   disp('Find an approximate complemetarity or maximal-feasible solution');
 else
   x=x(1:n)/s(n+1);
   s=s(1:n)/s(n+1);
   disp('The problem is infeasible');
 end;
 return
%
%  This program solves the monotone LCP.
%
%  See splcp.m
%
%  set eta = 1-gamma and gamma < 1.
%



