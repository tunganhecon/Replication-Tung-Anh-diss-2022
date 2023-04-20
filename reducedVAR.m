function [B,SIGMA,U,V] = reducedVAR(dependent_reg, X_reg, t, q)

B=inv(X_reg'*X_reg)*X_reg'*dependent_reg;
U=dependent_reg-B*X_reg;
SIGMA=U*U'/(t-24-24*q-1);	
V=B(:,1);
B=B(:,2:q*24+1);


