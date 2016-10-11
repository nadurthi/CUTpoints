function [c,f]=moment_4th_ND_eqns(x)
%there are 6 variables r1,r2,r3,w1,w2,w3
%first solve the last 3 eqns for ws and thenfirst 3 eqns for r
%n is the dimension of the system
r1=x(1);
r2=x(2);

f=r1^2+2*r2^2-r1^2*r2^2;
c=[];
end
