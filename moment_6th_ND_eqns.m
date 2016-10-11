function f=moment_6th_ND_eqns(x,n)
%there are 6 variables r1,r2,r3,w1,w2,w3
%first solve the last 3 eqns for ws and thenfirst 3 eqns for r
% global n
%n is the dimension of the system
r1=x(1);
r2=x(2);
r3=x(3);
if n<=6
w=[(7-(n-1))/r1^6;1/(2^n*r2^6);1/(2*r3^6)];
f=[2*r1^2,2^n*r2^2,4*(n-1)*r3^2;2*r1^4,2^n*r2^4,4*(n-1)*r3^4;0,2^n*r2^4,4*r3^4]*w-[1;3;1];
end
if n>=7 
    w=[(14-n)/(2*r1^6);(n-5)/(2^n*(n-3)*r2^6);1/(4*(n-3)*r3^6)];
    f=[2*r1^2,2^n*r2^2,4*(n-1)*(n-2)*r3^2;2*r1^4,2^n*r2^4,4*(n-1)*(n-2)*r3^4;0,2^n*r2^4,8*(n-2)*r3^4]*w-[1;3;1];
end
end