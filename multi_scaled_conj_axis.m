function X=multi_scaled_conj_axis(n,h,m)
% m is the number of times h is repeated in a point
g=general_conj_axis(n,n);
X=[];
for i=1:1:n
    p=g;
    p(:,i)=p(:,i)*h;
    X=vertcat(X,p);
end
