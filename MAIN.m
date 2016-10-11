
n=6;  % dimension
P=eye(n); % desired covariance matrix
mu=zeros(n,1); % desired mean

%% 4th order points
[X4,w4]=conjugate_dir_gausspts_4thmoments(mu,P);


%% 6th order points (dimension 2 to 9)
if (n>=2 && n<=9)
    [X6,w6]=conjugate_dir_gausspts_6moment(mu,P);
end

%% 8th order points (dimension 2 to 6)
if (n>=2 && n<=6)
    [X8,w8]=conjugate_dir_gausspts_8moment(mu,P);
end