# CUTpoints
Conjugate Unscented Transform points

The matlab files provided generate the CUT points. Please see MAIN.m to generate the points of moment orders 4,6 and 8 for the dimensions described in the script. Alternatively, the points are already generated and saved in .mat files(loaded directly into matlab) in the folder "CUT points- mat files". The Gaussian CUT points are in the folder "CUTpoints/CUT points- mat files/gaussian_mat/CUT/". The .mat files are named as "cut6_7D_gaussian.mat" which imples -->  upto 6th order moments in 7 dimensions.   

The .mat files for unform distribution have not been added. But they can be generated using the function
[X,w]=uniform_sigma_pts(bnd_low,bnd_lower,N)
where N is the moment order, like N=2,4,6,8. bnd_low is lower corner of hypercube and bnd_upper is the diagonally opposite upper corner like for example in 3D they are [-2,-2,-2] and [4,4,4] 


- Nagavenkat Adurthi
- Puneet Singla 
- Tarunraj Singh
