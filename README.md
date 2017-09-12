# FastTV-L0-Smoothing-Final-Version
Fast L0 Smoothing Final Version
C++ Version your can set residue term and ADMM or Penalty Method.
Matlab Version has two function:
[S,E]=Gradient_L0smoothing_ADMM(Im,omega,mu,pattern,residue) ADMM Algorithm of L0 Smoothing
[S,E]=Gradient_L0smoothing_Penalty_Method(Im,omega,pattern,residue) Penalty Method of L0 Smoothing
S is output image, E is energy set.
Im is input image, omega is parameter of \omega, mu is start parameter of start \mu
pattern: 1: print out energy and store in E
         2: print out image of each step
         otherwise ignore energy and reconstruction of each step (to accelerate)
residue: mod 1: residue norm is L21
         mod 2: residue norm is L22
