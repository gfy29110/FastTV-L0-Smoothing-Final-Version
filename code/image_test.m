Im = imread('pflower.jpg');
[S_ADMM,E1] = Gradient_L0smoothing_ADMM(Im, 200, 200, 2, 2);
[S_Penalty_Method,E2] = Gradient_L0smoothing_Penalty_Method(Im, 200, 2, 2);