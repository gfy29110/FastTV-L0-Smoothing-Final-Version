Im = imread('pflower.jpg');
[S_ADMM,E1] = Gradient_L0smoothing_ADMM(Im, 200, 200, 1, 2);
[S_Penalty_Method,E2] = Gradient_L0smoothing_Penalty_Method(Im, 200, 1, 2);
N1 = size(E1, 2);
N2 = size(E2, 2);
X1 = 0:N1-1;
X2 = 0:N2-1;
figure;
plot(X1, E1, '-*');
title('ADMM');
xlabel('step');
ylabel('energy');
figure;
plot(X2, E2, '-*');
title('Penalty Method');
xlabel('step');
ylabel('energy');