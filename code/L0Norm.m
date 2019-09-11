function [ norm ] = L0Norm( S_x,S_y,M,N )
X_logic = abs(S_x(1:M-1,:,:)) <= 0.001;
Y_logic = abs(S_y(:,1:N-1,:)) <= 0.001;
X_Logic = prod(X_logic,3);
Y_Logic = prod(Y_logic,3);
Logic = ~(X_Logic.*Y_Logic);
norm = sum(sum(Logic));
end

