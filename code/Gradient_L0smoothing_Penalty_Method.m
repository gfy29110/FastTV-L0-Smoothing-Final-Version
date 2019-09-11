function [ S,E ] = Gradient_L0smoothing_Penalty_Method( Im,omega,pattern,residue )
%FastTV in L_0 norm
t0=clock;
%% Default Value
if nargin==3
    pattern = 0;
    residue = 2;
end
if ~exist('omega','var')
    omega = 1e2;
end
if ~exist('mu','var')
    mu = 1e2;
end

%% Deal with Input
I = im2double(Im);
[M,N,D] = size(I);
win_x = [1,-1];
win_y = [1;-1];
I_x = zeros(M,N-1,D);
I_y = zeros(M-1,N,D);
for k=1:D
    I_x(:,:,k) = conv2(I(:,:,k),win_x,'valid');
    I_y(:,:,k) = conv2(I(:,:,k),win_y,'valid');
end

%% Assignment for Starter
S_x = I_x;
S_y = I_y;
Z_x = I_x;
Z_y = I_y;
E=[];

%% Memory Pre-allocated

%% Error Parameter
mu = 0.04;
mu_max = 1e6;
i = 0;

%% Energy

if(pattern == 1)
    Energy = L0Norm(I_x,I_y,M,N);
    E = [E,Energy];
    fprintf('%3dth step, Energy is %f\n',i,Energy);
end


%% ADMM
while  (mu < mu_max)%tentative
    %Fix S,U Update Z
    ZA = S_x(1:end-1,:,:);
    ZB = S_y(:,1:end-1,:);
    ZC = S_x(2:end,:,:) - S_y(:,2:end,:);
    Z_judge = (ZA.^2 + ZB.^2 + ZC.^2 - (ZA - ZB - ZC).^2./3) > 2/mu;
    Z_x(1:end-1,:,:) = Z_judge.*(2.*ZA + ZB + ZC)./3;
    Z_y(:,1:end-1,:) = Z_judge.*(ZA + 2.*ZB - ZC)./3;
    Z_x(end,:,:) = S_x(end,:,:);
    Z_y(:,end,:) = S_y(:,end,:);
    
    %Fix U,Z Update S
    if(residue == 1)
        r = 3;
    else
        r = 1;
    end
    for t = 1:r
        if(residue == 1)
            NS = omega./sqrt((S_x(2:end,:,:) - I_x(2:end,:,:)).^2 + (S_y(:,2:end,:) - I_y(:,2:end,:)).^2 + 1e-7);
        else
            NS = omega;
        end
        NS1 = NS + mu;
        NS2 = NS + 2*mu;
        NS3 = NS + 3*mu;
        MS = NS1.*NS3;
        SA = NS.*I_x(2:end,:,:) + mu*(Z_x(2:end,:,:));
        SB = NS.*I_y(:,2:end,:) + mu*(Z_y(:,2:end,:));
        SC = mu*(Z_x(1:end-1,:,:) - Z_y(:,1:end-1,:));
        S_x(2:end,:,:) = (SA.*NS2 + SB*mu + SC.*NS1)./(MS);
        S_y(:,2:end,:) = (SA*mu + SB.*NS2 - SC.*NS1)./(MS);
        S_x(1,:,:) = (Z_x(1,:,:));
        S_y(:,1,:) = (Z_y(:,1,:));
    end
    
    i = i + 1;
    %Energy
    if(pattern == 1)
        E1 = L0Norm(S_x,S_y,M,N);
        if(residue == 1)
            E2 = sum(sum(sum(sqrt((I_x(1:end-1,:,:) - S_x(1:end-1,:,:)).^2 + (I_y(:,1:end-1,:) - S_y(:,1:end-1,:)).^2))))*omega/2;
        else
            E2 = sum(sum(sum((I_x(1:end-1,:,:) - S_x(1:end-1,:,:)).^2 + (I_y(:,1:end-1,:) - S_y(:,1:end-1,:)).^2)))*omega/2;
        end
        Energy = E1 + E2;
        E = [E,Energy];
        fprintf('%3dth step, Energy is %f, E1:%f, E2:%f\n',i,Energy,E1,E2);
    else
        fprintf('.');
    end
    mu = mu * 2;
    
    if(pattern == 2)
        
        FFT_Reconstruct(I, S_x, S_y, omega, 'Penalty_Method', i);
    end
    
end

%% New Method to Reconsturct Image (FFT)
S = FFT_Reconstruct(I, S_x, S_y, omega, 'Penalty_Method', i);
%% 
fprintf('\n');
T = etime(clock,t0);
fprintf('Cost of time is %f\n',T);
end