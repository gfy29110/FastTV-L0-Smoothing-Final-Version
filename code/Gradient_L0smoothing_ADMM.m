function [ S,E ] = Gradient_L0smoothing_ADMM( Im,omega,mu,pattern,residue )
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
U_x = zeros(M,N-1,D);
U_y = zeros(M-1,N,D);
U_t = zeros(M-1,N-1,D);
E=[];

%% Memory Pre-allocated

%% Error Parameter
Dete_R = 1.1;              
Dete_S = 1.1; 
Dete_Z_x = zeros(M,N-1,D);
Dete_Z_y = zeros(M-1,N,D);
R_Stop_Iter = 1.0; 
S_Stop_Iter = 1.0; 
i = 0;

%% Energy

if(pattern == 1)
    Energy = L0Norm(I_x,I_y,M,N);
    E = [E,Energy];
    fprintf('%3dth step, Energy is %f\n',i,Energy);
end


%% ADMM
while  (Dete_R > R_Stop_Iter||Dete_S > S_Stop_Iter) &&(i<10)%tentative
    %Fix S,U Update Z
    ZA = S_x(1:end-1,:,:) + U_x(1:end-1,:,:);
    ZB = S_y(:,1:end-1,:) + U_y(:,1:end-1,:);
    ZC = S_x(2:end,:,:) - S_y(:,2:end,:) + U_t;
    Z_judge = (ZA.^2 + ZB.^2 + ZC.^2 - (ZA - ZB - ZC).^2./3) > 2/mu;
    Z_x(1:end-1,:,:) = Z_judge.*(2.*ZA + ZB + ZC)./3;
    Z_y(:,1:end-1,:) = Z_judge.*(ZA + 2.*ZB - ZC)./3;
    Z_x(end,:,:) = S_x(end,:,:) + U_x(end,:,:);
    Z_y(:,end,:) = S_y(:,end,:) + U_y(:,end,:);
    
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
        SA = NS.*I_x(2:end,:,:) + mu*(Z_x(2:end,:,:) - U_x(2:end,:,:));
        SB = NS.*I_y(:,2:end,:) + mu*(Z_y(:,2:end,:) - U_y(:,2:end,:));
        SC = mu*(Z_x(1:end-1,:,:) - Z_y(:,1:end-1,:) - U_t);
        S_x(2:end,:,:) = (SA.*NS2 + SB*mu + SC.*NS1)./(MS);
        S_y(:,2:end,:) = (SA*mu + SB.*NS2 - SC.*NS1)./(MS);
        S_x(1,:,:) = (Z_x(1,:,:) - U_x(1,:,:));
        S_y(:,1,:) = (Z_y(:,1,:) - U_y(:,1,:));
    end
    
    %Fix S,Z Update U
    R_x = S_x - Z_x;
    R_y = S_y - Z_y;
    R_t = S_x(2:M,:,:) - S_y(:,2:N,:) -Z_x(1:M-1,:,:) + Z_y(:,1:N-1,:);
    U_x = U_x + R_x;
    U_y = U_y + R_y;
    U_t = U_t + R_t;
    
    %Update Stop 
    i = i + 1;
    Dete_R = sqrt(sum(sum(sum(R_x.^2))) + sum(sum(sum(R_y.^2))) + sum(sum(sum(R_t.^2))));
    Dete_S = sqrt(sum(sum(sum((Dete_Z_x).^2))) + sum(sum(sum((Dete_Z_y).^2))));
    R_max = max(sqrt(sum(sum(sum(S_x.^2)))+sum(sum(sum(S_y.^2)))+sum(sum(sum((S_x(2:M,:,:) - S_y(:,2:N,:)).^2)))),sqrt(sum(sum(sum(Z_x.^2)))+sum(sum(sum(Z_y.^2)))+sum(sum(sum((Z_x(1:M-1,:,:)-Z_y(:,1:N-1,:)).^2)))));
    S_max = sum(sum(sum(U_x.^2)))+sum(sum(sum(U_y.^2)))+sum(sum(sum(U_t.^2)));
    R_Stop_Iter = 1e-5*sqrt(3*M*N-2*M-2*N+1) + 1e-5*R_max;
    S_Stop_Iter = 1e-5*sqrt(M*N) + 1e-5*S_max;
    %fprintf('%dth step, R:%f/%f\n',i,Dete_R,R_Stop_Iter);
    
    %Energy
    if(pattern == 1)
        E1 = L0Norm(S_x,S_y,M,N);
        if(residue == 1)
            E2 = sum(sum(sum(sqrt((I_x(1:end-1,:,:) - S_x(1:end-1,:,:)).^2 + (I_y(:,1:end-1,:) - S_y(:,1:end-1,:)).^2))))*omega/2;
        else
            E2 = sum(sum(sum((I_x(1:end-1,:,:) - S_x(1:end-1,:,:)).^2 + (I_y(:,1:end-1,:) - S_y(:,1:end-1,:)).^2)))*omega/2;
        end
        Energy = E1 + E2;
        E = [E,Energy]; %#ok<AGROW>
        fprintf('%3dth step, Energy is %f, E1:%f, E2:%f\n',i,Energy,E1,E2);
    else
        fprintf('.');
    end
    %Update mu
    if Dete_R > 10*Dete_S
        mu = mu*2;
    end
    if Dete_S > 10*Dete_R
        mu = mu/2;
    end
    
    %save part
    if(pattern == 2)
        FFT_Reconstruct(I, S_x, S_y, omega, 'ADMM', i);
    end
end

%% New Method to Reconsturct Image (FFT)
S = FFT_Reconstruct(I, S_x, S_y, omega, 'ADMM', i);
%% 
fprintf('\n');
T = etime(clock,t0);
fprintf('Cost of time is %f\n',T);
end