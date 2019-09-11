function [ S ] = FFT_Reconstruct( I,S_x,S_y, omega, method, i )
[M,N,D] = size(I);
sizeI2D = [M,N];
win_x = [1,-1];
win_y = [1;-1];
otfFx = psf2otf(win_x,sizeI2D);
otfFy = psf2otf(win_y,sizeI2D);
Normin1 = fft2(I);
Denormin2 = abs(otfFx).^2 + abs(otfFy ).^2;
if D>1
    Denormin2 = repmat(Denormin2,[1,1,D]);
end
Denormin   = 1 + omega*Denormin2;
h = [S_x, I(:,1,:) - I(:,end,:)];
v = [S_y; I(1,:,:) - I(end,:,:)];
Normin2 = [h(:,end,:) - h(:, 1,:), -diff(h,1,2)];
Normin2 = Normin2 + [v(end,:,:) - v(1, :,:); -diff(v,1,1)];
FS = (Normin1 + omega*fft2(Normin2))./Denormin;
S = real(ifft2(FS));
figure, imshow(S);
str = [method, '\','omega=',num2str(omega),',', 'step=', num2str(i),'.jpg'];
imwrite(S ,str);

end

