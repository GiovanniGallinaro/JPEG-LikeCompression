clearvars
close all

M = [5  78.55  15.01  3.40  1.59  20.75 
 15  89.76  13.14  3.66  1.83  24.70
 30  103.13  11.44  4.18  2.10  27.01 
 50  117.14  10.07  4.63  2.38  28.75 
 65  126.18  9.35  5.21  2.57  29.70
 80  133.44  8.84  5.45  2.71  30.42
 100  142.56  8.27  6.18  2.90  31.19];

J = [5  14.9  0.02  21.41  
 15  31.8  0.04  25.36 
 30  49.7  0.04  27.73  
 50  67.4  0.04  29.59   
 65  82.7  0.06  30.96 
 80  112  0.06  33.30];

quality = M(:, 1);
size = M(:, 2);
ratio = M(:, 3);
time = M(:, 4);
bpp = M(:, 5);
psnr = M(:, 6);

jpsnr = J(:, 4);
temp_psnr = psnr(1:end-1);

figure
plot(quality, psnr, 'o');
xticks(quality); yticks(psnr);
xlim([0 110]); ylim([19 32.5]);
ylabel("PSNR [dB]"); xlabel("Quality Factor");
grid
% saveas(gcf, "psnr.png")

figure
plot(quality, ratio, 'o');
xticks(quality); yticks(flip(ratio));
xlim([0 110]); ylim([7.5 16]);
ylabel("Compression ratio"); xlabel("Quality Factor");
grid
% saveas(gcf, "ratio.png")

figure
plot(quality(1:end-1), temp_psnr, 'o');
hold on
plot(quality(1:end-1), jpsnr, 'o');
hold off
xticks(quality(1:end-1)); yticks(sort(jpsnr));
xlim([0 95]); ylim([19 35]);
ylabel("PSNR [dB]"); xlabel("Quality Factor");
grid
legend('JPEG-like', 'Standard JPEG', 'Location', 'northwest')
saveas(gcf, "jpegcomp.png")

