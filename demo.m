close all; clear; clc;
%% loading data
load('target_2s.mat','cube','x');
target = cube; % target data, low SNR
clear cube;
load('ref1.mat','cube'); 
ref = cube; % reference data, high SNR
clear cube;
load('target_20s.mat','cube');
long_integration = cube;
clear cube;
%% despiking
target = despike2(target);
ref = despike2(ref);
long_integration = despike2(long_integration);
%% compute cosine distance of the mean spectra
T = mean(target, 2);
T = T - mean(T);
T = T/norm(T);
R = mean(ref, 2);
R = R - mean(R);
R = R/norm(R);
dist = 1 - pdist([T, R]', 'cosine');
%% comparison of CLRMA and ALRMA
[recon1, I1, N1] = CLRMA(target, ref, 400, 1:1337, 5, 100, 1e-5, 0.05);
[recon2, I2, N2] = ALRMA(target,400, 1:1337, 5, 100, 1e-5, 0.05);
% ALRMA for long integration time data as the ground truth
[GT, ~, ~] = ALRMA(long_integration,400, 1:1337, 5, 100, 1e-5, 0);

%% results on images
idx = 1286;
im1 = recon1(idx,:); im1 = reshape(im1, 400, [])'; im1 = im1(5:end-5,80:end-80); 
lmx1=min(min(im1)); hmx1=max(max(im1)); im1 = (im1-lmx1)/(hmx1-lmx1);

im2 = recon2(idx,:);  im2 = reshape(im2, 400, [])'; im2 = im2(5:end-5,80:end-80); 
lmx2=min(min(im2)); hmx2=max(max(im2)); im2 = (im2-lmx2)/(hmx2-lmx2);

im3 = GT(idx,:);  im3 = reshape(im3, 400, [])'; im3 = im3(5:end-5,80:end-80);
lmx3=min(min(im3)); hmx3=max(max(im3)); im3 = (im3-lmx3)/(hmx3-lmx3);
ssim1 = ssim(im1, im3); ssim2 = ssim(im2, im3);

im4 = target(idx,:);  im4 = reshape(im4, 400, [])'; im4 = im4(5:end-5,80:end-80);
lmx4=min(min(im4)); hmx4=max(max(im4)); im4= (im4-lmx4)/(hmx4-lmx4);
ssim4 = ssim(im4, im3);

figure;
subplot(221); imagesc(im1); colormap jet; title(['CLRMA,ssim=', num2str(ssim1)]);
subplot(222); imagesc(im2); colormap jet; title(['ALRMA,ssim=', num2str(ssim2)]);
subplot(223); imagesc(im4); colormap jet; title(['Noisy,ssim=', num2str(ssim4)]);
subplot(224); imagesc(im3); colormap jet; title('Groud truth, long exposure time: 20s/line');

%% results on spectra
nn = 8172;
spec_1 = recon1(:,nn); 
spec_2 = recon2(:,nn); 
spec_3 = GT(:,nn); 
spec_4 = target(:,nn);

figure;
subplot(221);
plot(x,spec_1,'r','LineWidth',1); title('CLRMA');
xlabel('Wavenumber/cm^-^1');ylabel('Intensity'); axis tight;
subplot(222); 
plot(x,spec_2,'g','LineWidth',1);title('ALRMA'); 
xlabel('Wavenumber/cm^-^1');ylabel('Intensity'); axis tight;
subplot(223);
plot(x,spec_3,'b','LineWidth',1);  title('Ground truth, long exposure time: 20s/line');
xlabel('Wavenumber/cm^-^1');ylabel('Intensity'); axis tight;
subplot(224); 
plot(x,spec_4,'black','LineWidth',1); title('Noisy'); 
xlabel('Wavenumber/cm^-^1');ylabel('Intensity'); axis tight;

%%
P = 0; Q = 0; L = 0;
for i = 1:size(recon1,2)
    [F,p] = freq_spectrum(recon1(:,i));
    [~,q] = freq_spectrum(recon2(:,i));
    [~,l] = freq_spectrum(target(:,i));
    P = P + p;
    Q = Q + q;
    L = L + l;
end
P = P/i; Q = Q/i; L = L/i;

snum = floor(0.1*length(F));
SP = sum(P(1:snum)); IP = sum(P) - SP; SIRP = SP/IP;
SQ = sum(Q(1:snum)); IQ = sum(Q) - SQ; SIRQ = SQ/IQ;
SL = sum(L(1:snum)); IL = sum(L) - SL; SIRL = SL/IL;

figure;
subplot(2,2,1); plot(F,L); title(['Raw, SIR= ',num2str(SIRL)]);
subplot(2,2,2); plot(F,P); title(['CLRMA, SIR= ',num2str(SIRP)]);
subplot(2,2,3); plot(F,Q); title(['ALRMA, SIR= ',num2str(SIRQ)]);
