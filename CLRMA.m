%% this function denoise the hyperspectral Raman data by SVD and Statistical analysis
function [Reconstructed, Index_selected, num] = CLRMA(data, reference, C,img_region, windowsize, count, thr1, thr2)
%% parameters 
% data: 2D matrix of the signal 
% reference: 2D matrix of the reference high snr data
% C: the columns of the imaging area
% img_region : select the wavenumber region to determine the smooth area
% windowsize: the window size to search the image to find local smooth area, 5 is recommended
% count: how many smooth area were searching, 100 as default
% thr1: determine how many svd components are in considerate, 1e-4 as default
% thr2: determine which svd components is selected, 1e-2 as default
% Reconstructed: the denoised signal using svd denoising
% Index_selected: the selected SVD component index
% num: how many top svd components are examined in the current process
%% config the default parameters
switch nargin
    case 3
        img_region = 1280:1290;
        windowsize = 5;
        count = 100;
        thr1 = 1e-4;
        thr2 = 1e-2;
    case 2
        img_region = 1280:1290;
        C = 400;
        windowsize = 5;
        count = 100;
        thr1 = 1e-4;
        thr2 = 1e-2;
    case 8
        disp(' ');
    otherwise
        error("Need 8 parameters: data,columns,imaging region, window size, count, thr1, thr2, but only %d is given!\n", nargin);
end
%% centralization
mcube = [data, reference];
MS = mean(mcube,2);
mcube = mcube - MS;
% data = data - MS;

%% find smooth local area on the msemap
[D,L] = size(data);
nim = mean(data(img_region,:),1);
nim = reshape(nim,C,[])';
if nargout~=1
    f1 = figure; figure(f1); subplot 221;
    imagesc(nim); colormap jet; title('Map');set(gca,'XTick',[],'YTick',[]);
end
gap = windowsize -1;
[R,~] = size(nim);
msemap = zeros(R-gap,C-gap); 
for i = 1:R-gap
    for j = 1:C-gap
        patch = nim(i:i+gap,j:j+gap);
        mv =  mean(mean(patch));
        patch = patch - mv;
        patch = patch.^2;
        patch = patch/(mv^2);
        msemap(i,j) = mean(mean(patch));
    end
end
% show the msemap
if nargout ~= 1
    figure(f1);subplot 223;
    imagesc(msemap);colormap jet; title('Mse Map');set(gca,'XTick',[],'YTick',[]);
end
vec = reshape(msemap,1,[]);
vec = sort(vec);
X = [];
Y = [];
for i = 1:count
    [x,y] = find(msemap == vec(i));
    X = [X,x(end)];
    Y = [Y,y(end)];
end
if nargout ~= 1
    figure(f1);subplot 223; hold on; 
    scatter(Y(1:ceil(count/10):end),X(1:ceil(count/10):end),80,'red','s','filled');
end
X = X + gap/2;
Y = Y + gap/2;
% show the anchors position on the image
if nargout ~= 1
    figure(f1);subplot 221; hold on;
    scatter(Y(1:ceil(count/10):end),X(1:ceil(count/10):end),80,'red','s','filled');
end
% Selected  smooth area
for i = 1:count
    sx(i,:) = X(i)-gap/2:X(i)+gap/2;
    sy(i,:) = Y(i)-gap/2:Y(i)+gap/2;
end
cube = reshape(data,D,C,[]); cube = permute(cube,[3 2 1]);
%%
if D<100
    M = 3;
else 
    M = 9;
end
for i = 1:count
    rec_spec(i,:) = cube2spec(cube,sx(i,:),sy(i,:)); % sx, sy区间的平均谱
    rec_spec(i,:) = wden(rec_spec(i,:),'sqtwolog','h','mln',5,'sym8');
%     rec_spec(i,:) = wdenoise(rec_spec(i,:),M,'DenoisingMethod','BlockJS');
end
%% SVD and stastical analysis
% determine how many top SVD components are in consideration
[U,S,V] = svd(mcube, 'econ');V = V'; V = V(:,1:L);
s = S(S~=0); ds = diff(log(s),2); ds = smooth(ds,9);
try
    [index, ~] = find(abs(ds)<=thr1); num = index(1)+2; 
catch
    num = length(s);
end

snrs = []; recon = repmat(MS, 1, L); 
for i = 1:num
    tmp = s(i)*U(:,i)*V(i,:);
    recon = recon + tmp;
    tmp_cube = reshape(recon,D,C,[]);
    tmp_cube = permute(tmp_cube,[3 2 1]);
    for k = 1:count
        newspec = cube2spec(tmp_cube,sx(k,:),sy(k,:));
        residual = rec_spec(k,:) - newspec';
        snrs(k,i) = snr(rec_spec(k,:),residual);
    end
end
% determine which SVD componets is selected to reconstruct the signal
% according to the mean SNRS matrix
% compute basic snr contribution of the mean spectrum
pad = zeros(count,1); 
for k = 1:count
    residual = rec_spec(k,:) - MS';
    pad(k) = snr(rec_spec(k,:),residual);
end
snrs = [pad, snrs]; dfsnr = [];
% disp(snrs);
for i = 1:length(X)
    dfsnr(i,:) = diff(snrs(i,:));
end
mdfsnr = mean(dfsnr); [Index_selected,~]=find(mdfsnr'>=thr2);
if nargout ~= 1
    figure(f1);
    subplot 222; 
    plot(mdfsnr, 'b*-');xlabel('SVs'); ylabel('SNR Contribution/dB'); title(['Mean curves of ', num2str(count), ' regions']);
    subplot 224; 
    plot(dfsnr'); xlabel('SVs'); ylabel('SNR Contribution/dB'); 
end
ss = zeros(size(S)); 
for j = 1:length(Index_selected)
    ss(Index_selected(j),Index_selected(j))= s(Index_selected(j));
end
Reconstructed = U*ss*V + MS;
end