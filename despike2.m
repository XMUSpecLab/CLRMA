%% desipike the data cube
function despiked = despike2(data,W,windowsize,thr)
%% parameters:
% D: spectral depth
% W: width
% data: 2D matrix
% thr: threshold 
% windowsize : medfilt window size, 3 is recommended
%% config parameters:
switch nargin
    case 1
        W = 400;
        windowsize = 3;
        thr = 1e-6;
    case 2
        windowsize = 3;
        thr = 1e-6;
    case 3
        error("4 parameters is needed, data, column pixels, window size, thr, but %d is given", nargin);
end

D = size(data,1);
cube = reshape(data,D,W,[]);
cube = permute(cube,[3,2,1]);
for i = 1:D
    medcube = medfilt2(cube(:,:,i),[windowsize,windowsize],'symmetric');
    err =  cube(:,:,i)- medcube;
    ee = reshape(err,[],1);
    % 分位数
    m = quantile(ee,ceil(0.99*length(ee)));
    ee(ee>=m(end)|ee<=m(1)) = [];
    pd = fitdist(ee,'Normal');
    lw = icdf(pd,thr);
    hw = icdf(pd,1-thr);
    [x,y] = find(err>=hw|err<=lw);
    if ~isempty(x)
        fprintf('第%d个波长共有%d个鬼峰\n',i,length(x));
        cube(x,y,i) = medcube(x,y);
    end
    clear x y;
end
cube = permute(cube,[3,2,1]);
despiked = reshape(cube,D,[]);
end