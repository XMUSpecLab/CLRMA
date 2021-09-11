%% desipike the data cube
function despiked = despike3(data,W,windowsize,thr)
%% parameters:
% D: spectral depth
% W: width
% data: 2D matrix
% thr: threshold 
% windowsize : medfilt window size, （3，3，3） is recommended
%% config parameters:
switch nargin
    case 1
        W = 400;
        windowsize = 3;
        thr = 1e-10;
    case 2
        windowsize = 3;
        thr = 1e-10;
    case 3
        error("4 parameters is needed, data, column pixels, window size, thr, but %d is given", nargin);
end

D = size(data,1);
cube = reshape(data,D,W,[]);
cube = permute(cube,[3,2,1]);

medcube = medfilt3(cube,[windowsize, windowsize, windowsize], 'symmetric');
err = cube - medcube;
%对err按平均谱进行缩放
scale = cube2spec(medcube,1:size(cube,1),1:W);

for i = 1:D
    err(:,:,i)=err(:,:,i)/scale(i);
end

ee = reshape(err,[],1);
m = quantile(ee,ceil(0.9999*length(ee)));
ee(ee>=m(end)|ee<=m(1)) = [];
pd = fitdist(ee,'Normal');
lw = icdf(pd,thr);
hw = icdf(pd,1-thr);

z = [];
for i = 1:D
    [x,y] = find(err(:,:,i)>=hw|err(:,:,i)<=lw);
    if ~isempty(x)
        fprintf('第%d个波长共有%d个鬼峰\n',i,length(x));
        cube(x,y,i) = medcube(x,y,i);
    end
end
cube = permute(cube,[3,2,1]);
despiked = reshape(cube,D,[]);
end