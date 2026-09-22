function out = bubble_night(rawDir, calFile, outDir, varargin)
%BUBBLE_NIGHT 面向等离子泡团的气辉整夜处理 → 地理画布上的"耗竭率图"
%
%   out = bubble_night(rawDir, calFile, outDir, 'Name', Value, ...)
%
% 与 airglow_geo_pipeline.m 的关系
%   那条线是"看图"用的: 8 bit 拉伸 + 邻帧背景, 结果适合人眼看波状结构,
%   但泡团是"气辉的百分之几~几十的耗竭", 8 bit 会把定量信息压掉。
%   本函数保留 float 域, 并针对泡团重做了背景处理。
%
% 算法次序（每一步都在对付一个实测确认的混淆项）
%   ① 减暗场: 天空 ~5400 DN 里基座占 3441 DN; 不扣则 D=(I-R0)/R0 被基座稀释。
%   ② 投影到 250 km 地理画布, z<=Zmax: 统一空间尺度, 一次做掉方位校正 rot。
%      z 越大 km/px 放大越快(30deg 时 0.75, 70deg 时 4.7, 88deg 时 36 km/px)。
%   ③ 除原始帧的**逐像素时间中位数** R0: 相机固定 ⇒ 渐晕/尘粒/固定亮点/条纹
%      等静态结构被精确除掉, 泡团"出现又消失"故除不掉。
%   ④ **方位中位数径向归一化**: 逐帧把以某点为中心(暗夜=天顶, 月夜=月面)的
%      径向轮廓除掉。用中位数而不是最小二乘拟合基 —— 实测最小二乘基在画布
%      边缘杠杆极强, 会把边缘处的泡团吸掉(注入 25% 只回收 6.6%); 换方位中位数后
%      回收 120%。对"方位局地"的泡团天然免疫, 对"方位对称"的径向轮廓精确抓住。
%   ⑤ 减一次逐像素时间中位数, 去掉残余静态项 → D(x,y,t) 单位 %。
%
% 参数（名值对）
%   'DarkFile'   暗场 PNG 路径（强烈建议给, 见①）
%   'Zmax'       投影限幅天顶角, 默认 70（度）。径向标定在 z<76 有实测锚点,
%                再往外是外推; 且 km/px 在 z>75 后迅速恶化。
%   'ResDeg'     画布分辨率, 默认 0.02 度 ≈ 2.22 km
%   'Step'       隔几帧取一帧, 默认 1
%   'MoonCenters' N×2 矩阵, 每帧月面在画布坐标(度)下的位置; 给了就用它当
%                径向归一化的圆心(月辉是以月面为中心的散射轮廓)。
%                可用 moon_canvas_pos.m 由 cal + 历表算出。
%   'MaskDeg'    圆心处要掩掉的半径(度), 默认 0（月面本身由饱和掩膜处理）
%   'Injector'   函数句柄 @(i, pr, ctx) -> pr。用于**合成泡团正对照**:
%                在第 i 帧投影后、归一化前把 pr 乘上 (1-delta)。
%                ctx 给出画布几何: xv, yv, Xg, Yg, valid, km_per_px, zMax。
%                见 bubble_synth_test.m。
%   'MakeFig'    是否出图, 默认 true
%   'Tag'        输出文件名后缀
%
% 输出 out 结构: dep(T,ny,nx single), times, xv, yv, zDeg, lat, lon, R0, thr
%
% 例
%   calFile = 'E:\qihui2\Processed_Output\calibration_params.mat';
%   out = bubble_night('E:\qihui2\原始数据\20250920', calFile, ...
%       'E:\qihui2\bubble_work\mat_out', ...
%       'DarkFile','E:\qihui2\mooncal_work\from_downloads\ODAZH_DCAI01_DFOA_AUX_STP_20260516115948_V01.00.PNG', ...
%       'Zmax', 70);

p = inputParser;
addParameter(p, 'Pattern', 'ODAZH_DCAI01_AROA_L0_STP_*.PNG');
addParameter(p, 'DarkFile', '');
addParameter(p, 'Zmax', 70);
addParameter(p, 'ResDeg', 0.02);
addParameter(p, 'Step', 1);
addParameter(p, 'MoonCenters', []);
addParameter(p, 'MaskDeg', 0);
addParameter(p, 'Injector', []);   % @(i, pr, ctx) -> pr   合成泡团正对照用
addParameter(p, 'SatLevel', 65000);
addParameter(p, 'NBin', 48);
addParameter(p, 'MakeFig', true);
addParameter(p, 'Tag', '');
addParameter(p, 'Verbose', true);
parse(p, varargin{:});
o = p.Results;

RE = 6371.0; H = 250.0;
S = load(calFile); cal = S.cal;
LAT0 = cal.obs_lat; LON0 = cal.obs_lon;
R_OBS = RE + 0.05; R_AIR = RE + H;

% 夜名（取 rawDir 最后一级目录名，末尾分隔符先去掉）
rd = rawDir;
while ~isempty(rd) && any(rd(end) == ['\' '/'])
    rd(end) = [];
end
[~, nightName] = fileparts(rd);
if isempty(nightName), nightName = 'night'; end

if isempty(o.DarkFile)
    warning('未给暗场: 基座 3441 DN 会把耗竭率稀释一半以上, 结果只能看形态不能看振幅。');
    dark = 0;
else
    dark = q_read_gray(o.DarkFile);    % 按位深统一到 16 bit 量级（与帧同一尺度）
    med = median(dark(:)); sd = std(dark(:));
    nfix = sum(dark(:) > med + 5*sd);
    dark = min(dark, med + 5*sd);      % 暗场热点不修会让减完出现负值, 归一化除零爆掉
    if o.Verbose, fprintf('暗场: 中位 %.0f DN, 修掉 %d 个 >5sigma 热点\n', med, nfix); end
end

% ---------- 画布几何（与 airglow_geo_pipeline 同一套约定）----------
cx = cal.zenith(1); cy = cal.zenith(2);
r_use = polyval(cal.rz_poly, o.Zmax);
gd = gamma_of_z(o.Zmax, R_OBS, R_AIR);
half = gd / sqrt(2) * 0.99;
xv = -half : o.ResDeg : half;
yv = xv;
[Xg, Yg] = meshgrid(xv, yv);
LATg = LAT0 + Yg;
LONg = LON0 + Xg / cosd(LAT0);
sinLat = sind(LATg); cosLat = cosd(LATg);
sinLat0 = sind(LAT0); cosLat0 = cosd(LAT0);
dLon = wrap180(LONg - LON0);
cosG = min(max(sinLat*sinLat0 + cosLat*cosLat0.*cosd(dLon), -1), 1);
Gamma = acos(cosG);
zDeg = atan2d(R_AIR.*sin(Gamma), R_AIR.*cos(Gamma) - R_OBS);
A = mod(atan2d(sind(dLon).*cosLat, cosLat0.*sinLat - sinLat0.*cosLat.*cosd(dLon)), 360);
A(Gamma < 1e-10) = 0;
r = polyval(cal.rz_poly, zDeg);
aa = deg2rad(A + cal.rot);
xs = cx + r .* sin(aa);
ys = cy - r .* cos(aa);
valid = isfinite(xs) & isfinite(ys) & (zDeg <= o.Zmax) & (r <= r_use) & ...
        xs >= 1 & xs <= 1024 & ys >= 1 & ys <= 1024;
[ny, nx] = size(xs);
if o.Verbose
    fprintf('画布 %dx%d | z<=%.0f deg | %.2f km/px | 半宽 %.3f deg (%.0f km)\n', ...
        nx, ny, o.Zmax, o.ResDeg*111.195, half, half*111.195);
end

% ---------- 第一遍: 投影 ----------
d = dir(fullfile(rawDir, o.Pattern));
if isempty(d), d = dir(fullfile(rawDir, '*.PNG')); end
names = {d.name};
files = names(1:o.Step:end);
T = numel(files);
A3 = nan(ny, nx, T, 'single');
times = NaT(T, 1, 'TimeZone', 'UTC');
satfrac = zeros(T, 1);
% 注入器上下文: 让 Injector 拿到画布几何(不必自己重算一套约定)
ctx = struct('xv', xv, 'yv', yv, 'Xg', Xg, 'Yg', Yg, 'valid', valid, ...
    'km_per_px', o.ResDeg * 111.195, 'zMax', o.Zmax, 'ny', ny, 'nx', nx);
nBad = 0;
for i = 1:T
    fp = fullfile(rawDir, files{i});
    im = q_read_gray(fp);                  % 按位深统一到 16 bit 量级
    if isempty(im)
        A3(:, :, i) = NaN;                 % 解码失败：该帧不可用
        satfrac(i) = 1;                    % 按"全饱和"记账，别让它假装正常
        nBad = nBad + 1;
        continue
    end
    sat = im >= o.SatLevel;
    im = im - dark;
    pr = interp2(1:size(im,2), 1:size(im,1), im, xs, ys, 'linear', NaN);
    bad = interp2(1:size(im,2), 1:size(im,1), double(sat), xs, ys, 'linear', 0);
    pr(bad > 0.3) = NaN;
    pr(~valid) = NaN;
    if ~isempty(o.MoonCenters)
        rr = hypot(Xg - o.MoonCenters(i,1), Yg - o.MoonCenters(i,2));
        pr(rr <= o.MaskDeg) = NaN;
    end
    if ~isempty(o.Injector)
        pr = o.Injector(i, pr, ctx);        % 合成泡团注入(正对照)
    end
    A3(:,:,i) = single(pr);
    satfrac(i) = mean(bad(valid) > 0.3);
    tok = regexp(files{i}, '(\d{8})(\d{6})', 'tokens');
    if ~isempty(tok)
        t = tok{1}{1}; times(i) = datetime([t(1:4) '-' t(5:6) '-' t(7:8) ' ' ...
            tok{1}{2}(1:2) ':' tok{1}{2}(3:4) ':' tok{1}{2}(5:6)], 'TimeZone', 'UTC');
    end
end
if o.Verbose
    fprintf('投影完成 %d 帧, 饱和像元率中位 %.3f%%\n', T, 100*median(satfrac));
    if nBad > 0
        fprintf('  ⚠ %d/%d 帧解码失败（已按"不可用"计入饱和率）\n', nBad, T);
    end
end

% ---------- 第二遍: R0 -> 方位中位数径向归一化 ----------
R0 = median(A3, 3, 'omitnan');
floorv = 0.15 * median(R0(valid), 'omitnan');
R0 = max(R0, floorv);
C3 = nan(ny, nx, T, 'single');
for i = 1:T
    c = double(A3(:,:,i)) ./ double(R0);
    if ~isempty(o.MoonCenters)
        ctr = o.MoonCenters(i, :);
    else
        ctr = [0 0];
    end
    c = radial_median_norm(c, Xg, Yg, valid, ctr, o.NBin);
    c(~valid) = NaN;
    C3(:,:,i) = single(c);
end
ref = median(C3, 3, 'omitnan');
dep = (C3 - ref) * 100;
dep(:, ~valid) = NaN;
v = dep(repmat(valid, [1 1 T]));
if o.Verbose
    fprintf('D: 中位 %+.2f%%  sigma %.2f%%  0.1 分位 %+.1f%%\n', ...
        median(v, 'omitnan'), std(v, 'omitnan'), pctl(v, 0.1));
end

out = struct('night', nightName, 'dep', dep, 'times', times, 'xv', xv, 'yv', yv, ...
    'zDeg', zDeg, 'lat', LATg, 'lon', LONg, 'R0', R0, 'valid', valid, ...
    'satfrac', satfrac, 'km_per_px', o.ResDeg*111.195, 'z_max', o.Zmax, ...
    'names', {files});

if nargin >= 3 && ~isempty(outDir)
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    f = fullfile(outDir, ['bubble' o.Tag '_dep.mat']);
    save(f, '-struct', 'out', '-v7.3');
    if o.Verbose, fprintf('已写 -> %s\n', f); end
end

if o.MakeFig && ~isempty(outDir)
    npanel = min(12, T);
    win = min(8, max(1, floor(T / 4)));
    idx = unique(round(linspace(1, T, npanel)));
    ncol = 4; nrow = ceil(numel(idx) / ncol);
    figure('Position', [50 50 420*ncol 340*nrow], 'Visible', 'off');
    for k = 1:numel(idx)
        i = idx(k);
        subplot(nrow, ncol, k);
        m = median(dep(:,:,max(1,i-win):min(T,i+win)), 3, 'omitnan');
        imagesc(xv, yv, flipud(m), [-25 25]); axis xy; colorbar;
        title(sprintf('%s UT  min %.0f%%', datestr(times(i), 'HH:MM'), min(m(:))));
    end
    sgtitle(sprintf('耗竭率图 %s', rawDir), 'Interpreter', 'none');
    f = fullfile(outDir, ['bubble' o.Tag '_montage.png']);
    exportgraphics(gcf, f); close(gcf);
    if o.Verbose, fprintf('已写 -> %s\n', f); end
end
end

% ======================================================================
function g = gamma_of_z(z_deg, R_OBS, R_AIR)
% 地心角(度) <- 天顶角(度): 视线与 H 球壳交点几何（单调, 二分）
lo = 0; hi = acosd(R_OBS / R_AIR) * (1 - 1e-12);
for k = 1:200
    mid = 0.5 * (lo + hi);
    if atan2d(R_AIR*sind(mid), R_AIR*cosd(mid) - R_OBS) < z_deg, lo = mid; else, hi = mid; end
end
g = 0.5 * (lo + hi);
end

function d = wrap180(d)
d = mod(d + 180, 360) - 180;
end

function p = pctl(x, pct)
% 避免依赖 Statistics Toolbox
x = sort(x(isfinite(x)));
if isempty(x), p = NaN; return; end
k = max(1, min(numel(x), round(pct/100 * numel(x))));
p = x(k);
end

function C = radial_median_norm(C, Xg, Yg, valid, ctr, nbin)
% 方位中位数径向归一化: 逐环取中位数当背景, 除之。
% 为什么不用最小二乘拟合基: 见文件头④ —— 拟合基在画布边缘杠杆太强, 会吃掉泡团。
X = Xg - ctr(1); Y = Yg - ctr(2);
rf = hypot(X, Y);
rf = rf / max(rf(:));
bi = min(max(floor(rf * nbin) + 1, 1), nbin);
m0 = valid & isfinite(C);
prof = nan(nbin, 1);
for b = 1:nbin
    mb = m0 & (bi == b);
    if sum(mb(:)) >= 60
        prof(b) = median(C(mb));
    end
end
good = isfinite(prof);
if sum(good) < 3, return; end
prof = interp1(find(good), prof(good), (1:nbin)', 'linear', 'extrap');
C = C ./ prof(bi);
end
