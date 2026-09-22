function res = bubble_search(depFile, outDir, varargin)
%BUBBLE_SEARCH  在耗竭率立方(_dep.mat)里搜索"持久、南北向拉长"的耗竭事件
%
%   res = bubble_search('E:\qihui2\bubble_work\bubble_dep.mat', 'E:\qihui2\bubble_work');
%   res = bubble_search(depFile, outDir, 'Thr', 10, 'MinFrames', 5, 'MakeFig', false)
%
% 输入
%   depFile  bubble_night.m 写出的 *_dep.mat（含 dep/times/xv/yv/valid/...）
%   outDir   输出目录：events.csv 与三张图
%
% 为什么判据要这么设（每条都对应一个实测踩过的坑）
% ① **绝对阈值下限 8%，不能用"3σ"**。平滑场的 σ 是**云的空间起伏**（实测 13~39%），
%    不是噪声；拿 3σ 当阈值等于 40~116%，永不判出。
% ② **三维时间连通（(帧,y,x) 26 邻域）且要求 ≥3 帧**。云的单帧亮/暗斑块连不起来；
%    180 s 采样下真泡团每帧只挪约 2.4 px，天然连通。这样"持久性"本身就成了判据。
% ③ **二阶矩给拉长比与主轴方位**。等离子泡团应南北向拉长（主轴方位接近 0/180）。
%
% 参数（名值对）
%   'Thr'        阈值(%)，默认 [] = 用 AbsFloor
%   'AbsFloor'   绝对阈值下限(%)，默认 8
%   'Smooth'     高斯平滑 σ(画布像元)，默认 5
%   'MinArea'    最小面积(画布像元)，默认 1000（≈4900 km²）
%   'MinFrames'  事件最短持续帧数，默认 3
%   'Dt'         帧间隔(s)，默认 [] = 从时间轴**实测**（Step>1 时自动识别，不会算错漂移率）
%   'BandDeg'    keogram/Hovmöller 取的中心窄带半宽(deg)，默认 0.35
%   'MakeFig'    是否出图，默认 true
%   'Verbose'    打印，默认 true
%
% 输出 res: sigma / thr / bandDeg / events(结构体数组) / mask / sm / eventsCsv

p = inputParser;
addParameter(p, 'Thr', []);
addParameter(p, 'AbsFloor', 8);
addParameter(p, 'Smooth', 5);
addParameter(p, 'MinArea', 1000);
addParameter(p, 'MinFrames', 3);
addParameter(p, 'Dt', []);
addParameter(p, 'BandDeg', 0.35);
addParameter(p, 'MakeFig', true);
addParameter(p, 'Verbose', true);
parse(p, varargin{:});
o = p.Results;

S = load(depFile);
if ~isfield(S, 'dep')
    error('bubble_search:noDep', '%s 里没有 dep 变量。', depFile);
end
dep = S.dep;                       % (ny, nx, T) single
valid = logical(S.valid);
xv = S.xv(:)'; yv = S.yv(:)';
[ny, nx, T] = size(dep);
nv = sum(valid(:));
resDeg = double(S.km_per_px) / 111.195;      % 画布 deg/像元
KM = 111.195;
night = get_night(S);
if isfield(S, 'times'), times = S.times; else, times = NaT(T, 1, 'TimeZone', 'UTC'); end

% 帧间隔：默认从时间轴**实测**（隔帧取样时 Dt 不是 180 s，写死会把漂移率算错几倍）
dtSec = o.Dt;
if isempty(dtSec)
    if numel(times) > 1 && ~any(isnat(times))
        dtSec = median(seconds(diff(times)));
    else
        dtSec = 180;
    end
end
if o.Verbose && abs(dtSec - 180) > 1
    fprintf('  实测帧间隔 %.0f s（不是 180 s，说明用了 Step>1）\n', dtSec);
end

if o.Verbose
    fprintf('  [%s] %d 帧 | 画布 %dx%d | z<=%.0f deg | %.3f deg/px (%.2f km/px)\n', ...
        night, T, nx, ny, S.z_max, resDeg, S.km_per_px);
end

% ---------- ① 逐帧高斯平滑（NaN 先填 0，之后用 valid 掩膜剔掉）----------
sm = zeros(ny, nx, T, 'single');
for i = 1:T
    f = double(dep(:, :, i));
    f(~isfinite(f)) = 0;
    sm(:, :, i) = single(imgaussfilt(f, o.Smooth));
end
v3 = repmat(valid, [1 1 T]);
sigma = std(sm(v3), 'omitnan');      % 仅供报告: 这是空间起伏, 不是噪声
thr = o.Thr;
if isempty(thr), thr = o.AbsFloor; end

% ---------- ② 阈值 + 三维连通(26 邻域) ----------
mask = (sm < -thr) & v3;
lab = uint32(bwlabeln(mask, ones(3, 3, 3, 'logical')));
nlab = double(max(lab(:)));

idxAll = find(mask);
% 事件上限：预分配避免增长告警。实测一夜事件数在个位数量级，500 足够。
MAXEV = 500;
evs = repmat(struct('night', '', 'ut0', '', 'ut1', '', 'i0', [], 'i1', [], ...
    'nframe', [], 'dur_min', [], 'area_km2', [], 'x_deg', [], 'y_deg', [], ...
    'elong', [], 'angle_deg', [], 'minD_pct', [], 'drift_px', [], 'drift_ms', []), ...
    MAXEV, 1);
nev = 0;
if nlab > 0
    [yy, xx, tt] = ind2sub([ny nx T], idxAll);
    lk = double(lab(idxAll));
    area = accumarray([lk, tt], 1, [nlab, T]);

    for k = 1:nlab
        fr = find(area(k, :) > 0);
        if numel(fr) < o.MinFrames, continue; end
        medArea = median(area(k, fr));
        if medArea < o.MinArea, continue; end
        big = fr(area(k, fr) >= 0.5 * o.MinArea);
        if isempty(big), continue; end
        if nev >= MAXEV, break; end

        % 只对本标签的体素算子集（省内存）
        sel = (lk == k);
        yy_k = yy(sel); xx_k = xx(sel); tt_k = tt(sel);
        el = []; ang = []; cxs = []; cys = [];
        for q = 1:numel(big)
            b = big(q);
            mb = (tt_k == b);
            if sum(mb) < 4, continue; end
            [e, a2, cyk, cxk] = moments2(yy_k(mb), xx_k(mb));
            el(end+1) = e; ang(end+1) = a2; cys(end+1) = cyk; cxs(end+1) = cxk; %#ok<AGROW>
        end
        if isempty(el), continue; end

        i0 = fr(1); i1 = fr(end);
        if numel(big) >= 3
            pf = polyfit(big(:), cxs(:), 1);
            driftPx = pf(1);
        else
            driftPx = NaN;
        end
        driftMs = driftPx * resDeg * KM * 1000 / dtSec;
        minD = min(sm(idxAll(sel)));

        nev = nev + 1;
        evs(nev) = struct( ...
            'night', night, ...
            'ut0', hhmm(times, i0), ...
            'ut1', hhmm(times, i1), ...
            'i0', i0, 'i1', i1, 'nframe', numel(fr), ...
            'dur_min', (i1 - i0 + 1) * dtSec / 60, ...
            'area_km2', medArea * (resDeg * KM)^2, ...
            'x_deg', xv(1) + (median(cxs) - 1) * resDeg, ...
            'y_deg', yv(1) + (median(cys) - 1) * resDeg, ...
            'elong', median(el), 'angle_deg', median(ang), ...
            'minD_pct', minD, 'drift_px', driftPx, 'drift_ms', driftMs);
    end
    evs = evs(1:nev);
    if nev > 0
        [~, ord] = sort([evs.nframe], 'descend');
        evs = evs(ord);
    end
end

maskFrac = mean(sum(sum(mask, 1), 2) / max(nv, 1));
if o.Verbose
    fprintf(['  空间σ=%.1f%% (云/梯度起伏, 非噪声)  阈值 %.1f%%  ' ...
        '耗竭面积率 %.2f%%  事件 %d 个 (>= %d 帧, >= %d px)\n'], ...
        sigma, thr, 100 * maskFrac, numel(evs), o.MinFrames, o.MinArea);
    for k = 1:min(5, numel(evs))
        e = evs(k);
        fprintf(['      %s~%s UT  %d 帧 (%.0f min)  x=%+.2f y=%+.2f deg  %.0f km2  ' ...
            '拉长 %.2f 主轴 %.0f deg  最深 %.1f%%  漂移 %.0f m/s\n'], ...
            e.ut0, e.ut1, e.nframe, e.dur_min, e.x_deg, e.y_deg, ...
            e.area_km2, e.elong, e.angle_deg, e.minD_pct, e.drift_ms);
    end
end

% ---------- ③ 输出 ----------
res = struct('sigma', sigma, 'thr', thr, 'bandDeg', o.BandDeg, ...
    'events', evs, 'mask', mask, 'sm', sm, 'nv', nv, 'resDeg', resDeg, ...
    'night', night, 'xv', xv, 'yv', yv, 'times', times);

if ~isempty(outDir)
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    csvp = fullfile(outDir, [night '_events.csv']);
    fid = fopen(csvp, 'w', 'n', 'UTF-8');
    fprintf(fid, ['night,ut0,ut1,nframe,dur_min,area_km2,x_deg,y_deg,' ...
        'elong,angle_deg,minD_pct,drift_px,drift_ms\n']);
    for k = 1:numel(evs)
        e = evs(k);
        fprintf(fid, '%s,%s,%s,%d,%.0f,%.0f,%.3f,%.3f,%.2f,%.1f,%.2f,%.2f,%.0f\n', ...
            e.night, e.ut0, e.ut1, e.nframe, e.dur_min, e.area_km2, ...
            e.x_deg, e.y_deg, e.elong, e.angle_deg, e.minD_pct, e.drift_px, e.drift_ms);
    end
    fclose(fid);
    res.eventsCsv = csvp;
    if o.Verbose, fprintf('  事件清单 -> %s  (%d 条)\n', csvp, numel(evs)); end

    if o.MakeFig
        res.figMontage   = fig_montage(res, outDir);
        res.figKeogram   = fig_keogram(res, outDir);
        res.figHovmoller = fig_hovmoller(res, outDir);
        if o.Verbose
            fprintf('  %s\n  %s\n  %s\n', res.figMontage, res.figKeogram, res.figHovmoller);
        end
    end
end
end

% ==========================================================================
function n = get_night(S)
% 夜名：优先用 bubble_night 写下的 night 字段；老文件则从文件名里解析日期
if isfield(S, 'night') && ~isempty(S.night)
    n = S.night;
    return
end
n = 'unknown';
if isfield(S, 'names') && ~isempty(S.names)
    tok = regexp(S.names{1}, '(\d{8})', 'tokens');
    if ~isempty(tok), n = tok{1}{1}; end
end
end

% ==========================================================================
function s = hhmm(times, i)
% NaT 安全的时间格式（datestr(NaT) 会报错）
if isnat(times(i))
    s = sprintf('#%d', i);
else
    s = datestr(times(i), 'HH:MM');
end
end

% ==========================================================================
function [elong, angDeg, cy, cx] = moments2(ys, xs)
% 二阶矩 -> (拉长比, 主轴方位 deg, 质心行, 质心列)
cy = mean(ys); cx = mean(xs);
wy = ys - cy; wx = xs - cx;
C = [mean(wx .* wx), mean(wx .* wy); mean(wx .* wy), mean(wy .* wy)];
[V, D] = eig(C);
[~, j] = max(diag(D));
elong = sqrt(max(D(end, end), 1e-9) / max(D(1, 1), 1e-9));
major = V(:, j);
angDeg = mod(atan2d(major(1), major(2)), 180);
end

% ==========================================================================
function p = fig_montage(res, outDir)
T = size(res.sm, 3);
win = min(8, max(1, floor(T / 4)));
npanel = min(12, T);
idx = unique(round(linspace(max(win, 1), max(T - win, 1), npanel)));
ncol = 4; nrow = ceil(numel(idx) / ncol);
figure('Position', [40 40 400 * ncol 320 * nrow], 'Visible', 'off');
for k = 1:numel(idx)
    i = idx(k);
    lo = max(1, i - win); hi = min(T, i + win);
    m = median(res.sm(:, :, lo:hi), 3, 'omitnan');
    subplot(nrow, ncol, k);
    imagesc(res.xv, res.yv, flipud(m), [-25 25]); axis xy; colorbar;
    hold on;
    contour(res.xv, res.yv, flipud(double(res.mask(:, :, i))), [0.5 0.5], ...
        'g', 'LineWidth', 0.7);
    hold off;
    title(sprintf('%s UT  min %.0f%%', hhmm(res.times, i), min(m(:))), ...
        'FontSize', 8);
    set(gca, 'FontSize', 6);
end
sgtitle(sprintf('%s 逐时段耗竭图（%.1f h 窗；绿线=候选，阈值 -%.1f%%）', ...
    res.night, 2 * win * 3 / 60, res.thr), 'Interpreter', 'none');
p = fullfile(outDir, [res.night '_dep_montage.png']);
exportgraphics(gcf, p); close(gcf);
end

% ==========================================================================
function p = fig_keogram(res, outDir)
% 时间-纬度图: 取 |东向| 中心窄带的南北剖面, 漂移的耗竭带 -> 斜条
sel = abs(res.xv) <= res.bandDeg;
kg = squeeze(median(res.sm(:, sel, :), 2, 'omitnan'));      % (ny, T)
if size(kg, 1) ~= numel(res.yv), kg = kg.'; end
figure('Position', [60 60 900 360], 'Visible', 'off');
th = (0:size(kg, 2) - 1) * 3 / 60;
imagesc(th, res.yv, kg); set(gca, 'YDir', 'normal');
caxis([-25 25]); colorbar; colormap(rdbu());
xlabel('UT 小时（自起始帧）'); ylabel('北向 [deg]');
title(sprintf('%s 时间-纬度图 (|东向|<%.2f deg 中位)', res.night, res.bandDeg), ...
    'FontSize', 10, 'Interpreter', 'none');
p = fullfile(outDir, [res.night '_keogram.png']);
exportgraphics(gcf, p); close(gcf);
end

% ==========================================================================
function p = fig_hovmoller(res, outDir)
% 时间-东向图: 取 |北向| 中心窄带的东西剖面
sel = abs(res.yv) <= res.bandDeg;
kg = squeeze(median(res.sm(sel, :, :), 1, 'omitnan'));      % (nx, T)
if size(kg, 1) ~= numel(res.xv), kg = kg.'; end
th = (0:size(kg, 2) - 1) * 3 / 60;
figure('Position', [60 60 900 360], 'Visible', 'off');
imagesc(th, res.xv, kg); set(gca, 'YDir', 'normal');
colorbar; colormap(rdbu());
xlabel('UT 小时（自起始帧）'); ylabel('东向 [deg]');
title(sprintf('%s 时间-东向图 (|北向|<%.2f deg 中位)', res.night, res.bandDeg), ...
    'FontSize', 10, 'Interpreter', 'none');
p = fullfile(outDir, [res.night '_hovmoller.png']);
exportgraphics(gcf, p); close(gcf);
end

% ==========================================================================
function c = rdbu()
% 蓝-白-红 发散色标（不依赖额外工具箱）
n = 64;
r = [linspace(0.13, 1, n/2), ones(1, n/2)]';
g = [linspace(0.30, 1, n/2), linspace(1, 0.42, n/2)]';
b = [ones(1, n/2), linspace(1, 0.13, n/2)]';
c = [r, g, b];
end
