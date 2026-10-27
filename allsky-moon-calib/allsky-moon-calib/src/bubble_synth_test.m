function res = bubble_synth_test(night, varargin)
%BUBBLE_SYNTH_TEST  合成泡团正对照 —— 回答"我们到底能看见多深的泡团"
%
%   res = bubble_synth_test('20250920');
%   res = bubble_synth_test('20250920', 'Amp', [0.25 0.10], 'Step', 40);
%   res = bubble_synth_test('20260326', 'MoonCenters', mc, 'MaskDeg', 0.8);
%
% 原理（受控注入）
%   往**真实帧**里注入一个参数完全已知的合成泡团：南北向拉长的超椭圆
%   （默认半轴 120 x 250 km，边缘软过渡 15 km，幂次 4），以 30 m/s 东向漂移，
%   在夜的中段存在数小时。然后走**同一条处理链**，量它回收了多少。
%   因为注入量已知，"回收率"就把整条链的保真度与噪声水平一次称出来。
%
% 为什么对照必须"配对、同帧"
%   dep 虽然减了逐像元时间中位数，但"整幅同时变亮/变暗"这一模式仍留在里面。
%   若拿**另一时段**的帧当对照，而注入窗恰好覆盖子夜、对照帧落在黄昏/黎明，
%   对照区自己就 +14%，会被误记成回收，报出 >100% 的荒谬回收率。
%   本函数在**同一帧**里取对照：把核区以画布原点（天顶）旋转 180 度，
%   得到"同帧、同半径、同面积、方位相反"的 patch ——
%   逐帧全局模式与径向轮廓被同时消掉，剩下的只有真正的方位局地差异。
%
% 参数（名值对）
%   'Amp'        注入振幅数组（占气辉的比例），默认 [0.25 0.10 0.05]
%   'Akm','Bkm'  超椭圆半轴(km)，默认 120 / 250
%   'EdgeKm'     边缘软过渡(km)，默认 15
%   'P'          超椭圆幂次，默认 4
%   'DriftMs'    东向漂移(m/s)，默认 30
%   'I0','I1'    注入起止帧号（1-based，闭开）；给空则用 Frac 推
%   'Frac'       注入窗占整夜的比例，默认 [0.25 0.70]
%   'X0','Y0'    注入中心初始位置(画布 deg)，默认 -2.40 / -1.20
%   'Step'       隔帧取样，默认 1（调试可给大值）
%   'MoonCenters','MaskDeg'   月夜：逐帧月面画布坐标与月盘掩膜半径，同 bubble_night
%   'MakeFig'    出图，默认 true
%   'OutDir'     输出目录，默认 cfg.bubRoot\synth
%   'Verbose'    打印，默认 true
%
% 输出 res: 结构体数组（每个振幅一条），字段
%   amp, recovered(净回收 %), expected(掩膜内平均注入 %), ratio(回收率),
%   sigma(逐帧散布 %), limit(3σ 单帧极限 %), limitStack(3σ 时窗平均极限 %),
%   det(是否可探测), ncore(核区像元数), areaKm2(核区面积)

p = inputParser;
addParameter(p, 'Amp', [0.25 0.10 0.05]);
addParameter(p, 'Akm', 120);
addParameter(p, 'Bkm', 250);
addParameter(p, 'EdgeKm', 15);
addParameter(p, 'P', 4);
addParameter(p, 'DriftMs', 30);
addParameter(p, 'I0', []);
addParameter(p, 'I1', []);
addParameter(p, 'Frac', [0.25 0.70]);
addParameter(p, 'X0', -2.40);
addParameter(p, 'Y0', -1.20);
addParameter(p, 'Step', 1);
addParameter(p, 'MoonCenters', []);
addParameter(p, 'MaskDeg', 0);
addParameter(p, 'MakeFig', true);
addParameter(p, 'OutDir', '');
addParameter(p, 'Verbose', true);
parse(p, varargin{:});
o = p.Results;

cfg = q_cfg();
KM = 111.195;
if isempty(o.OutDir), o.OutDir = fullfile(cfg.bubRoot, 'synth'); end
if ~exist(o.OutDir, 'dir'), mkdir(o.OutDir); end
rawDir = fullfile(cfg.rawRoot, night);
if ~exist(rawDir, 'dir')
    error('bubble_synth_test:noNight', '找不到夜目录: %s', rawDir);
end

% 先数帧数：注入窗 I0/I1 要按它折算，且 Injector 在投影前就要求 i0/i1 已定
d = dir(fullfile(rawDir, cfg.pattern));
if isempty(d), d = dir(fullfile(rawDir, '*.PNG')); end
if isempty(d), d = dir(fullfile(rawDir, '*.png')); end
nAll = numel(d);
T = numel(1:o.Step:nAll);
if T < 6
    error('bubble_synth_test:tooFew', '夜 %s 只有 %d 帧(Step=%d)，正对照需要更多帧。', ...
        night, T, o.Step);
end

resCell = cell(numel(o.Amp), 1);
for ka = 1:numel(o.Amp)
    amp = o.Amp(ka);
    i0 = o.I0; i1 = o.I1;
    if isempty(i0) || isempty(i1)
        i0 = max(2, round(o.Frac(1) * T));
        i1 = min(T + 1, round(o.Frac(2) * T) + 1);
    end
    if i1 <= i0 + 2
        if o.Verbose
            fprintf('  [跳过] 振幅 %.0f%%: 注入窗太短 (i0=%d i1=%d T=%d)\n', ...
                100 * amp, i0, i1, T);
        end
        continue
    end
    B = struct('amp', amp, 'Akm', o.Akm, 'Bkm', o.Bkm, 'EdgeKm', o.EdgeKm, ...
        'P', o.P, 'DriftMs', o.DriftMs, 'x0', o.X0, 'y0', o.Y0, ...
        'dt', cfg.dtSec, 'i0', i0, 'i1', i1);
    inj = @(i, pr, ctx) inject_one(i, pr, ctx, B);
    tag = sprintf('_syn%02d', round(amp * 100));

    if o.Verbose
        fprintf('\n>>> 注入振幅 %.0f%%\n', 100 * amp);
        fprintf(['    中心 x=%+.2f y=%+.2f deg = %+.0f %+.0f km   ' ...
            '尺寸 %.0f x %.0f km\n'], B.x0, B.y0, B.x0 * KM, B.y0 * KM, ...
            2 * B.Akm, 2 * B.Bkm);
        fprintf('    第 %d ~ %d 帧, 共 %.1f h,  漂移 %.0f m/s 东向\n', ...
            i0, i1, (i1 - i0) * cfg.dtSec / 3600, B.DriftMs);
    end
    out = bubble_night(rawDir, cfg.calFile, o.OutDir, ...
        'DarkFile', cfg.darkFile, 'Zmax', cfg.zMax, 'ResDeg', cfg.resDeg, ...
        'Step', o.Step, 'Injector', inj, 'MakeFig', false, ...
        'MoonCenters', o.MoonCenters, 'MaskDeg', o.MaskDeg, ...
        'Tag', tag, 'Verbose', o.Verbose);

    e = evaluate(out, B, o.Verbose);
    if ~isempty(e)
        if o.MakeFig
            e.fig = figure_one(out, B, e, o.OutDir, night, tag);
        else
            e.fig = '';
        end
        resCell{ka} = e;
    end
    clear out
end

res = [resCell{:}];
if ~isempty(res)
    res = res(~isnan([res.amp]));
end

if o.Verbose && ~isempty(res)
    fprintf('\n%s\n', repmat('=', 1, 96));
    fprintf('正对照汇总（%s）\n', night);
    fprintf('%s\n', repmat('=', 1, 96));
    fprintf('%6s %10s %11s %8s %8s %10s %11s %8s\n', ...
        '注入', '净回收', '掩膜内注入', '回收率', 'σ', '3σ单帧', '3σ时窗', '判决');
    for k = 1:numel(res)
        e = res(k);
        fprintf('%5.0f%% %9.2f%% %10.2f%% %7.0f%% %7.2f%% %9.2f%% %10.2f%% %8s\n', ...
            100 * e.amp, e.recovered, e.expected, 100 * e.ratio, e.sigma, ...
            e.limit, e.limitStack, tern(e.det, '可探测', '不可探测'));
    end
end

% ---------- 落盘 ----------
if ~exist(cfg.logRoot, 'dir'), mkdir(cfg.logRoot); end
save(fullfile(o.OutDir, [night '_synth.mat']), 'res', 'o');
fid = fopen(fullfile(cfg.logRoot, ['正对照_' night '.txt']), 'w', 'n', 'UTF-8');
fprintf(fid, '正对照 %s  (画布 z<=%.0f deg, %.3f deg/px)\n', night, cfg.zMax, cfg.resDeg);
fprintf(fid, '%6s %10s %11s %8s %8s %10s %11s\n', ...
    '注入%', '净回收%', '掩膜内注入%', '回收率%', 'sigma%', '3sig单帧%', '3sig时窗%');
for k = 1:numel(res)
    e = res(k);
    fprintf(fid, '%6.0f %10.2f %11.2f %8.0f %8.2f %10.2f %11.2f\n', ...
        100 * e.amp, e.recovered, e.expected, 100 * e.ratio, e.sigma, e.limit, e.limitStack);
end
fclose(fid);
if o.Verbose
    fprintf('\n已写 -> %s\n', fullfile(o.OutDir, [night '_synth.mat']));
    fprintf('已写 -> %s\n', fullfile(cfg.logRoot, ['正对照_' night '.txt']));
end
end

% ==========================================================================
function d = bfield(i, ctx, B)
% 该帧的耗竭比例 delta(x,y) ∈ [0,1]（0 = 无耗竭）；不在注入窗内返回 []
KM = 111.195;
if i < B.i0 || i >= B.i1
    d = [];
    return
end
t = (i - B.i0) * B.dt;
x0 = B.x0 + B.DriftMs * t / 1000 / KM;          % 东向漂移 -> deg
dx = (ctx.xv(:)' - x0) * KM;                    % 1 x nx  [km]
dy = (ctx.yv(:) - B.y0) * KM;                   % ny x 1  [km]
q = (abs(dx) ./ B.Akm) .^ B.P + (abs(dy) ./ B.Bkm) .^ B.P;
sd = B.Akm * (1 - q) / B.P;                     % 有符号边距的一阶近似 [km]
% sd > 0 在超椭圆**内部** ⇒ 用 +tanh: 内部 -> 1, 外部 -> 0
% （符号写反会变成"椭圆外整体变暗", 回收率恒为 0 —— 这是踩过的坑）
d = B.amp * 0.5 * (1 + tanh(sd / B.EdgeKm));
end

% ==========================================================================
function pr = inject_one(i, pr, ctx, B)
d = bfield(i, ctx, B);
if isempty(d), return; end
pr = pr .* (1 - d);
end

% ==========================================================================
function e = evaluate(out, B, verbose)
dep = out.dep;                      % (ny, nx, T)
valid = logical(out.valid);
xv = out.xv; yv = out.yv;
T = size(dep, 3);
ctx = struct('xv', xv, 'yv', yv);

if B.i1 > T + 1
    B.i1 = T + 1;
end
if B.i1 <= B.i0 + 2
    e = [];
    if verbose, fprintf('  [评估失败] 注入窗太短 (i0=%d i1=%d T=%d)\n', B.i0, B.i1, T); end
    return
end

imid = floor((B.i0 + B.i1) / 2);
dmid = bfield(imid, ctx, B);
coreIm = (dmid > 0.5 * B.amp) & valid;
bh = coreIm(end:-1:1, end:-1:1) & valid;      % 以天顶为原点 => 反极 = 双轴翻转

if sum(coreIm(:)) < 100 || sum(bh(:)) < 100
    e = [];
    if verbose, fprintf('  [评估失败] 核区/反极区太小\n'); end
    return
end

ins = []; expc = []; areas = [];
for i = B.i0:min(B.i1, T)
    d = bfield(i, ctx, B);
    if isempty(d), continue; end
    c = (d > 0.5 * B.amp) & valid;
    if sum(c(:)) < 100, continue; end
    bp = c(end:-1:1, end:-1:1) & valid;
    if sum(bp(:)) < 100, continue; end
    di = double(dep(:, :, i));
    ins(end+1) = mean(di(c)) - mean(di(bp));          %#ok<AGROW>
    expc(end+1) = -100 * mean(d(c));                  %#ok<AGROW>
    areas(end+1) = sum(c(:));                         %#ok<AGROW>
end

ctrl = [1:6:max(B.i0 - 26, 1), (B.i1 + 26):6:T];
ctrl = ctrl(ctrl >= 1 & ctrl <= T);
ons = [];
for q = 1:numel(ctrl)
    di = double(dep(:, :, ctrl(q)));
    ons(end+1) = mean(di(coreIm)) - mean(di(bh));     %#ok<AGROW>
end

if isempty(ins) || isempty(ons)
    e = [];
    if verbose, fprintf('  [评估失败] ins=%d ons=%d\n', numel(ins), numel(ons)); end
    return
end
if numel(ons) < 5 && verbose
    fprintf(['  [注意] 对照帧只有 %d 个（隔帧取样太大时会把对照期挤没）。\n' ...
        '         σ 与 3σ 极限不可信，请用 Step=1 重跑正对照。\n'], numel(ons));
end

rec = mean(ins);
base = median(ons);
expect = mean(expc);
net = rec - base;
if numel(ins) > 1 && numel(ons) > 1
    sd = sqrt((numel(ins) * var(ins) + numel(ons) * var(ons)) / (numel(ins) + numel(ons)));
else
    sd = 0;
end
ratio = net / expect;
lim1 = 3 * sd;
limN = 3 * sd / sqrt(max(numel(ins), 1));
ncore = round(median(areas));
areaKm2 = ncore * (out.km_per_px)^2;

e = struct('amp', B.amp, 'recovered', net, 'expected', expect, 'ratio', ratio, ...
    'sigma', sd, 'limit', lim1, 'limitStack', limN, ...
    'det', abs(net) > lim1, 'ncore', ncore, 'areaKm2', areaKm2, ...
    'nframe', numel(ins));

if verbose
    fprintf('\n%s\n', repmat('=', 1, 92));
    fprintf('正对照结果（注入振幅 %.1f%%，核区 %d 画布像元 = %.0f km2）\n', ...
        100 * B.amp, ncore, areaKm2);
    fprintf('%s\n', repmat('=', 1, 92));
    fprintf('  注入期核区净 D     = %+.2f%%   （掩膜内平均注入 %+.2f%%）\n', rec, expect);
    fprintf('  对照期同位置净 D   = %+.2f%%   (n=%d)\n', base, numel(ons));
    fprintf('  净回收振幅         = %+.2f%%  ->  回收率 %.1f%%\n', net, 100 * ratio);
    fprintf('  逐帧散布 σ         = %.2f%%  （配对核区均值，%d 帧）\n', sd, numel(ins));
    fprintf('  3σ 单帧极限        = %.2f%% 耗竭（面积 %.0f km2）\n', lim1, areaKm2);
    fprintf('  3σ 时窗平均极限    = %.2f%% 耗竭\n', limN);
    fprintf('  判定               = %s\n', tern(abs(net) > lim1, '可探测', '单帧不可探测'));
end
end

% ==========================================================================
function p = figure_one(out, B, e, outDir, night, tag)
sm = out.dep;
T = size(sm, 3);
win = 8;
idx = unique(round(linspace(max(win, 1), max(T - win, 1), 12)));
figure('Position', [40 40 1500 700], 'Visible', 'off');

subplot(2, 3, [1 2 4 5]);
m = median(sm(:, :, :), 3, 'omitnan');
imagesc(out.xv, out.yv, flipud(m), [-25 25]); axis xy; colorbar;
hold on
dmid = bfield(floor((B.i0 + B.i1) / 2), struct('xv', out.xv, 'yv', out.yv), B);
contour(out.xv, out.yv, flipud(double(dmid > 0.5 * B.amp)), [0.5 0.5], 'g', 'LineWidth', 1);
hold off
title(sprintf('%s  整夜中位耗竭图 + 注入轮廓（%.0f%%）', night, 100 * B.amp), ...
    'Interpreter', 'none');

subplot(2, 3, 3);
imagesc(out.xv, out.yv, flipud(dmid) * 100); axis xy; colorbar;
title('注入的真值 δ(x,y)  [%]');
set(gca, 'FontSize', 7);

subplot(2, 3, 6);
i = idx(min(numel(idx), 6));
imagesc(out.xv, out.yv, flipud(sm(:, :, i)), [-25 25]); axis xy; colorbar;
title(sprintf('单帧耗竭图（第 %d 帧）', i));
set(gca, 'FontSize', 7);

sgtitle(sprintf(['%s 正对照  注入 %.0f%%  ->  净回收 %.2f%%  (回收率 %.0f%%)   ' ...
    'σ=%.2f%%  3σ单帧=%.1f%%'], night, 100 * B.amp, e.recovered, ...
    100 * e.ratio, e.sigma, e.limit), 'Interpreter', 'none');
p = fullfile(outDir, [night tag '_control.png']);
exportgraphics(gcf, p); close(gcf);
end

% ==========================================================================
function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
