function out = q5_moon(nights, varargin)
%Q5_MOON  ⑤ 月夜专线：月面圆心归一化 + 月盘掩膜，并标定月夜的实际探测极限
%
%   q5_moon('20260326')
%   q5_moon({'20260323','20260325','20260326','20260327'}, 'Amp', [0.25 0.10])
%   q5_moon('20260326', 'Synth', false)      % 只出耗竭图与搜索，不做正对照
%
% 为什么月夜要特殊处理
%   月辉是**以月面为中心**的散射轮廓。对"以天顶为圆心"的径向归一化来说,
%   它是一个方位不对称的大梯度, 拟合不掉, 会把泡团的对比度一起压掉。
%   把径向归一化的**圆心从"天顶"搬到"月面"**, 月辉就退化成纯径向轮廓,
%   一次除掉; 若月盘本身落进视场, 再用 MaskDeg 把饱和月盘挖掉。
%
% 月夜能看多深（实测，4 夜注入正对照）
%   月辉把耗竭对比度**稀释约 1.7~1.9 倍** ⇒ 月夜只能看较深的泡团（≳20~25%）。
%   注意: 月夜的**第一限制仍然是云, 不是月光**（云指数 22~24% 的那两夜回收率掉到 0）。
%   本函数最后跑的 bubble_synth_test 会给**当前这夜**的真实极限, 不要套用别的夜。
%
% 参数（名值对）
%   'Synth'    是否做注入正对照(定探测极限)，默认 true
%   'Amp'      注入振幅数组，默认 [0.25 0.10]
%   'Search'   是否做事件搜索，默认 true
%   'Step'     隔帧取样，默认 1
%   'Zmax'     投影限幅天顶角(deg)，默认取 cfg.zMax
%   'MakeFig'  出图，默认 true
%   'Verbose'  默认 true
%
% 输出 out: 结构体数组，每夜一条（night / probes / synth / events）

p = inputParser;
addParameter(p, 'Synth', true);
addParameter(p, 'Amp', [0.25 0.10]);
addParameter(p, 'Search', true);
addParameter(p, 'Step', 1);
addParameter(p, 'Zmax', []);
addParameter(p, 'MakeFig', true);
addParameter(p, 'Verbose', true);
parse(p, varargin{:});
o = p.Results;

cfg = q_cfg();
if nargin < 1 || isempty(nights)
    error('q5_moon:noNight', '请给观测夜，例如 q5_moon(''20260326'')。');
end
if ischar(nights), nights = {nights}; end
zMax = o.Zmax; if isempty(zMax), zMax = cfg.zMax; end

fprintf('\n########## ⑤ 月夜专线 ##########\n');
fprintf('做法: 径向归一化圆心 = 月面; 月盘入视场时按 MoonMaskDeg 挖掉。\n');
fprintf('参照: 月辉把耗竭对比度稀释约 1.7~1.9 倍, 月夜只适合找较深的泡团。\n');

% ---- 先做处理 + 搜索（q3_bubble 会自动判月夜并上月亮圆心）----
r3 = q3_bubble(nights, 'Rank', false, 'Search', o.Search, 'Step', o.Step, ...
    'Zmax', zMax, 'MoonAuto', true, 'MakeFig', o.MakeFig, 'Verbose', o.Verbose);

out = repmat(struct('night', '', 'events', 0, 'synth', []), numel(nights), 1);
for k = 1:numel(r3)
    out(k).night = r3(k).night;
    out(k).events = r3(k).events;
end

% ---- 逐夜定探测极限 ----
if o.Synth
    fprintf('\n########## 月夜探测极限（注入正对照）##########\n');
    for k = 1:numel(nights)
        n = nights{k};
        rawDir = fullfile(cfg.rawRoot, n);
        if exist(rawDir, 'dir') ~= 7, continue; end
        fprintf('\n----- %s -----\n', n);
        % 让正对照走**同一套**月亮圆心处理，否则量的是"忘了扣月辉"的极限
        [mc, maskDeg] = moon_centers_full(rawDir, cfg, o.Step, zMax);
        if isempty(mc)
            fprintf('  这一夜月亮不在天上，月夜专线退化为暗夜处理。\n');
        end
        s = bubble_synth_test(n, 'Amp', o.Amp, 'Step', o.Step, ...
            'MoonCenters', mc, 'MaskDeg', maskDeg, 'Zmax', zMax, ...
            'MakeFig', o.MakeFig, 'Verbose', o.Verbose);
        out(k).synth = s;
    end
end

% ---- 月夜判读表 ----
fprintf('\n%s\n', repmat('=', 1, 92));
fprintf('月夜判读\n');
fprintf('%s\n', repmat('=', 1, 92));
fprintf('%-10s %8s %10s | %s\n', '夜', '事件数', '25%%回收率', '判读');
for k = 1:numel(out)
    if isempty(out(k).night), continue; end
    rc = NaN;
    if ~isempty(out(k).synth)
        i25 = find(abs([out(k).synth.amp] - 0.25) < 1e-6, 1);
        if ~isempty(i25), rc = out(k).synth(i25).ratio; end
    end
    if isnan(rc)
        note = '未做正对照';
    elseif rc >= 0.4
        note = '月辉下仍能回收注入 -> 深泡团可判';
    else
        note = '回收率低：该夜云/散射太强，不建议用于泡团判定';
    end
    fprintf('%-10s %8d %9.0f%% | %s\n', out(k).night, out(k).events, 100 * rc, note);
end
fprintf('\n提醒: 判读月夜结果前先看该夜的云指数; 月夜的第一限制是云不是月光。\n');
end

% ==========================================================================
function [mc, maskDeg] = moon_centers_full(rawDir, cfg, step, zMax)
% 与 q3_bubble 内的口径保持一致（这里独立实现，避免跨文件耦合）
mc = []; maskDeg = 0;
d = dir(fullfile(rawDir, cfg.pattern));
if isempty(d), d = dir(fullfile(rawDir, '*.PNG')); end
if isempty(d), return; end
files = {d.name};
files = files(1:step:end);
[cen, alt] = moon_canvas_pos(files, cfg.calFile, 'Zmax', zMax, 'Site', cfg.site);
if ~any(alt > 0), return; end
mc = cen;
if any(alt > (90 - zMax)), maskDeg = 0.9; end
end
