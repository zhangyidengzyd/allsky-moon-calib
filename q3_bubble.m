function out = q3_bubble(nights, varargin)
%Q3_BUBBLE  ③ 泡团一条龙：清晰度排序 → 逐夜处理 → 事件搜索
%
%   q3_bubble('20250920')
%   q3_bubble({'20251119','20251020','20250824'}, 'Step', 2)
%   q3_bubble('all', 'Rank', false)          % rawRoot 下所有夜，不排序
%   q3_bubble('20260326', 'MoonAuto', true)  % 月夜：自动上"月面圆心归一化"
%
% 三步各自的原理
%   ① clarity_rank  —— 先挑干净的夜。云在耗竭率图上给出 10~60% 的斑块，与泡团极像；
%                      实测多数"暗夜"其实整夜有薄云。云指数 = 逐像元时间标准差的中位数。
%   ② bubble_night  —— 核心处理，5 步：减暗场 → 投影到 250 km 画布(一次做掉方位校正)
%                      → 除逐像元时间中位数(除静止结构) → 方位中位数径向归一化
%                      (除随时间变的径向平滑结构) → 减时间中位数 ⇒ D=(I-R)/R (%)
%   ③ bubble_search —— 8% 绝对阈值 + (帧,y,x) 三维连通(≥3 帧) + 二阶矩看南北向拉长
%
% 月夜的处理（本函数自动做）
%   月辉是**以月面为中心**的散射轮廓，对"以天顶为圆心"的归一化来说是方位不对称的大梯度。
%   把径向归一化的圆心从"天顶"搬到"月面"，它就退化成纯径向轮廓被一次除掉；
%   若月盘本身落进视场，再用 'MaskDeg' 把饱和月盘挖掉。
%
% 参数（名值对）
%   'Rank'      先跑清晰度排序，默认 true
%   'Clarity'   外部已算好的清晰度表（struct 数组，含 night/cloud 字段）。
%               给它可以避免重复计算（启动脚本 start_process 就是这么用的）；
%               与 'Rank' 同时给时优先用本参数。
%   'Search'    处理完做事件搜索，默认 true
%   'Step'      隔帧取样，默认 1
%   'Zmax'      投影限幅天顶角(deg)，默认取 cfg.zMax (70)
%   'ResDeg'    画布分辨率，默认取 cfg.resDeg
%   'DarkFile'  暗场路径，默认取 cfg.darkFile；置 '' 则不扣
%   'MoonAuto'  自动判月夜并上"月面圆心归一化 + 月盘掩膜"，默认 true
%   'MoonMaskDeg' 月盘掩膜半径(deg)，默认 0.9
%   'AbsFloor'  搜索的绝对阈值下限(%)，默认 8
%   'Thr'       搜索阈值(%)，默认 [] = 用 AbsFloor
%   'MinFrames' 事件最短持续帧数，默认 3
%   'MinArea'   事件最小面积(画布像元)，默认 1000
%   'MakeFig'   出图，默认 true
%   'Verbose'   默认 true
%
% 输出 out: 结构体数组，每夜一条
%   night / depFile / clarity / moon(是否上月亮圆心) / events / eventsCsv / fig*

p = inputParser;
addParameter(p, 'Rank', true);
addParameter(p, 'Clarity', []);
addParameter(p, 'Search', true);
addParameter(p, 'Step', 1);
addParameter(p, 'Zmax', []);
addParameter(p, 'ResDeg', []);
addParameter(p, 'DarkFile', []);
addParameter(p, 'MoonAuto', true);
addParameter(p, 'MoonMaskDeg', 0.9);
addParameter(p, 'AbsFloor', 8);
addParameter(p, 'Thr', []);
addParameter(p, 'MinFrames', 3);
addParameter(p, 'MinArea', 1000);
addParameter(p, 'MakeFig', true);
addParameter(p, 'Verbose', true);
parse(p, varargin{:});
o = p.Results;

cfg = q_cfg();
if nargin < 1 || isempty(nights)
    error('q3_bubble:noNight', '请给观测夜，例如 q3_bubble(''20250920'')。');
end
if ischar(nights)
    if strcmpi(nights, 'all')
        d = dir(cfg.rawRoot);
        d = d([d.isdir] & ~ismember({d.name}, {'.', '..'}));
        nights = sort({d.name});
    else
        nights = {nights};
    end
end
zMax = o.Zmax; if isempty(zMax), zMax = cfg.zMax; end
resDeg = o.ResDeg; if isempty(resDeg), resDeg = cfg.resDeg; end
darkFile = o.DarkFile; if isempty(darkFile), darkFile = cfg.darkFile; end
if isempty(darkFile) && isempty(o.DarkFile)
    warning('q3_bubble:noDark', '未配置暗场：耗竭率会被 3441 DN 基座稀释一半以上。');
end
if ~exist(cfg.bubRoot, 'dir'), mkdir(cfg.bubRoot); end
if ~exist(cfg.logRoot, 'dir'), mkdir(cfg.logRoot); end

% ---------- ① 清晰度排序 ----------
clar = [];
if o.Rank
    fprintf('\n########## ① 清晰度排序 ##########\n');
    cr = clarity_rank(nights, 'NFrame', 12, 'MakeFig', o.MakeFig, 'Verbose', o.Verbose);
    clar = cr.rows;
end

out = repmat(struct('night', '', 'depFile', '', 'clarity', NaN, ...
    'moon', false, 'events', 0, 'eventsCsv', ''), numel(nights), 1);

for k = 1:numel(nights)
    n = nights{k};
    rawDir = fullfile(cfg.rawRoot, n);
    if ~exist(rawDir, 'dir')
        if o.Verbose, fprintf('\n[%s] 目录不存在，跳过\n', n); end
        continue
    end
    fprintf('\n########## ② 泡团处理  %s ##########\n', n);

    % ---------- 月夜：逐帧月面画布坐标 ----------
    mc = []; maskDeg = 0; isMoon = false;
    if o.MoonAuto
        [mc, maskDeg, isMoon] = moon_centers(rawDir, cfg, o, zMax);
    end
    if o.Verbose && isMoon
        fprintf('月夜判定: 是（月面圆心归一化 + 月盘掩膜 %.2f deg）\n', maskDeg);
    end

    % bubble_night 自己会把 _dep.mat 写到 cfg.bubRoot\bubble<tag>_dep.mat
    tag = ['_' n];
    depFile = fullfile(cfg.bubRoot, ['bubble' tag '_dep.mat']);
    bubble_night(rawDir, cfg.calFile, cfg.bubRoot, ...
        'DarkFile', darkFile, 'Zmax', zMax, 'ResDeg', resDeg, 'Step', o.Step, ...
        'MoonCenters', mc, 'MaskDeg', maskDeg, 'MakeFig', o.MakeFig, ...
        'Tag', tag, 'Verbose', o.Verbose);

    out(k).night = n;
    out(k).depFile = depFile;
    out(k).moon = isMoon;
    % 云指数两个来源：① 自己刚算的 clar；② 调用方（start_process）算好传进来的表。
    % 都不给时保持 NaN —— 汇总表里会打印成 '-' 而不是误导性的 0。
    tbl = clar;
    if isempty(tbl), tbl = o.Clarity; end
    if ~isempty(tbl) && isstruct(tbl) && isfield(tbl, 'night')
        ix = find(strcmp({tbl.night}, n), 1);
        if ~isempty(ix) && isfield(tbl, 'cloud'), out(k).clarity = tbl(ix).cloud; end
    end

    % ---------- ③ 事件搜索 ----------
    if o.Search
        fprintf('\n########## ③ 搜索  %s ##########\n', n);
        r = bubble_search(depFile, cfg.bubRoot, 'Thr', o.Thr, ...
            'AbsFloor', o.AbsFloor, 'MinFrames', o.MinFrames, ...
            'MinArea', o.MinArea, 'MakeFig', o.MakeFig, 'Verbose', o.Verbose);
        out(k).events = numel(r.events);
        if isfield(r, 'eventsCsv'), out(k).eventsCsv = r.eventsCsv; end
    end
end

% ---------- 汇总 ----------
fprintf('\n%s\n', repmat('=', 1, 84));
fprintf('泡团搜索汇总\n');
fprintf('%s\n', repmat('=', 1, 84));
fprintf('%-10s %10s %8s %8s  %s\n', '夜', '云指数%', '月夜', '事件数', '说明');
tot = 0;
for k = 1:numel(out)
    if isempty(out(k).night), continue; end
    tot = tot + out(k).events;
    note = '';
    if out(k).events == 0, note = '无符合判据的候选（云占主导时属正常）'; end
    if isfinite(out(k).clarity), clTxt = sprintf('%10.2f', out(k).clarity);
    else, clTxt = sprintf('%10s', '-'); end
    fprintf('%-10s %s %8s %8d  %s\n', out(k).night, clTxt, ...
        tern(out(k).moon, '是', '否'), out(k).events, note);
end
fprintf('\n合计事件 %d 个。判据: 平滑后 D < -%.0f%%，三维连通 >= %d 帧，面积 >= %d px。\n', ...
    tot, choose(o.Thr, o.AbsFloor), o.MinFrames, o.MinArea);
fprintf('提醒: 泡团应**南北向拉长**且随时间**平移**；无方向性的斑块纹理是云。\n');
end

% ==========================================================================
function [mc, maskDeg, isMoon] = moon_centers(rawDir, cfg, o, zMax)
% 逐帧月面在画布坐标下的位置；顺带判定月盘是否落进视场以决定掩膜半径。
mc = []; maskDeg = 0; isMoon = false;
d = dir(fullfile(rawDir, cfg.pattern));
if isempty(d), d = dir(fullfile(rawDir, '*.PNG')); end
if isempty(d), return; end
files = {d.name};
files = files(1:o.Step:end);
[cen, alt, ~, illum] = moon_canvas_pos(files, cfg.calFile, 'Zmax', zMax, 'Site', cfg.site);
if ~any(alt > 0), return; end
isMoon = true;
mc = cen;
% 月盘是否在视场内: 画布覆盖 z<=zMax，而 z = 90 - alt
if any(alt > (90 - zMax))
    maskDeg = o.MoonMaskDeg;
end
if o.Verbose
    fprintf('  月亮: 高度 %+.1f ~ %+.1f deg, 照亮比 %.2f, 在地平上 %d/%d 帧\n', ...
        min(alt), max(alt), median(illum), sum(alt > 0), numel(alt));
end
end

% ==========================================================================
function c = choose(a, b)
if isempty(a), c = b; else, c = a; end
end

% ==========================================================================
function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
