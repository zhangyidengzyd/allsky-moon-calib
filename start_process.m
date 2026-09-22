function out = start_process(nights, varargin)
%START_PROCESS  ★启动脚本 ②  图像处理 + 地理投影★
%
%   start_process('20250920')
%   start_process({'20250920','20251119'}, 'Step', 2)
%   start_process('all', 'TopN', 5)              % 全部夜里挑云最少的前 5 夜
%   start_process('20250920', 'Project', false)  % 只做泡团定量，不出投影图
%   start_process('20250920', 'Bubble', false)   % 只出投影图，不做泡团定量
%
% ■ 两步各自是什么
%   ① q2_project  图像处理 + 地理投影（8 bit 出图，**为看形态**）
%      读帧 → 亮度拉伸 → 邻帧滑动平均背景 → 按 cal.rot 建逆映射直接采到地理网格。
%      ★ 方位校正在这一步就做掉了，不需要再单独旋转一次。
%      输出 <cfg.projRoot>\<夜>\*_geo.png
%   ② q3_bubble   泡团定量（float 域，**出数据**）
%      减暗场 → 投影到 250 km 画布 → 除逐像元时间中位数 → 方位中位数径向归一化
%      → 减时间中位数 ⇒ D = (I-R)/R (%) → 8% 绝对阈值 + (帧,y,x) 三维连通找事件。
%      输出 <cfg.bubRoot>\bubble_<夜>_dep.mat 与 events.csv
%
% ⚠ 8 bit 会把百分之几的耗竭压掉，所以 ① 只能看形态；要振幅一律看 ②。
% ⚠ 月夜由 q3_bubble 自动识别并把径向归一化圆心从天顶搬到月面，无需手动干预。
%   若月盘落进视场，再按 MoonMaskDeg 把饱和月盘挖掉。
%
% ■ 换原始图像文件夹
%   只改 q_cfg.m 第 21 行 cfg.rawRoot。本脚本从 q_cfg 取全部路径，别处不用动。
%
% ■ 参数（名值对）
%   'Step'      隔帧取样，默认 1
%   'Project'   是否做 ① 投影出图，默认 true
%   'Bubble'    是否做 ② 泡团定量，默认 true
%   'Search'    ② 里是否做事件搜索，默认 true
%   'Rank'      是否先按云指数排序，默认 true
%   'TopN'      >0 时先排序、只跑云最少的前 N 夜，默认 0 = 不限
%   'ResDeg'    画布分辨率(deg)，默认取 cfg.resDeg (0.02)
%   'Zmax'      投影限幅天顶角(deg)，默认取 cfg.zMax
%   'GridMode'  ① 投影的网格模式：'geo'（默认，等经纬度）或 'km'（等地面距离）
%   'SpanDeg'   ① 画布总跨度(deg)，默认 11.98 = 马欣论文表 3.1 的 600x600 @0.02°/px；[] = 自动
%   'Axes'      ① 是否另存带经纬度坐标轴 + 色标的 <名>_geo_ax.png，默认 true
%   'GridStep'  经纬网线间隔(deg)，默认 2
%   'BgWindow'  邻帧背景半窗，默认 10；给 0 表示不做背景
%   'MakeFig'   出图，默认 true（只作用于 ② 与清晰度排序；① 的投影图始终会出）
%   'Diary'     全过程写日志文件，默认 true
%   'Verbose'   默认 true
%
% 输出 out: nights / proj / bubble / table / logFile / elapsedSec

p = inputParser;
addParameter(p, 'Step', 1);
addParameter(p, 'Project', true);
addParameter(p, 'Bubble', true);
addParameter(p, 'Search', true);
addParameter(p, 'Rank', true);
addParameter(p, 'TopN', 0);
addParameter(p, 'ResDeg', []);
addParameter(p, 'Zmax', []);
addParameter(p, 'GridMode', 'geo');
addParameter(p, 'SpanDeg', 11.98);
addParameter(p, 'Axes', true);
addParameter(p, 'GridStep', 2);
addParameter(p, 'BgWindow', 10);
addParameter(p, 'MakeFig', true);
addParameter(p, 'Diary', true);
addParameter(p, 'Verbose', true);
parse(p, varargin{:});
o = p.Results;

t0 = tic;

% ==========================================================================
% 0  定位工程根（脚本在哪儿，工程根就在哪儿 —— 不写死盘符）
% ==========================================================================
here = fileparts(mfilename('fullpath'));
if ~isempty(here)
    cd(here);
    addpath(here);
end
if exist('q_cfg.m', 'file') ~= 2
    error('start_process:noCfg', '找不到 q_cfg.m —— 本脚本必须与它放在同一目录。');
end

% ==========================================================================
% 0.1  日志（diary 落盘；函数退出/报错时自动关闭）
% ==========================================================================
cfg = q_cfg();
logFile = '';
if o.Diary
    if exist(cfg.logRoot, 'dir') ~= 7, mkdir(cfg.logRoot); end
    logFile = fullfile(cfg.logRoot, ...
        ['start_process_' datestr(now, 'yyyymmdd_HHMMSS') '.txt']);
    diary(logFile); diary on;
    % 这个变量必须存活到函数结束（含中途报错），diary 才会被关掉
    cleanupDiary = onCleanup(@() diary('off'));
end

% ==========================================================================
fprintf('\n%s\n', repmat('#', 1, 78));
fprintf('#  启动脚本 ②  图像处理 + 地理投影   start_process\n');
fprintf('#  时间 %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('#  工程根 %s\n', here);
fprintf('%s\n', repmat('#', 1, 78));

fprintf('\n########## 0  路径体检 ##########\n');
q_cfg('check');

% 定标是投影的前提，缺了直接给结论，不浪费一轮读图
fprintf('\n########## 0.1  定标依赖 ##########\n');
if exist(cfg.calFile, 'file') ~= 2
    error('start_process:noCal', ...
        ['没有管线定标文件，投影无法进行：\n      %s\n' ...
         '  请先跑 start_calib（或 q1_calibrate）生成它。'], cfg.calFile);
end
S = load(cfg.calFile);
if isfield(S, 'cal')
    cal = S.cal;
    fprintf('  定标 : %s\n', cfg.calFile);
    fprintf('         zenith(1-based) = (%.3f, %.3f)   rot = %+.4f deg\n', ...
        double(cal.zenith(1)), double(cal.zenith(2)), double(cal.rot));
    fprintf('         视场限幅 z <= %.0f deg（km/px：z=30→0.64，z=70→1.34，z=88→3.12）\n', ...
        cfg.zMax);
else
    error('start_process:badCal', '定标文件缺少 cal 变量：%s', cfg.calFile);
end

% ==========================================================================
% 1  解析 "处理哪几夜"
% ==========================================================================
if nargin < 1 || isempty(nights), nights = 'all'; end
if ischar(nights) && strcmpi(nights, 'all'), nights = {}; end
if ischar(nights), nights = {nights}; end
if isempty(nights)
    d = dir(cfg.rawRoot);
    d = d([d.isdir] & ~ismember({d.name}, {'.', '..'}));
    nights = sort({d.name});
end
if isempty(nights)
    error('start_process:noNight', 'rawRoot 下没有日期目录：%s', cfg.rawRoot);
end
fprintf('\n########## 1  处理清单 ##########\n');
fprintf('  共 %d 夜：%s\n', numel(nights), strjoin(nights, ', '));

% ==========================================================================
% 2  清晰度排序（只算一次；TopN 也用它来挑夜）
% ==========================================================================
clar = [];
if o.Bubble && o.Rank
    fprintf('\n########## 2  清晰度排序（云指数越小越干净）##########\n');
    cr = clarity_rank(nights, 'NFrame', 12, 'MakeFig', o.MakeFig, ...
        'Verbose', o.Verbose);
    clar = cr.rows;
end

if o.TopN > 0 && numel(nights) > 1
    if isempty(clar)
        cr = clarity_rank(nights, 'NFrame', 12, 'MakeFig', o.MakeFig, ...
            'Verbose', o.Verbose);
        clar = cr.rows;
    end
    m = min(o.TopN, numel(clar));
    nights = {clar(1:m).night};
    fprintf('\n########## 2.1  按 TopN 取云最少的前 %d 夜 ##########\n', m);
    for k = 1:m
        fprintf('  %2d. %-10s 云指数 %6.2f%%\n', k, clar(k).night, clar(k).cloud);
    end
end

% ==========================================================================
% 3  逐夜：① 投影出图  →  ② 泡团定量
% ==========================================================================
projOut = repmat(struct('night', '', 'geoDir', '', 'procDir', '', 'nFrames', 0), ...
    numel(nights), 1);
bubOut = repmat(struct('night', '', 'depFile', '', 'clarity', NaN, ...
    'moon', false, 'events', 0, 'eventsCsv', ''), numel(nights), 1);

for k = 1:numel(nights)
    n = nights{k};
    fprintf('\n%s\n', repmat('#', 1, 78));
    fprintf('#  [%d/%d]  %s\n', k, numel(nights), n);
    fprintf('%s\n', repmat('#', 1, 78));

    rawDir = fullfile(cfg.rawRoot, n);
    if exist(rawDir, 'dir') ~= 7
        fprintf('[跳过] 目录不存在：%s\n', rawDir);
        continue
    end

    % ---------- ① 图像处理 + 地理投影（8 bit，看形态） ----------
    if o.Project
        try
            % q2_project 的参数集合与 q3_bubble 不同：它没有 MakeFig
            r2 = q2_project(n, 'ResDeg', o.ResDeg, 'GridStep', o.GridStep, ...
                'BgWindow', o.BgWindow, 'GridMode', o.GridMode, ...
                'SpanDeg', o.SpanDeg, 'Axes', o.Axes, 'Verbose', o.Verbose);
            projOut(k).night = n;
            projOut(k).geoDir = r2(1).geo;
            projOut(k).procDir = r2(1).proc;
            projOut(k).nFrames = r2(1).nFrames;
        catch ME
            fprintf('[投影失败] %s : %s\n', n, ME.message);
        end
    end

    % ---------- ② 泡团定量（float 域，出数据） ----------
    if o.Bubble
        try
            % 清晰度已在第 2 步算过，这里用 'Clarity' 传进去，避免重复计算
            r3 = q3_bubble(n, 'Rank', false, 'Clarity', clar, 'Search', o.Search, ...
                'Step', o.Step, 'Zmax', o.Zmax, 'ResDeg', o.ResDeg, ...
                'MoonAuto', true, 'MakeFig', o.MakeFig, 'Verbose', o.Verbose);
            bubOut(k).night = n;
            bubOut(k).depFile = r3(1).depFile;
            bubOut(k).clarity = r3(1).clarity;
            bubOut(k).moon = r3(1).moon;
            bubOut(k).events = r3(1).events;
            if isfield(r3(1), 'eventsCsv')
                bubOut(k).eventsCsv = r3(1).eventsCsv;
            end
        catch ME
            fprintf('[泡团处理失败] %s : %s\n', n, ME.message);
        end
    end
end

% ==========================================================================
% 4  汇总
% ==========================================================================
fprintf('\n%s\n', repmat('=', 1, 100));
fprintf('汇总（%d 夜）\n', numel(nights));
fprintf('%s\n', repmat('=', 1, 100));
fprintf('%-10s %8s %6s %9s %8s  %s\n', '夜', '云指数%', '月夜', '投影图数', '事件数', '输出');
totEv = 0;
totFig = 0;
for k = 1:numel(nights)
    if isempty(projOut(k).night) && isempty(bubOut(k).night), continue; end
    totEv = totEv + bubOut(k).events;
    totFig = totFig + projOut(k).nFrames;
    cl = bubOut(k).clarity;
    if ~isfinite(cl) && ~isempty(clar)
        ix = find(strcmp({clar.night}, nights{k}), 1);
        if ~isempty(ix), cl = clar(ix).cloud; end
    end
    if isfinite(cl), clTxt = sprintf('%8.2f', cl); else, clTxt = '       -'; end
    fprintf('%-10s %s %6s %9d %8d  %s\n', nights{k}, clTxt, ...
        tern(bubOut(k).moon, '是', '否'), projOut(k).nFrames, bubOut(k).events, ...
        choose(projOut(k).geoDir, cfg.bubRoot));
end
fprintf('\n投影图合计 %d 张，事件合计 %d 个。\n', totFig, totEv);
fprintf('判据：平滑后 D < -8%%，三维连通 >= 3 帧，面积 >= 1000 px。\n');
fprintf('提醒：泡团应**南北向拉长**且随时间**平移**；无方向性的斑块纹理是云。\n');

% ==========================================================================
out = struct();
out.nights = nights;
out.proj = projOut;
out.bubble = bubOut;
out.clarity = clar;
out.logFile = logFile;
out.elapsedSec = toc(t0);

fprintf('\n%s\n', repmat('=', 1, 78));
fprintf('图像处理 + 地理投影结束，用时 %.1f s\n', out.elapsedSec);
fprintf('  投影图（看形态，8 bit） : %s\\<夜>\\*_geo.png\n', cfg.projRoot);
fprintf('  耗竭率（定量，float）   : %s\\bubble_<夜>_dep.mat\n', cfg.bubRoot);
fprintf('  事件清单                : %s\\*_events.csv\n', cfg.bubRoot);
if ~isempty(logFile)
    fprintf('  日志                    : %s\n', logFile);
end
fprintf('下一步：q5_moon(''20260326'') 处理月夜；q4_verify 随时自检。\n');
fprintf('%s\n', repmat('=', 1, 78));
end

% ==========================================================================
function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

% ==========================================================================
function c = choose(a, b)
if isempty(a), c = b; else, c = a; end
end
