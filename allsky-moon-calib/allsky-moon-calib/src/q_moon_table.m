function T = q_moon_table(nights, varargin)
%Q_MOON_TABLE  逐夜日月位置表 —— 判断哪些夜是"真暗夜"、哪些是月夜
%
%   q_moon_table                 % rawRoot 下所有夜
%   q_moon_table({'20260326'})
%   q_moon_table('all', 'Step', 10)
%
% 做什么
%   对每一夜：从**帧文件名**里解析时间戳（本链所有帧名都含 yyyymmddHHMMSS），
%   取起/中/止三个时刻，用 moon_altaz_ref.m（Meeus 47 章，含站心视差）算月亮的
%   地平高度、方位、照亮比。据此判该夜是 月夜 / 暗夜。
%
% 为什么不用现成的 txt
%   旧版本读的是 Python 从帧内 tEXt 导出的 times_all.txt。本函数**直接读文件名的
%   时间戳**，因为已实测确认：文件名与帧内 CreationTime 逐帧一致且唯一，
%   整夜跨度 12:02–21:11 ⇒ 是 **UT**（对应当地 20:02–次日 05:11，正好一夜）。
%   这样纯 MATLAB 链不再依赖任何 Python 中间产物。
%
% 参数（名值对）
%   'Step'     抽样步长（只影响统计"月亮在地平上的帧数"），默认 5
%   'Write'    是否写文本文件到 cfg.logRoot，默认 true
%   'Verbose'  打印，默认 true
%
% 输出 T: 结构体数组，每夜一条
%   night / t0 / tm / t1 / alt(3x1) / az(3x1) / illum / moonUp(比例) / isMoonNight

p = inputParser;
addParameter(p, 'Step', 5);
addParameter(p, 'Write', true);
addParameter(p, 'Verbose', true);
parse(p, varargin{:});
o = p.Results;

cfg = q_cfg();
if nargin < 1 || isempty(nights) || (ischar(nights) && strcmpi(nights, 'all'))
    d = dir(cfg.rawRoot);
    d = d([d.isdir] & ~ismember({d.name}, {'.', '..'}));
    nights = sort({d.name});
elseif ischar(nights)
    nights = {nights};
end

lon = cfg.site(1); lat = cfg.site(2); hgt = cfg.site(3);
T = repmat(struct('night', '', 't0', '', 'tm', '', 't1', '', ...
    'alt', zeros(3, 1), 'az', zeros(3, 1), 'illum', NaN, ...
    'moonUp', NaN, 'isMoonNight', false, 'n', 0), numel(nights), 1);
nkeep = 0;

for k = 1:numel(nights)
    n = nights{k};
    dd = fullfile(cfg.rawRoot, n);
    if exist(dd, 'dir') ~= 7, continue; end
    fl = dir(fullfile(dd, cfg.pattern));
    if isempty(fl), fl = dir(fullfile(dd, '*.PNG')); end
    if isempty(fl), fl = dir(fullfile(dd, '*.png')); end
    if isempty(fl), continue; end
    names = {fl.name};

    tt = NaT(numel(names), 1, 'TimeZone', 'UTC');
    ok = false(numel(names), 1);
    for q = 1:numel(names)
        tok = regexp(names{q}, '(\d{8})(\d{6})', 'tokens');
        if isempty(tok), continue; end
        s = tok{1}{1};
        tt(q) = datetime([s(1:4) '-' s(5:6) '-' s(7:8) ' ' ...
            tok{1}{2}(1:2) ':' tok{1}{2}(3:4) ':' tok{1}{2}(5:6)], 'TimeZone', 'UTC');
        ok(q) = true;
    end
    tt = tt(ok);
    if isempty(tt), continue; end
    tt = sort(tt);
    nAll = numel(ok);

    idx = unique(round(linspace(1, numel(tt), 3)));
    ts = tt(idx);
    [alt, az, ~, ill] = moon_altaz_ref(ts, lon, lat, hgt);

    sub = tt(1:o.Step:end);
    [altAll, ~, ~, illAll] = moon_altaz_ref(sub, lon, lat, hgt); %#ok<ASGLU>
    frac = mean(altAll > 0);

    nkeep = nkeep + 1;
    T(nkeep).night = n;
    T(nkeep).t0 = datestr(ts(1), 'yyyy-mm-dd HH:MM');
    T(nkeep).tm = datestr(ts(2), 'yyyy-mm-dd HH:MM');
    T(nkeep).t1 = datestr(ts(3), 'yyyy-mm-dd HH:MM');
    T(nkeep).alt = alt(:);
    T(nkeep).az = az(:);
    T(nkeep).illum = ill(2);
    T(nkeep).moonUp = frac;
    T(nkeep).isMoonNight = frac > 0.05;
    T(nkeep).n = nAll;
end
T = T(1:nkeep);

if o.Verbose
    fprintf('\n%s\n', repmat('=', 1, 100));
    fprintf('逐夜日月位置（站点 %.3fE %.3fN %.0fm，时间戳按 UT）\n', lon, lat, hgt * 1000);
    fprintf('%s\n', repmat('=', 1, 100));
    fprintf('%-10s %6s %16s %20s %20s %16s %7s %8s %s\n', ...
        '夜', '帧数', '起始(UT)', '月高 起/中/止 (deg)', '月方 起/中/止 (deg)', ...
        '照亮比', '月在地平', '类别', '');
    for k = 1:numel(T)
        fprintf('%-10s %6d %16s %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %16.3f %6.0f%% %8s\n', ...
            T(k).night, T(k).n, T(k).t0, T(k).alt(1), T(k).alt(2), T(k).alt(3), ...
            T(k).az(1), T(k).az(2), T(k).az(3), T(k).illum, ...
            100 * T(k).moonUp, tern(T(k).isMoonNight, '月夜', '暗夜'));
    end
    nm = sum([T.isMoonNight]);
    fprintf('\n共 %d 夜：月夜 %d，暗夜 %d。\n', numel(T), nm, numel(T) - nm);
    if nm > 0
        fprintf('月夜清单: %s\n', strjoin({T([T.isMoonNight]).night}, ', '));
    end
end

if o.Write
    if ~exist(cfg.logRoot, 'dir'), mkdir(cfg.logRoot); end
    fp = fullfile(cfg.logRoot, '逐夜日月位置.txt');
    fid = fopen(fp, 'w', 'n', 'UTF-8');
    fprintf(fid, '夜,帧数,起始UT,月高起,月高中,月高止,月方起,月方中,月方止,照亮比,月在地平,类别\n');
    for k = 1:numel(T)
        fprintf(fid, '%s,%d,%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.3f,%.3f,%s\n', ...
            T(k).night, T(k).n, T(k).t0, T(k).alt(1), T(k).alt(2), T(k).alt(3), ...
            T(k).az(1), T(k).az(2), T(k).az(3), T(k).illum, T(k).moonUp, ...
            tern(T(k).isMoonNight, '月夜', '暗夜'));
    end
    fclose(fid);
    if o.Verbose, fprintf('已写 -> %s\n', fp); end
end
end

% ==========================================================================
function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
