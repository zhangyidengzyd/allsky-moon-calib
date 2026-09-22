function res = clarity_rank(nights, varargin)
%CLARITY_RANK  逐夜"天空清晰度"排序 —— 找泡团之前必须先挑没有云的夜
%
%   res = clarity_rank();                    % 排 rawRoot 下所有夜
%   res = clarity_rank({'20250920','20251119'});
%   res = clarity_rank('all', 'NFrame', 16, 'MakeFig', true);
%
% 为什么需要这一步
%   泡团是气辉里 10~30% 的耗竭。云在半透明时同样给出 10~60% 的亮暗斑块,
%   而且形态更大、更慢, 在耗竭率图上极易与泡团混淆。实测 12 个"暗夜"里有
%   大半整夜被薄云覆盖(原始帧上肉眼可见丝状纹理)。所以先排序, 再在干净的
%   夜里找泡团 —— 这比在任何夜上硬跑有效得多。
%
% 判据（不需要完整处理链, 只要相对可比）
%   ① 逐帧减暗场, 只保留视场内像元;
%   ② 按半径分环取中位当背景并除之 —— 除掉渐晕/径向轮廓等**方位对称**结构;
%      云是方位局地的, 不会被除掉, 所以云留了下来;
%   ③ 对每个像元求整夜**时间标准差** —— 静态固定图案被消掉,
%      剩下的是**随时间变化**的结构, 云是主要来源;
%   ④ cloud = 该场在半径内像元上的**中位数**(%)。数值越小 = 天越干净。
%
% 参数（名值对）
%   'NFrame'   每夜抽多少帧，默认 12
%   'Rmax'     视场内半径(px)，默认 500
%   'NBin'     径向分环数，默认 40
%   'MakeFig'  是否出排序图，默认 true
%   'Verbose'  打印，默认 true
%   'CsvName'  输出 csv 文件名，默认 'clarity.csv'
%
% 输出 res: 结构体数组, 按 cloud 升序 —— res(1) 是最干净的夜

p = inputParser;
addParameter(p, 'NFrame', 12);
addParameter(p, 'Rmax', 500);
addParameter(p, 'NBin', 40);
addParameter(p, 'MakeFig', true);
addParameter(p, 'Verbose', true);
addParameter(p, 'CsvName', 'clarity.csv');
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

% 天顶像点：定标文件里是 1-based，本函数用 0-based
S = load(cfg.calFile);
X0 = double(S.cal.zenith(1)) - 1;
Y0 = double(S.cal.zenith(2)) - 1;

% 暗场
if isempty(cfg.darkFile)
    warning('clarity_rank:noDark', '未给暗场, 云指数会受基座影响(但相对排序仍可用)。');
    dark = 0;
else
    dark = q_read_gray(cfg.darkFile);   % 按位深统一到 16 bit 量级（与帧同一尺度）
end

Rmax = o.Rmax; NBin = o.NBin;
[yy, xx] = ndgrid(0:1023, 0:1023);
rr = hypot(xx - X0, yy - Y0);
m = rr <= Rmax;
rb = min(max(floor(rr / Rmax * NBin) + 1, 1), NBin);

rows = repmat(struct('night', '', 'n', 0, 'cloud', NaN, ...
    'cloud90', NaN, 'level', NaN), numel(nights), 1);
keep = false(numel(nights), 1);
if o.Verbose
    fprintf('\n===== 逐夜清晰度 =====\n');
end
for k = 1:numel(nights)
    n = nights{k};
    dd = fullfile(cfg.rawRoot, n);
    if ~exist(dd, 'dir')
        if o.Verbose, fprintf('  [%s] 目录不存在, 跳过\n', n); end
        continue
    end
    fs = dir(fullfile(dd, cfg.pattern));
    if isempty(fs), fs = dir(fullfile(dd, '*.PNG')); end
    if isempty(fs), fs = dir(fullfile(dd, '*.png')); end
    if numel(fs) < 5
        if o.Verbose, fprintf('  [%s] 帧太少(%d), 跳过\n', n, numel(fs)); end
        continue
    end
    idx = unique(round(linspace(1, numel(fs), min(o.NFrame, numel(fs)))));
    A = zeros(1024, 1024, numel(idx), 'single');
    lev = zeros(numel(idx), 1);
    used = 0;
    for q = 1:numel(idx)
        im = q_read_gray(fullfile(dd, fs(idx(q)).name));   % 按位深统一到 16 bit 量级
        if isempty(im), continue; end
        im = im - dark;
        used = used + 1;
        lev(used) = median(im(m));
        prof = ones(NBin, 1);
        for b = 1:NBin
            mb = m & (rb == b);
            if sum(mb(:)) > 200
                prof(b) = max(median(im(mb)), 1.0);
            end
        end
        A(:, :, used) = single(im ./ prof(rb));
    end
    if used < 5
        if o.Verbose, fprintf('  [%s] 可用帧不足, 跳过\n', n); end
        continue
    end
    A = A(:, :, 1:used);
    tstd = std(A, 0, 3);                       % 逐像元时间标准差
    v = tstd(m);
    rows(k).night = n;
    rows(k).n = used;
    rows(k).cloud = median(v) * 100;
    rows(k).cloud90 = prctile(v, 90) * 100;
    rows(k).level = median(lev(1:used));
    keep(k) = true;
    if o.Verbose
        fprintf('  [%s] 云指数 %6.2f%%  (p90 %6.2f%%)  中位电平 %.0f DN\n', ...
            n, rows(k).cloud, rows(k).cloud90, rows(k).level);
    end
end
rows = rows(keep);
[~, ord] = sort([rows.cloud], 'ascend');
rows = rows(ord);

if o.Verbose
    fprintf('\n%s\n', repmat('=', 1, 62));
    fprintf('清晰度排序（越前越干净, 越适合找泡团）\n');
    fprintf('%s\n', repmat('=', 1, 62));
    for k = 1:numel(rows)
        fprintf('%3d. %s   云指数 %6.2f%%   p90 %6.2f%%   电平 %5.0f DN\n', ...
            k, rows(k).night, rows(k).cloud, rows(k).cloud90, rows(k).level);
    end
end

% ---------- 输出 ----------
if ~exist(cfg.bubRoot, 'dir'), mkdir(cfg.bubRoot); end
csvp = fullfile(cfg.bubRoot, o.CsvName);
fid = fopen(csvp, 'w', 'n', 'UTF-8');
fprintf(fid, 'night,n,cloud,cloud90,level\n');
for k = 1:numel(rows)
    fprintf(fid, '%s,%d,%.3f,%.3f,%.1f\n', rows(k).night, rows(k).n, ...
        rows(k).cloud, rows(k).cloud90, rows(k).level);
end
fclose(fid);
res.rows = rows;
res.csv = csvp;
res.fig = '';
if o.Verbose, fprintf('\n已写 -> %s\n', csvp); end

if o.MakeFig && ~isempty(rows)
    fg = fullfile(cfg.bubRoot, 'figures');
    if ~exist(fg, 'dir'), mkdir(fg); end
    figure('Position', [60 60 900 max(300, 22 * numel(rows) + 140)], 'Visible', 'off');
    barh([rows.cloud]);
    set(gca, 'YDir', 'reverse', 'YTick', 1:numel(rows), ...
        'YTickLabel', {rows.night}, 'FontSize', 8);
    xlabel('云指数 (%)   —— 越小越干净');
    title('逐夜天空清晰度排序（泡团搜索应优先取顶部几夜）');
    grid on; box off;
    res.fig = fullfile(fg, 'clarity_rank.png');
    exportgraphics(gcf, res.fig); close(gcf);
    if o.Verbose, fprintf('已写 -> %s\n', res.fig); end
end
end
