function out = q2_project(nights, varargin)
%Q2_PROJECT  ② 地理投影一条龙：原图 → 方位校正 → 250 km 地理画布
%
%   q2_project('20250920')
%   q2_project({'20250920','20251119'}, 'ResDeg', 0.02, 'GridStep', 2)
%   q2_project('20250920', 'BgWindow', 0)      % 不做邻帧背景，看原始投影
%
% 原理
%   ① 单帧读取 + 亮度拉伸（8 bit，为**人眼看形态**服务）；
%   ② 滑动窗口平均背景（只缓存 2*BgWindow+1 帧，避免整夜缓存把 MATLAB 撑爆）；
%   ③ 存未旋转的 *_proc.png；
%   ④ 按 cal.rot 建**逆映射**，直接从未旋转图采样到地理网格 → *_geo.png。
%      关键：这一步**同时把方位校正做掉了** —— 不需要再单独旋转一次。
%   画布几何：以测站为原点的等距经纬网格；天顶角 z 由视线与 250 km 球壳求交得到。
%
% ⚠ 输出是 8 bit 拉伸图，只能看形态。要**定量的耗竭率**请用 q3_bubble（bubble_night.m），
%   它保留 float 域；8 bit 会把百分之几的耗竭信息压掉。
%
% 参数（名值对）
%   'ResDeg'    画布分辨率(deg)，默认取 cfg.resDeg (0.02)
%   'GridMode'  'geo'（默认，等经纬度）或 'km'（等地面距离）
%   'SpanDeg'   画布总跨度(deg)，默认 11.98 = 599x0.02，即与马欣论文表 3.1 完全相同的
%               600x600 @0.02°/px；给 12 得 601x601；给 [] 表示自动取视场内接正方形
%   'Axes'      是否另存带经纬度坐标轴 + 色标的 <名>_geo_ax.png，默认 true
%   'GridStep'  经纬网线间隔(deg)，默认 2
%   'BgWindow'  邻帧背景半窗，默认 10；给 0 表示不做背景
%   'Margin'    视场边缘留白(px)，默认 15
%   'Threshold'/'TargetMean'/'MaxGain'  拉伸参数，默认 128 / 140 / 4.5
%   'SaveProc'  是否存未旋转的 proc 图，默认 true
%   'Verbose'   默认 true
%
% 关于画幅（为什么默认 11.98° / 'geo'）
%   马欣论文表 3.1：投影图 600x600、0.02°/px，覆盖 13.5~25.48°N / 103.22~115.2°E
%   —— 经纬跨度都是 11.98°（= 599 x 0.02），说明那是一张**等经纬度**网格（'geo'）。
%   本函数默认 'geo' + SpanDeg 11.98，画布 600x600 @0.02°，以本测站 (109.133E, 19.526N)
%   为中心 -> 13.536~25.516°N / 103.143~115.123°E。与论文残余的 ≤0.09° 差异全部来自
%   测站坐标不同（论文站 ≈109.2E, 19.49N），不是投影误差。
%   ⚠ 高度不同：论文用 H = 300 km，本链用 cfg.hShell = 250 km。画布**范围**与 H 无关，
%     但同一纹理落在画布里的**位置**会随 H 变。要与论文逐点比对必须统一 H。
%
% 输出 out: 结构体数组，每夜一条（night / raw / proc / geo / nFrames / geoLon / geoLat / spanKm）

p = inputParser;
addParameter(p, 'ResDeg', []);
addParameter(p, 'GridMode', 'geo');
addParameter(p, 'SpanDeg', 11.98);
addParameter(p, 'Axes', true);
addParameter(p, 'GridStep', 2);
addParameter(p, 'BgWindow', 10);
addParameter(p, 'Margin', 15);
addParameter(p, 'Threshold', 128);
addParameter(p, 'TargetMean', 140);
addParameter(p, 'MaxGain', 4.5);
addParameter(p, 'SaveProc', true);
addParameter(p, 'Verbose', true);
parse(p, varargin{:});
o = p.Results;

cfg = q_cfg();
if nargin < 1 || isempty(nights)
    error('q2_project:noNight', '请给观测夜，例如 q2_project(''20250920'')。');
end
if ischar(nights), nights = {nights}; end
resDeg = o.ResDeg;
if isempty(resDeg), resDeg = cfg.resDeg; end
if ~exist(cfg.projRoot, 'dir'), mkdir(cfg.projRoot); end

% 管线会用 cal.obs_lat/obs_lon 覆盖入参 —— 这里读出来，好让输出的经纬范围与图一致
calAux = struct();
try
    Sc = load(cfg.calFile);
    if isfield(Sc, 'cal'), calAux = Sc.cal; end
catch
    calAux = struct();
end

out = repmat(struct('night', '', 'raw', '', 'proc', '', 'geo', '', ...
    'nFrames', 0, 'geoLon', [], 'geoLat', [], 'spanKm', []), numel(nights), 1);

for k = 1:numel(nights)
    n = nights{k};
    rawDir = fullfile(cfg.rawRoot, n);
    if ~exist(rawDir, 'dir')
        error('q2_project:noDir', '找不到夜目录: %s', rawDir);
    end
    procDir = fullfile(cfg.projRoot, [n '_proc']);
    geoDir = fullfile(cfg.projRoot, n);

    if o.Verbose
        fprintf('\n########## ② 地理投影  %s ##########\n', n);
        fprintf('原始帧   : %s\n', rawDir);
        fprintf('中间产物 : %s\n', procDir);
        fprintf('投影输出 : %s\n', geoDir);
        fprintf('画幅     : %s 网格, 跨度 %s deg, %.3f deg/px\n', ...
            o.GridMode, tern_str(isempty(o.SpanDeg), '自动', sprintf('%.2f', o.SpanDeg)), resDeg);
    end

    airglow_geo_pipeline(rawDir, procDir, cfg.calFile, cfg.hShell, ...
        cfg.site(2), cfg.site(1), ...
        'Pattern', cfg.pattern, 'ResDeg', resDeg, 'GridStep', o.GridStep, ...
        'GridMode', o.GridMode, 'SpanDeg', o.SpanDeg, 'Zmax', cfg.zMax, ...
        'Axes', o.Axes, ...
        'BgWindow', o.BgWindow, 'Margin', o.Margin, ...
        'Threshold', o.Threshold, 'TargetMean', o.TargetMean, ...
        'MaxGain', o.MaxGain, 'SaveProc', o.SaveProc, 'GeoFolder', geoDir);

    d = dir(fullfile(geoDir, '*_geo.png'));
    out(k).night = n;
    out(k).raw = rawDir;
    out(k).proc = procDir;
    out(k).geo = geoDir;
    out(k).nFrames = numel(d);
    % 画布经纬范围（与管线同约定：节点式 -half:res:half；测站坐标以 cal 为准）
    if ~isempty(o.SpanDeg)
        if isfield(calAux, 'obs_lat'), la0 = calAux.obs_lat; lo0 = calAux.obs_lon;
        else, la0 = cfg.site(2); lo0 = cfg.site(1); end
        hv = -o.SpanDeg/2 : resDeg : o.SpanDeg/2;
        if strcmpi(o.GridMode, 'km'), lov = lo0 + hv / cosd(la0);
        else, lov = lo0 + hv; end
        out(k).geoLon = [lov(1) lov(end)];
        out(k).geoLat = [la0 + hv(1) la0 + hv(end)];
        out(k).spanKm = [111.195 * cosd(la0) * (lov(end)-lov(1)), ...
                         111.195 * (hv(end)-hv(1))];
        if o.Verbose
            fprintf('画布经纬 : lon %.3f~%.3f °E, lat %.3f~%.3f °N  (东西 %.0f km x 南北 %.0f km)\n', ...
                out(k).geoLon(1), out(k).geoLon(2), out(k).geoLat(1), out(k).geoLat(2), ...
                out(k).spanKm(1), out(k).spanKm(2));
        end
    end
    if o.Verbose
        fprintf('完成: %d 张投影图 -> %s\n', numel(d), geoDir);
    end
end
end

% ==========================================================================
function s = tern_str(c, a, b)
if c, s = a; else, s = b; end
end
