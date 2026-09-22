function airglow_geo_pipeline(rawFolder, outFolder, calFile, H, obs_lat, obs_lon, varargin)
%AIRGLOW_GEO_PIPELINE 气辉预处理 + 地理投影（精简低内存版）
%
% 只保留核心流程：
%   1) 单帧读取与亮度拉伸；
%   2) 滑动窗口平均背景（仅缓存 2*BgWindow+1 帧）；
%   3) 保存未旋转 *_proc.png；
%   4) 按 cal.rot 建立逆映射，直接从未旋转图采样并保存 *_geo.png；
%   5) 可选：另存带**经纬度坐标轴 + 色标**的 *_geo_ax.png（'Axes', true）。
%
% 删除原版的整晚图像缓存、三维背景栈、keogram、调试 MAT、验证图和多种旋转模式，
% 避免长时间运行时内存持续增长导致 MATLAB 闪退。
%
% 画布几何（'GridMode' / 'SpanDeg'）
%   'km'  ：等地面距离方格（x 按 1/cos(lat0) 换算经度）——正方形千米网，等距。
%   'geo' ：等经纬度方格（1° 经度 = 1° 纬度，不按 cos 换算）——与文献常规图一致。
%           马欣论文表 3.1 的 600x600 @0.02°/px（13.5~25.48°N / 103.22~115.2°E）
%           就是 'geo' 网格：经纬跨度相等（11.98° = 599*0.02）。
%   'SpanDeg' 不给时 = 自动取视场内接正方形；给定则用它（如 12 度）并检查四角是否
%           仍落在定标视场内，超出会 warning。
%   'Zmax'  天顶角限幅(deg)，不给 = 不额外限幅。★ 不加这一项时，"自动"画布**只**被
%           cal.radius（可见盘边）截住，而那可能远在径向标定的实测锚点之外：
%           一期盘边≈87 deg（实测只到 76 deg）⇒ 自动画布会一直铺到 z≈84 deg。
%           传 cfg.zMax 后：
%             · SpanDeg=[]（自动）—— 半径再压到 r(Zmax)，画布真的限在 z≦Zmax；
%             · 显式给了 SpanDeg —— **半径不裁**（既有输出保持不变），只在四角越过
%               r(Zmax) 时给 spanBeyondZmax 告警。
%
% 示例：
%   airglow_geo_pipeline('E:\qihui2\原始数据\20260319', ...
%       'E:\qihui2\Output\out_0319', ...
%       'E:\qihui2\Processed_Output\calibration_params.mat', ...
%       250, 19.5, 109.2, 'GeoFolder', 'E:\qihui2\投影结果\20260319');
%
%   % 与马欣论文一致的画幅（12°见方、0.02°/px、带经纬度标注）
%   airglow_geo_pipeline(raw, proc, calFile, 250, 19.526, 109.133, ...
%       'GridMode', 'geo', 'SpanDeg', 12, 'Axes', true, 'GeoFolder', geo);

    %% ---------- 参数 ----------
    p = inputParser;
    addParameter(p, 'Pattern', 'ODAZH_DCAI01_AROA_L0_STP_*.png');
    addParameter(p, 'BgWindow', 10);
    addParameter(p, 'Threshold', 128);
    addParameter(p, 'TargetMean', 140);
    addParameter(p, 'MaxGain', 4.5);
    addParameter(p, 'Margin', 15);
    addParameter(p, 'ResDeg', 0.02);
    addParameter(p, 'GridMode', 'km');       % 'km' 或 'geo'
    addParameter(p, 'GridStep', 2);
    addParameter(p, 'SpanDeg', []);          % 画布总跨度(deg)；[] = 自动内接正方形
    addParameter(p, 'Zmax', []);             % 天顶角限幅(deg)；[] = 不额外限幅（只由 cal.radius 定）
    addParameter(p, 'Axes', false);          % 是否另存带经纬度坐标轴/色标的 *_geo_ax.png
    addParameter(p, 'SaveProc', true);
    addParameter(p, 'GeoFolder', '');
    parse(p, varargin{:});
    o = p.Results;

    if ~any(strcmpi(o.GridMode, {'km', 'geo'}))
        error('GridMode 只能是 ''km'' 或 ''geo''。');
    end
    if ~isempty(o.SpanDeg) && (~isscalar(o.SpanDeg) || ~isfinite(o.SpanDeg) || o.SpanDeg <= 0)
        error('SpanDeg 必须为正标量（度）。');
    end
    if ~isempty(o.Zmax) && (~isscalar(o.Zmax) || ~isfinite(o.Zmax) || o.Zmax <= 0 || o.Zmax >= 90)
        error('Zmax 必须为 (0, 90) 内的正标量（度）。');
    end
    if ~isscalar(H) || H <= 0 || obs_lat < -90 || obs_lat > 90 || ...
            obs_lon < -180 || obs_lon > 180
        error('H 或测站经纬度无效。');
    end
    o.BgWindow = max(0, round(o.BgWindow));

    %% ---------- 读取标定 ----------
    cal = load_calibration(calFile);
    if isfield(cal, 'obs_lat') && isfield(cal, 'obs_lon') && ...
            isscalar(cal.obs_lat) && isscalar(cal.obs_lon) && ...
            isfinite(cal.obs_lat) && isfinite(cal.obs_lon)
        obs_lat = cal.obs_lat;
        obs_lon = cal.obs_lon;
    end
    cx = cal.zenith(1);
    cy = cal.zenith(2);
    theta_deg = wrap180(cal.rot);      % 当前矩阵约定：正角为逆时针
    rz_eval = make_rz_eval(cal.rz_poly);

    %% ---------- 文件与图像几何 ----------
    fileList = dir(fullfile(rawFolder, o.Pattern));
    if isempty(fileList)
        fileList = dir(fullfile(rawFolder, upper(o.Pattern)));
    end
    nFrames = numel(fileList);
    if nFrames == 0
        error('未找到匹配文件：%s', fullfile(rawFolder, o.Pattern));
    end
    if ~exist(outFolder, 'dir'), mkdir(outFolder); end
    geoFolder = o.GeoFolder;
    if isempty(geoFolder), geoFolder = fullfile(outFolder, 'geo'); end
    if ~exist(geoFolder, 'dir'), mkdir(geoFolder); end

    img0 = imread(fullfile(fileList(1).folder, fileList(1).name));
    [h, w, ~] = size(img0);
    r_boundary = min([cx-1, w-cx, cy-1, h-cy]) - 2;
    r_use = max(1, min([cal.radius - o.Margin, r_boundary]));
    r_from_radius = r_use;                 % 只用盘边算出的半径（仅用于提示/告警）
    % ★ 天顶角限幅（2026-09-17 新增）。不传 Zmax 时行为与以前完全一致。
    %   动机：cal.radius 是**可见盘边**（一期≈87 deg），而径向标定的实测锚点只到
    %   alt 14 deg（z=76 deg）⇒ "自动"画布会越出标定范围、并落进渐晕外环。
    %   cfg.zMax(70) 本来就是"投影限幅天顶角"的意思，这里让它真正生效。
    %   ⚠ 只作用于 SpanDeg=[] 的"自动"分支：显式给了 SpanDeg 的既有输出保持不变，
    %     那种情况只在四角检查里给告警。
    if ~isempty(o.Zmax) && isempty(o.SpanDeg)
        r_zmax = rz_eval(o.Zmax);
        if isfinite(r_zmax) && r_zmax > 0
            r_use = min(r_use, r_zmax);
        end
    end
    [X, Y] = meshgrid(1:w, 1:h);
    mask = hypot(X-cx, Y-cy) <= r_use;
    clear X Y img0;

    theta = deg2rad(theta_deg);
    T1 = [1 0 0; 0 1 0; -cx -cy 1];
    T2 = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
    T3 = [1 0 0; 0 1 0;  cx  cy 1];
    tform_rot = affine2d(T1 * T2 * T3);

    fprintf('>>> %d 帧 | 图像 %dx%d | 天顶(%.2f, %.2f) | 旋转 %.2f° | 半径 %.1f\n', ...
        nFrames, w, h, cx, cy, theta_deg, r_use);
    if ~isempty(o.Zmax) && isempty(o.SpanDeg) && r_use < r_from_radius - 0.5
        fprintf('>>> 已按 Zmax=%.1f° 限幅: 半径 %.1f -> %.1f px（只用盘边 cal.radius=%.1f 会铺到 z≈%.1f°）\n', ...
            o.Zmax, r_from_radius, r_use, cal.radius, invert_rz_poly(rz_eval, r_from_radius));
    end
    fprintf('>>> 测站坐标：lat=%.4f°, lon=%.4f°\n', obs_lat, obs_lon);

    %% ---------- 预计算投影网格 ----------
    Re = 6371;
    R_obs = Re + 0.05;
    R_air = Re + H;
    coslat = cosd(obs_lat);
    if abs(coslat) < 1e-6
        error('测站纬度过于接近极点。');
    end

    d_ground = ground_dist_from_r(r_use, rz_eval, R_obs, R_air);
    if isempty(o.SpanDeg)
        half = d_ground / sqrt(2) * 0.99;   % 内接正方形，避免黑边
    else
        half = o.SpanDeg / 2;
        % 四角是否仍在定标视场内：先算出角点的地心角，与 d_ground 比
        if strcmpi(o.GridMode, 'km')
            dLonHalf = half / coslat;       % km 模式的经度被拉宽 1/cos(lat0)
        else
            dLonHalf = half;
        end
        g_corner = acosd(min(1, ...
            sind(obs_lat) * sind(obs_lat + half) + ...
            cosd(obs_lat) * cosd(obs_lat + half) * cosd(dLonHalf)));
        if g_corner > d_ground
            warning('airglow_geo_pipeline:spanTooBig', ...
                ['SpanDeg=%.2f° 的四角（地心角 %.2f°）超出现有定标视场 %.2f°，' ...
                 '四角会被裁成黑边。最大可用 SpanDeg ≈ %.2f°。'], ...
                o.SpanDeg, g_corner, d_ground, ...
                2 * half * (d_ground / g_corner) * 0.99);
        elseif ~isempty(o.Zmax)
            % 显式 SpanDeg 时刻意不裁半径（保持既有输出），但四角若越过 zMax 的
            % "可用上限"，就得告警：那些像素落在标定锚点之外、且在渐晕外环里。
            r_lim = min(r_from_radius, rz_eval(o.Zmax));
            d_lim = ground_dist_from_r(r_lim, rz_eval, R_obs, R_air);
            if isfinite(d_lim) && g_corner > d_lim
                warning('airglow_geo_pipeline:spanBeyondZmax', ...
                    ['SpanDeg=%.2f° 的四角（地心角 %.2f°）超出 zMax=%.1f° 的可用上限' ...
                     '（对应地心角 %.2f°）；四角像素在标定锚点之外、且落在渐晕外环里，' ...
                     '定量结果不可用（画布本身未裁，仅告警）。'], ...
                    o.SpanDeg, g_corner, o.Zmax, d_lim);
            end
        end
    end

    x_vec = -half:o.ResDeg:half;        % 西 -> 东
    y_vec = -half:o.ResDeg:half;        % 南 -> 北
    [Xg, Yg] = meshgrid(x_vec, y_vec);
    LATg = min(max(obs_lat + Yg, -89.999), 89.999);
    if strcmpi(o.GridMode, 'km')
        LONg = obs_lon + Xg / coslat;
    else
        LONg = obs_lon + Xg;
    end

    sinLat = sind(LATg);
    cosLat = cosd(LATg);
    sinLat0 = sind(obs_lat);
    cosLat0 = cosd(obs_lat);
    dLon = wrap180(LONg - obs_lon);
    cosGamma = min(max(sinLat .* sinLat0 + cosLat .* cosLat0 .* cosd(dLon), -1), 1);
    Gamma = acos(cosGamma);
    zDeg = atan2d(R_air .* sin(Gamma), R_air .* cos(Gamma) - R_obs);
    A_geo = initial_bearing_deg(obs_lat, obs_lon, LATg, LONg);
    A_geo(Gamma < 1e-10) = 0;

    r = rz_eval(zDeg);
    xPix = cx + r .* sind(A_geo);
    yPix = cy - r .* cosd(A_geo);
    valid = isfinite(r) & isfinite(xPix) & isfinite(yPix) & ...
        (zDeg <= 90) & (r >= 0) & (r <= r_use) & ...
        (xPix >= 1) & (xPix <= w) & (yPix >= 1) & (yPix <= h);
    if ~any(valid(:))
        error('投影网格没有有效像素，请检查 rz_poly、radius 和测站坐标。');
    end

    % 坐标轴向量：与 geo/ax_img 的行列一一对应，均单调递增
    lon_axis = obs_lon + x_vec / coslat;
    if strcmpi(o.GridMode, 'geo'), lon_axis = obs_lon + x_vec; end
    lat_axis = obs_lat + y_vec;

    gmap = struct('mode', o.GridMode, 'x0', x_vec(1), 'y_max', y_vec(end), ...
        'res', o.ResDeg, 'coslat', coslat, 'obs_lat', obs_lat, 'obs_lon', obs_lon, ...
        'lon_range', [min(LONg(1,:)) max(LONg(1,:))], ...
        'lat_range', [min(LATg(:,1)) max(LATg(:,1))]);
    fprintf('>>> 投影画布 %dx%d | 有效像素 %.1f%% | lon %.3f～%.3f°, lat %.3f～%.3f°\n', ...
        numel(x_vec), numel(y_vec), 100*mean(valid(:)), ...
        gmap.lon_range(1), gmap.lon_range(2), gmap.lat_range(1), gmap.lat_range(2));
    fprintf('>>> 网格 %s | 跨度 %.2f°E x %.2f°N = %.0f x %.0f km | %.3f°/px = %.2f km/px（南北）| 经纬标注 %s\n', ...
        o.GridMode, gmap.lon_range(2)-gmap.lon_range(1), ...
        gmap.lat_range(2)-gmap.lat_range(1), ...
        111.195 * cosd(obs_lat) * (gmap.lon_range(2)-gmap.lon_range(1)), ...
        111.195 * (gmap.lat_range(2)-gmap.lat_range(1)), ...
        o.ResDeg, o.ResDeg * 111.195, tern_txt(o.Axes, '开', '关'));
    clear Xg Yg LATg LONg sinLat cosLat dLon cosGamma Gamma zDeg A_geo r;

    %% ---------- 流式处理：只保留滑动窗口内的帧 ----------
    W = o.BgWindow;
    bufCap = max(1, 2*W + 1);
    frameBuf = cell(bufCap, 1);
    frameIdx = zeros(bufCap, 1);

    % 预读第一帧所需的前半窗口
    for k = 1:min(nFrames, W+1)
        frameBuf{k} = load_stretched_frame(fullfile(fileList(k).folder, fileList(k).name), mask, o);
        frameIdx(k) = k;
    end

    for i = 1:nFrames
        % 移除窗口左侧旧帧，加入窗口右侧新帧
        old = frameIdx > 0 & frameIdx < i-W;
        frameBuf(old) = {[]};
        frameIdx(old) = 0;
        newIdx = i + W;
        if newIdx <= nFrames && ~any(frameIdx == newIdx)
            slot = find(frameIdx == 0, 1);
            if isempty(slot)
                [~, slot] = min(frameIdx);
            end
            frameBuf{slot} = load_stretched_frame(fullfile(fileList(newIdx).folder, fileList(newIdx).name), mask, o);
            frameIdx(slot) = newIdx;
        end

        pos = find(frameIdx == i, 1);
        if isempty(pos) || isempty(frameBuf{pos})
            continue;
        end
        cur = frameBuf{pos};

        % 平均背景：逐帧累加，不创建 h×w×window 三维数组
        acc = zeros(h, w, 'single');
        cnt = 0;
        for k = 1:bufCap
            if frameIdx(k) > 0 && frameIdx(k) ~= i && ~isempty(frameBuf{k})
                acc = acc + single(frameBuf{k});
                cnt = cnt + 1;
            end
        end
        if cnt > 0
            bg = acc / cnt;
        else
            bg = zeros(h, w, 'single');
        end

        % 背景扣除与 0.5%～99% 对比度归一化
        img_bg = single(cur) - bg;
        v = img_bg(mask);
        v_lo = pct_manual(v, 0.5);
        v_hi = pct_manual(v, 99);
        dyn = v_hi - v_lo;
        proc = zeros(h, w, 'uint8');
        if dyn > 1
            proc(mask) = uint8(round(max(0, min(255, (v - v_lo) * 255 / dyn))));
        end

        [~, name, ~] = fileparts(fileList(i).name);
        if o.SaveProc
            imwrite(proc, fullfile(outFolder, [name, '_proc.png']));
        end

        % 由旋转后坐标逆映射到原始图像，只进行一次双线性插值
        [x_src, y_src] = transformPointsInverse(tform_rot, xPix, yPix);
        src_valid = valid & ...
            x_src >= 1 & x_src <= w & ...
            y_src >= 1 & y_src <= h;

        geo = zeros(size(xPix), 'single');
        geo(src_valid) = interp2(single(proc), ...
            x_src(src_valid), y_src(src_valid), 'linear');
        geo = max(0, min(255, geo));
        geo_img = flipud(uint8(round(geo)));
        RGB = draw_grid_lite(geo_img, gmap, o.GridStep);
        imwrite(RGB, fullfile(geoFolder, [name, '_geo.png']));

        % 带经纬度坐标轴与色标的那一版（论文/汇报直接用）
        if o.Axes
            write_geo_axes(fullfile(geoFolder, [name, '_geo_ax.png']), ...
                uint8(round(geo)), lon_axis, lat_axis, obs_lat, obs_lon, ...
                name, H, o.GridStep);
        end

        if i == 1 || mod(i, 10) == 0 || i == nFrames
            fprintf('  [%d/%d] %s\n', i, nFrames, name);
        end
        clear cur acc bg img_bg v proc x_src y_src src_valid geo geo_img RGB;
    end

    fprintf('完成 → 预处理图：%s | 投影图：%s\n', outFolder, geoFolder);
end

%% ======================== 标定与文件读取 ========================

function cal = load_calibration(calFile)
    if isstruct(calFile)
        cal = calFile;
    else
        S = load(calFile);
        cal = [];
        fn = fieldnames(S);
        for k = 1:numel(fn)
            if isstruct(S.(fn{k})) && isfield(S.(fn{k}), 'zenith')
                cal = S.(fn{k});
                break;
            end
        end
        if isempty(cal) && isfield(S, 'zenith')
            cal = S;
        end
    end
    req = {'zenith', 'rot', 'radius', 'rz_poly'};
    for k = 1:numel(req)
        if ~isfield(cal, req{k})
            error('标定结构体缺少字段：%s', req{k});
        end
    end
    if numel(cal.zenith) < 2 || ~isscalar(cal.rot) || ~isscalar(cal.radius) || cal.radius <= 0
        error('标定参数 zenith/rot/radius 无效。');
    end
end

function rz_eval = make_rz_eval(rz_poly)
    rz = rz_poly(:).';
    if numel(rz) < 2 || ~all(isfinite(rz))
        error('标定参数 rz_poly 无效。');
    end
    if abs(rz(end-1)) > 50
        rz_eval = @(zdeg) polyval(rz, deg2rad(zdeg));   % 弧度系数
    else
        rz_eval = @(zdeg) polyval(rz, zdeg);            % 角度系数
    end
end

function im8 = load_stretched_frame(fpath, mask, o)
    img = q_read_gray(fpath) / 65535;      % 位深对齐(8bit x257) 并归一到 [0,1]
    if isempty(img)
        im8 = [];
        return;
    end
    img = single(img);
    tmp8 = uint8(round(max(0, min(1, img)) * 255));
    mean8 = mean(tmp8(:));
    if mean8 <= o.Threshold
        gain = min(o.TargetMean / max(mean8, 1), o.MaxGain);
        img = min(1, img * gain);
    end
    im8 = uint8(round(max(0, min(1, img)) * 255));
    im8(~mask) = 0;
end

%% ======================== 投影与网格辅助 ========================

function d = wrap180(d)
    d = mod(d + 180, 360) - 180;
end

function dg = ground_dist_from_r(r_px, rz_eval, R_obs, R_air)
%GROUND_DIST_FROM_R  像面半径 r(px) -> 地面地心角距离(deg)。
% 原先是主函数里的一段内联代码；2026-09-17 抽出，好让"zMax 告警"复用同一套物理。
z_m = invert_rz_poly(rz_eval, r_px);
gh  = acos(R_obs / R_air);
if z_m >= 90 - 1e-8
    dg = rad2deg(gh);                      % 越过地平 ⇒ 取地平距离
else
    dg = rad2deg(fzero(@(g) atan2d(R_air .* sin(g), R_air .* cos(g) - R_obs) - z_m, ...
                 [0, gh * (1 - 1e-10)]));
end
end

function z_max = invert_rz_poly(rz_eval, r_target)
    z_grid = linspace(0, 90, 9001);
    r_grid = reshape(rz_eval(z_grid), size(z_grid));
    good = isfinite(r_grid);
    if ~any(good)
        error('rz_poly 在 0–90° 范围内没有有效值。');
    end
    dr = r_grid - r_target;
    idx = find(good(1:end-1) & good(2:end) & dr(1:end-1).*dr(2:end) <= 0, 1, 'first');
    if isempty(idx)
        if all(r_grid(good) < r_target)
            z_max = 90;
            return;
        end
        error('无法由 rz_poly 反解有效半径 %.2f px。', r_target);
    end
    if abs(dr(idx)) < 1e-10
        z_max = z_grid(idx);
    elseif abs(dr(idx+1)) < 1e-10
        z_max = z_grid(idx+1);
    else
        z_max = fzero(@(z) rz_eval(z) - r_target, [z_grid(idx), z_grid(idx+1)]);
    end
    z_max = min(max(z_max, 0), 90);
end

function A = initial_bearing_deg(lat0, lon0, lat, lon)
    dLon = wrap180(lon - lon0);
    A = atan2d(sind(dLon) .* cosd(lat), ...
        cosd(lat0) .* sind(lat) - sind(lat0) .* cosd(lat) .* cosd(dLon));
    A = mod(A, 360);
end

function RGB = draw_grid_lite(geo_img, gmap, step)
    RGB = repmat(geo_img, [1 1 3]);
    [nrows, ncols, ~] = size(RGB);
    if strcmpi(gmap.mode, 'km')
        colfun = @(L) ((L - gmap.obs_lon) * gmap.coslat - gmap.x0) / gmap.res + 1;
    else
        colfun = @(L) ((L - gmap.obs_lon) - gmap.x0) / gmap.res + 1;
    end
    rowfun = @(B) (gmap.y_max - (B - gmap.obs_lat)) / gmap.res + 1;

    lon_lines = ceil(gmap.lon_range(1)/step)*step : step : floor(gmap.lon_range(2)/step)*step;
    lat_lines = ceil(gmap.lat_range(1)/step)*step : step : floor(gmap.lat_range(2)/step)*step;
    for L = lon_lines
        c = round(colfun(L));
        if c >= 1 && c <= ncols
            rows = 1:2:nrows;
            RGB(rows, c, 1) = 255;
            RGB(rows, c, 2) = 220;
            RGB(rows, c, 3) = 0;
        end
    end
    for B = lat_lines
        r = round(rowfun(B));
        if r >= 1 && r <= nrows
            cols = 1:2:ncols;
            RGB(r, cols, 1) = 255;
            RGB(r, cols, 2) = 220;
            RGB(r, cols, 3) = 0;
        end
    end

    % 测站红色十字
    sc = round(colfun(gmap.obs_lon));
    sr = round(rowfun(gmap.obs_lat));
    if sr >= 1 && sr <= nrows && sc >= 1 && sc <= ncols
        arm = 10;
        cs = max(1, sc-arm):min(ncols, sc+arm);
        rs = max(1, sr-arm):min(nrows, sr+arm);
        RGB(sr, cs, 1) = 255; RGB(sr, cs, 2) = 60; RGB(sr, cs, 3) = 60;
        RGB(rs, sc, 1) = 255; RGB(rs, sc, 2) = 60; RGB(rs, sc, 3) = 60;
    end
end

function p = pct_manual(v, pct)
    v = sort(v(:));
    if isempty(v)
        p = 0;
        return;
    end
    k = max(1, min(numel(v), round(pct/100 * numel(v))));
    p = v(k);
end

function s = tern_txt(c, a, b)
    if c, s = a; else, s = b; end
end

%% =====================================================================
%  带经纬度坐标轴 + 色标的投影图
%  img 与 lonv/latv 一一对应：img(i,j) 在 (lonv(j), latv(i))，两者都单调递增。
%  图窗句柄用 persistent 复用，避免逐帧 new/close 造成内存与句柄泄漏。
%% =====================================================================
function write_geo_axes(outFile, img, lonv, latv, lat0, lon0, name, H, stepDeg)
    persistent fg
    if isempty(fg) || ~ishghandle(fg)
        fg = figure('Visible', 'off', 'Color', 'w', 'Position', [60 60 900 720]);
    end
    clf(fg);
    ha = axes('Parent', fg);

    imagesc(ha, lonv, latv, img);
    axis(ha, 'image');
    set(ha, 'YDir', 'normal');
    colormap(ha, gray(256));
    hold(ha, 'on');
    plot(ha, lon0, lat0, 'r+', 'MarkerSize', 16, 'LineWidth', 1.8);
    hold(ha, 'off');

    xt = ceil(lonv(1)/stepDeg)*stepDeg : stepDeg : floor(lonv(end)/stepDeg)*stepDeg;
    yt = ceil(latv(1)/stepDeg)*stepDeg : stepDeg : floor(latv(end)/stepDeg)*stepDeg;
    set(ha, 'XTick', xt, 'YTick', yt, 'TickDir', 'out', 'FontSize', 9, ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'Layer', 'top');
    grid(ha, 'on');
    box(ha, 'on');
    xlabel(ha, '经度 (°E)', 'FontSize', 11);
    ylabel(ha, '纬度 (°N)', 'FontSize', 11);
    title(ha, sprintf('%s UT   |   H = %.0f km   |   红十字符号 = 测站 (%.3f°E, %.3f°N)', ...
        fmt_stamp(name), H, lon0, lat0), 'FontSize', 10);
    cb = colorbar(ha);
    cb.Label.String = '相对亮度（8 bit 拉伸）';
    cb.Label.FontSize = 10;

    try
        exportgraphics(fg, outFile, 'Resolution', 130);
    catch
        print(fg, outFile, '-dpng', '-r130');
    end
end

function s = fmt_stamp(name)
    tok = regexp(name, '(\d{8})(\d{6})', 'tokens', 'once');
    if isempty(tok)
        s = name;
        return
    end
    d = tok{1}; t = tok{2};
    s = sprintf('%s-%s-%s %s:%s:%s', d(1:4), d(5:6), d(7:8), t(1:2), t(3:4), t(5:6));
end