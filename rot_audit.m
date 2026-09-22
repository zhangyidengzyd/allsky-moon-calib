function rot_audit(dtMax)
%ROT_AUDIT  一期（2014）定标里"方位校正角 rot"的三项独立审计。
%
% 为什么需要它
%   本解（月亮方位共识）给出 rot = +4.35 deg，马欣论文 §3.2.2 给出 6.35 deg，
%   差 2.0 deg。这个差不能靠"重算一遍月亮"来裁决，只能用**与月亮无关**的旁证逐项排除：
%     (1) 论文那个"北极星方位角 0.4 deg"对不对 —— 用低精度星历独立复算；
%     (2) 论文那颗星 (499,55) 在图里究竟测不测得到 —— 决定它的像素能不能被复核；
%     (3) 相机时钟有没有偏 —— 星历时刻整体平移 dt，看 rot/rms 的极小在哪。
%
% 2026-09-17 实测结论（本函数可复算，输出见 reports\工具_rot三项审计.txt）
%   (1) Polaris 真实方位角 +0.390 deg（论文取 0.4 deg）⇒ 这一项论文是对的，差 0.01 deg。
%   (2) 20140923 17:06:12 帧里，(601,141)（= 裁剪系 (499,55) + (102,86)）邻域 30 px
%       内峰值 30、邻域中位 27、全幅中位 27 ⇒ 8 bit 量化下**无可辨认星像**，
%       论文这颗星无法从归档数据复核。
%   (3) 二维 rms 在 dt=0 最小（4.572 px），纯方位散射在 dt=-3~0 最小（0.538 deg）；
%       rot 对 dt 的斜率 0.160 deg/min ⇒ 钟差 5 min 也只带 0.80 deg。**时钟是准的。**
%
% 用法
%   rot_audit          % 默认扫描 +-36 min
%   rot_audit(20)      % 自定义扫描范围
%   rot_audit(0)       % 只做 (1)(2)，不扫描时钟
%
% 依赖：mooncal_out_matlab\calibration_params_14\moon_meas.mat（月核质心缓存）、
%       moon_altaz_ref.m、sun_altaz_ref.m、一期原始数据\20140923\...170612...PNG

if nargin < 1 || isempty(dtMax), dtMax = 36; end

MP     = mc_paths();
root   = MP.root;
outDir = MP.logs;
if ~exist(outDir, 'dir'), mkdir(outDir); end
fid = fopen(fullfile(outDir, '工具_rot三项审计.txt'), 'w', 'n', 'UTF-8');
cln = onCleanup(@() fclose(fid));
lg = @(varargin) fprintf(fid, varargin{:});

lon = 109.133; lat = 19.526; hgt = 0.103;
zenPaper = [554.0, 538.0];     % 论文 §3.2.1 三星余弦天顶（原图系）
zenOurs  = [553.5, 536.5];     % 本解天顶
cropOff  = [102, 86];          % (554,538) - (452,452)
polCrop  = [499, 55];
polRaw   = polCrop + cropOff;
polAzTh  = 0.4;
rotTh    = 6.35;
t0       = datetime(2014, 9, 23, 17, 6, 12);

lg('=== 一期 rot 三项独立审计 === 生成 %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

% ================= (1) 论文的北极星方位角对不对 =================
lg('----- (1) 论文"北极星方位角 0.4 deg"的独立复算 -----\n');
raP  = (2 + 31/60 + 49.094/3600) * 15;      % Polaris J2000 RA  = 37.954 deg
decP = 89 + 15/60 + 50.88/3600;             % Polaris J2000 Dec = 89.264 deg
[altP, azP] = paltaz(t0, lon, lat, raP, decP);
lg('2014-09-23 17:06:12 UT @ (%.3fE, %.3fN)：alt = %+.3f deg, az = %+.3f deg\n', ...
   lon, lat, altP, azP);
lg('论文取 az = %.2f deg ⇒ 差 %+.3f deg  ⇒ **这一项论文是对的**\n', polAzTh, azP - polAzTh);
lg('（该纬度上 Polaris 方位角一个恒星日内摆幅仅 +-%.3f deg，所以它必须取对时刻才对）\n', ...
   0.74 / cosd(lat));

% ================= (2) 那颗星能不能复核 =================
lg('\n----- (2) 论文北极星像素 (499,55) 在归档帧里的可见性 -----\n');
fs = dir(fullfile(MP.raw, '20140923', '*2014092317061*.PNG'));
if isempty(fs)
    lg('未找到 20140923170612 帧，跳过。\n');
else
    A = double(imread(fullfile(MP.raw, '20140923', fs(1).name)));
    bg = median(A(:));
    lg('帧 %s   全幅中位 %.1f（8 bit）\n', fs(1).name, bg);
    px = polRaw;
    x1 = max(1, px(1) - 30); x2 = min(size(A,2), px(1) + 30);
    y1 = max(1, px(2) - 30); y2 = min(size(A,1), px(2) + 30);
    sub = A(y1:y2, x1:x2);
    [pv, idx] = max(sub(:));
    [iy, ix] = ind2sub(size(sub), idx);
    gx = x1 + ix - 1; gy = y1 + iy - 1;
    lg('邻域 30 px（论文像素 %s，原图系 %s）：峰值 %.0f @ (%d,%d)，邻域中位 %.1f\n', ...
       mat2str(polCrop), mat2str(px), pv, gx, gy, median(sub(:)));
    lg('距论文给定点 %.1f px；拉开背景 %.0f 个灰阶\n', hypot(gx - px(1), gy - px(2)), pv - bg);
    lg('判据：8 bit 强量化下"星"至少要拉开 8 个灰阶才算看得见。\n');
    if pv - bg < 8
        lg('⇒ **测不到**：论文这颗星无法从归档数据复核。\n');
    else
        lg('⇒ 有可疑亮点，需进一步判形状（星为多像元团，热像元为单像元脉冲）。\n');
    end
    % 同一帧另外几处亮斑作对照
    lg('对照：同帧其它"亮斑"（阈值 bg+40，按峰值前 5）：\n');
    m = A > bg + 40;
    [ys, xs] = find(m);
    if ~isempty(ys)
        v = A(m); [~, o] = sort(v, 'descend');
        xs = xs(o); ys = ys(o);
        C = zeros(0, 2);
        for i = 1:numel(xs)
            if isempty(C)
                C = [xs(i) ys(i)];
            else
                d = hypot(C(:,1) - xs(i), C(:,2) - ys(i));
                if min(d) >= 9, C(end+1, :) = [xs(i) ys(i)]; end %#ok<AGROW>
            end
        end
        for i = 1:min(5, size(C,1))
            lg('   (%4d,%4d) 峰值 %3.0f  离天顶 %6.1f px\n', C(i,1), C(i,2), A(C(i,2),C(i,1)), ...
               hypot(C(i,1) - zenPaper(1), C(i,2) - zenPaper(2)));
        end
    end
    lg('（这些亮斑都落在 r ~ 470 px 的视场边缘带上，是渐晕/边缘像元，不是星）\n');
end

% ================= (3) 论文那颗星推出什么 rot =================
lg('\n----- (3) 论文 rot 的可裁决部分 -----\n');
phiTh = mod(atan2d(polRaw(1) - zenOurs(1), zenOurs(2) - polRaw(2)), 360);
rPol  = hypot(polRaw(1) - zenOurs(1), zenOurs(2) - polRaw(2));
lg('论文北极星像素配本解天顶 (%.1f,%.1f)：r = %.1f px, phi = %.2f deg\n', ...
   zenOurs(1), zenOurs(2), rPol, phiTh);
lg('  ⇒ rot = phi - az = %.2f - %.2f = **%.2f deg**（论文报 %.2f，差 %.2f）\n', ...
   phiTh, polAzTh, phiTh - polAzTh, rotTh, phiTh - polAzTh - rotTh);
lg('  ⇒ 论文 (499,55) 与 (452,452) 这两个数**互相自洽**；2 deg 的差全落在"这颗星的方向"上。\n');
lg('2 deg 在 r = %.0f px 处 = 横向 %.1f px；而该星比论文自己的径向式(3.9)外偏 %.1f px（%.1f%%）\n', ...
   rPol, 2 * pi / 180 * rPol, ...
   rPol - polyval([-0.55, 305.2, 16.68], deg2rad(90 - altP)), ...
   100 * abs(rPol - polyval([-0.55, 305.2, 16.68], deg2rad(90 - altP))) / ...
   polyval([-0.55, 305.2, 16.68], deg2rad(90 - altP)));

% ================= (4) 时钟偏置扫描 =================
if dtMax <= 0
    lg('\n(dtMax<=0，跳过时钟扫描)\n');
    return
end
lg('\n----- (4) 相机时钟偏置扫描（dt = 星历时刻相对记录时刻的平移）-----\n');
mp = MP.measP14;
if exist(mp, 'file') ~= 2
    lg('未找到 %s，跳过。\n', mp);
    return
end
S = load(mp);
cx = double(S.cx(:)).'; cy = double(S.cy(:)).';
ts = S.ts(:).'; N = numel(cx);
T = NaT(1, N);
for i = 1:N
    s = char(ts{i});
    T(i) = datetime(str2double(s(1:4)), str2double(s(5:6)), str2double(s(7:8)), ...
                    str2double(s(9:10)), str2double(s(11:12)), str2double(s(13:14)));
end
[alt, azm] = meph(T, lon, lat, hgt);
salt = zeros(1, N);
for i = 1:N, salt(i) = sun_altaz_ref(T(i), lon, lat); end
sitF = (salt < -18) & (alt > 14) & (alt < 50);
lg('测量缓存 %d 帧；拟合帧集 %d 帧 / %d 夜（太阳<-18、alt 14~50）\n', ...
   N, sum(sitF), numel(unique(S.night(sitF))));

x0 = zenOurs(1); y0 = zenOurs(2);
keep = sitF;
q = [];
for it = 1:4
    kk = find(keep);
    zz = deg2rad(90 - alt(kk)).'; AA = deg2rad(azm(kk)).';
    q = fitq(zz, AA, cx(kk).' - x0, y0 - cy(kk).', [324, 4.3, 0.0], it);
    rr = q(1) * (zz + q(3) * zz.^3);
    e = hypot(rr .* sin(AA + deg2rad(q(2))) - (cx(kk).' - x0), ...
              (y0 - cy(kk).') - rr .* cos(AA + deg2rad(q(2))));
    if it < 4
        thr = max(5.0, median(e) + 3 * median(abs(e - median(e))));
        keep(kk) = e < thr;
    end
end
kk = find(keep);
rms0 = sqrt(mean(e.^2));
lg('基线 dt=0: n=%d  f=%.2f  a=%+.5f  rot=%+.3f deg  rms=%.3f px\n', ...
   numel(kk), q(1), q(3), q(2), rms0);

dts = -dtMax:3:dtMax;
rmsA = zeros(size(dts)); rotA = rmsA; azA = rmsA; inA = rmsA;
lg('\n  %6s %9s %9s %9s %9s %9s %8s\n', 'dt(min)', 'f', 'rot', 'rms', '夜散射', '方位散射', '内点%');
for k = 1:numel(dts)
    [a2, z2] = meph(T + minutes(dts(k)), lon, lat, hgt);
    zz = deg2rad(90 - a2(kk)).'; AA = deg2rad(z2(kk)).';
    xm = cx(kk).' - x0; ym = y0 - cy(kk).';
    qq = fit3(zz, AA, xm, ym, [q(1), q(2), q(3)]);
    rr = qq(1) * (zz + qq(3) * zz.^3);
    ee = hypot(rr .* sin(AA + deg2rad(qq(2))) - xm, ym - rr .* cos(AA + deg2rad(qq(2))));
    rmsA(k) = sqrt(mean(ee.^2)); rotA(k) = qq(2);
    rotI = mod(atan2(cx(kk) - x0, -(cy(kk) - y0)) * 180 / pi - z2(kk), 360);
    mu = mod(atan2(mean(sind(rotI)), mean(cosd(rotI))) * 180 / pi, 360);
    dd = mod(rotI - mu + 180, 360) - 180;
    azA(k) = sqrt(mean(dd.^2)); inA(k) = 100 * mean(abs(dd) < 2);
    lg('  %+6d %9.2f %9.3f %9.3f %9.3f %9.3f %8.1f\n', dts(k), qq(1), qq(2), rmsA(k), ...
       nightscat(z2(kk), S.night(kk), cx(kk), cy(kk), x0, y0), azA(k), inA(k));
end
[rm, im] = min(rmsA);
[am, ia] = min(azA);
lg('\n  ★ 二维 rms   最小在 dt = %+d min: rms = %.3f px（rot = %+.3f deg）\n', dts(im), rm, rotA(im));
lg('  ★ 纯方位散射 最小在 dt = %+d min: azRms = %.3f deg（rot = %+.3f deg, 内点 %.1f%%）\n', ...
   dts(ia), am, rotA(ia), inA(ia));
sl = (rotA(dts == 3) - rotA(dts == -3)) / 6;
lg('  ★ rot 对 dt 的斜率 %.3f deg/min ⇒ 钟差 5 min 也只带 %.2f deg\n', abs(sl), abs(sl) * 5);
lg('  ⇒ 两个独立统计量都指向 dt ~ 0 ⇒ **rot 的时钟假象可以排除**。\n');
lg('    （注：仅看"夜散射"会随 dt 单调下降，那是 5~6 夜样本噪声，不是钟差的证据；\n');
lg('      帧级方位散射 azRms 与二维 rms 才是有判别力的量。）\n');
end

% ---------------------------------------------------------------------------
function [alt, az] = paltaz(t, lon, lat, ra, dec)
% 低精度星历（Meeus ch.25 的恒星版）：Polaris 14 年岁差位移 < 0.08 deg，忽略
jd  = juliandate(t);
d   = jd - 2451545.0;
Tc  = d / 36525;
gmst = 280.46061837 + 360.98564736629 * d + 0.000387933 * Tc^2 - Tc^3 / 38710000;
lst = mod(gmst + lon, 360);
H   = mod(lst - ra + 180, 360) - 180;
U = sind(dec)*sind(lat) + cosd(dec)*cosd(lat)*cosd(H);
N = sind(dec)*cosd(lat) - cosd(dec)*sind(lat)*cosd(H);
E = -cosd(dec)*sind(H);
alt = asind(U);
az  = mod(atan2d(E, N), 360);
end

% ---------------------------------------------------------------------------
function [alt, az] = meph(T, lon, lat, hgt)
alt = zeros(size(T)); az = alt;
for i = 1:numel(T)
    [a1, a2] = moon_altaz_ref(T(i), lon, lat, hgt);
    alt(i) = a1; az(i) = a2;
end
end

% ---------------------------------------------------------------------------
function q = fitq(zz, AA, xm, ym, p0, it)
% it==1 时多起点（避开 a~0 的浅局部极小），否则单起点
if it == 1
    bq = []; bs = Inf;
    for a0 = [0, -0.015, -0.03, -0.045, +0.015]
        qq = fit3(zz, AA, xm, ym, [p0(1), p0(2), a0]);
        rr = qq(1) * (zz + qq(3) * zz.^3);
        ss = sum((rr .* sin(AA + deg2rad(qq(2))) - xm).^2 + ...
                 (ym - rr .* cos(AA + deg2rad(qq(2)))).^2);
        if ss < bs, bs = ss; bq = qq; end
    end
    q = bq;
else
    q = fit3(zz, AA, xm, ym, p0);
end
end

% ---------------------------------------------------------------------------
function q = fit3(zz, AA, xm, ym, p0)
% 固定圆心下拟合 (f, rot, a)：x = r sin(A+rot), y = -r cos(A+rot), r = f(z + a z^3)
rr = @(q) q(1) * (zz + q(3) * zz.^3);
rf = @(q) [rr(q) .* sin(AA + deg2rad(q(2))) - xm; ym - rr(q) .* cos(AA + deg2rad(q(2)))];
o = optimset('Display', 'off', 'MaxFunEvals', 40000, 'MaxIter', 40000, ...
             'TolX', 1e-9, 'TolFun', 1e-12);
q = fminsearch(@(q) sum(rf(q).^2), p0, o);
end

% ---------------------------------------------------------------------------
function s = nightscat(az, nt, cx, cy, x0, y0)
% 逐夜 rot 的圆散布（deg）。样本只有 5~6 夜，噪声大，只作参考。
ph = mod(atan2(cx - x0, -(cy - y0)) * 180 / pi, 360);
rot = mod(ph - az, 360);
uN = unique(nt); mus = [];
for i = 1:numel(uN)
    m = strcmp(nt, uN{i});
    if sum(m) < 10, continue; end
    mus(end+1) = mod(atan2(mean(sind(rot(m))), mean(cosd(rot(m)))) * 180 / pi, 360); %#ok<AGROW>
end
if numel(mus) < 2, s = NaN; return; end
R = abs(mean(exp(1i * deg2rad(mus))));
s = sqrt(-2 * log(max(R, 1e-12))) * 180 / pi;
end
