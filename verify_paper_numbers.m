function acc = verify_paper_numbers()
%VERIFY_PAPER_NUMBERS  逐条核对论文里的数值与随包标定文件是否一致。
%
%   verify_paper_numbers            % 逐条打印 + 末尾汇总
%   acc = verify_paper_numbers()    % acc.n / acc.bad / acc.miss
%
% 做什么
%   论文表 1（三个时期的采用解）、式 (5)（Δrot/Δy0 的斜率）、以及"剔 20131021 后
%   合 10 夜 rot = 4.283° ± 0.204°"这条跨期结论，全部回溯到包内 data/calib/*.mat。
%   断言口径 = **论文的印刷精度**（例如表 1 的 x0 印到 2 位小数 ⇒ 容差 0.005），
%   所以本脚本回答的是"论文里印的那个数，能不能从随包文件逐位复算出来"。
%
% ★ 本脚本只读 data/calib/ 下的 .mat，不读原始帧，因此**不需要 240 GB 存档数据**
%   就能跑。这是本包对"审稿阶段可复现"的核心承诺。
%
% 注意
%   · 表 2（与马欣 2021 的对照）与表 3（与文献漂移速度对照）不在此断言范围内：
%     前者来自 calib_2014 的 shape_cmp / 三星余弦对照，后者来自文献，不是 .mat 里的量。
%     表 2 相关的量在本脚本末尾以"参考输出"打印，供人工对照。
%   · 式 (6) 的符号约定与代码不一致（见 docs/公式与代码对照.md 的 ⚠），本脚本不覆盖。

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % <pkg>
addpath(fullfile(root, 'src'));
MP = mc_paths();

acc.n = 0; acc.bad = 0; acc.miss = 0;

fprintf('\n%s\n', repmat('=', 1, 78));
fprintf('  论文数值 <- 随包标定文件  逐条核对\n');
fprintf('  数据源 %s\n', fullfile(MP.calib, '*.mat'));
fprintf('%s\n', repmat('=', 1, 78));

% ==========================================================================
%  表 1  三个时期的采用解
% ==========================================================================
T = { ...
  '二期 2026-03',  MP.calP2Native; ...
  '一期 2014-09',  MP.calP14Native; ...
  '一期 2013-10',  MP.calP13Native};

EXPECT = {
%   x0       y0       f        a          rot      rms     nFrm  nNgt
    538.87,  476.64,  422.24, -0.0998,    0.39,    2.11,   349,  5;   % 二期
    553.50,  536.50,  324.42, -0.0330,    4.35,    4.57,   316,  6;   % 一期 2014
    553.50,  536.50,  333.90, -0.0385,    4.20,    2.01,   194,  5};  % 一期 2013

TOL = [0.005, 0.005, 0.005, 0.00005, 0.005, 0.005, 0, 0];

for k = 1:size(EXPECT, 1)
    fprintf('\n--- 表 1  %s   [%s]\n', T{k, 1}, relpath(MP, T{k, 2}));
    S = load(T{k, 2});
    if ~isfield(S, 'calMat')
        fprintf('  [!!] 文件里没有 calMat：%s\n', T{k, 2}); acc.miss = acc.miss + 1; continue
    end
    m = S.calMat;
    vals = [m.x0, m.y0, m.f, m.a, m.rot, m.rms_px, m.n_frames, m.n_nights];
    nm   = {'天顶像点 x0', '天顶像点 y0', '径向标度 f', '非线性 a', ...
            '方位定向角 rot', '二维残差 rms', '帧数', '夜数'};
    for j = 1:numel(vals)
        acc = chk(acc, sprintf('%s %s', T{k, 1}, nm{j}), vals(j), EXPECT{k, j}, TOL(j));
    end
    if isfield(m, 'z_rim_deg')
        fprintf('       参考：盘边天顶角 z_rim = %.2f°（论文 §3.3 讨论用，非表 1 项）\n', m.z_rim_deg);
    end
end

% ==========================================================================
%  式 (5)  Δrot / Δy0
% ==========================================================================
fprintf('\n--- 式 (5)  Δrot/Δy0（圆心偏 1 px 引起的 rot 偏置）\n');
S14 = load(MP.calP14Native); S13 = load(MP.calP13Native);
acc = chk(acc, '式(5) 2014-09  Δrot/Δy0', S14.calMat.rot_per_y0, +0.14, 0.005);
acc = chk(acc, '式(5) 2013-10  Δrot/Δy0', S13.calMat.rot_per_y0, -0.17, 0.005);

% ==========================================================================
%  表 1 的"钉死"标注与跨期结论
% ==========================================================================
fprintf('\n--- 表 1 脚注 / §3.2 跨期结论\n');
if isfield(S13.calMat, 'center_source')
    hasFix = ~isempty(strfind(S13.calMat.center_source, 'CenterFix'));  %#ok<STREMP>
    acc = chklogical(acc, '2013 的圆心由 2014 解约束（表 1 "钉死"）', hasFix, true);
else
    fprintf('  [!!] 2013 解缺 center_source 字段，无法确认"钉死"\n'); acc.miss = acc.miss + 1;
end
dRot = abs(S13.calMat.rot - S14.calMat.rot);
acc = chk(acc, '2013 与 2014 的 rot 之差', dRot, 0.15, 0.005);

% 剔 20131021 后合 10 夜
[r10, s10, n10] = merge_rot(S14.calMat, S13.calMat, '20131021');
acc = chk(acc, '合 10 夜（剔 20131021）rot', r10, 4.283, 0.0005);
acc = chk(acc, '合 10 夜（剔 20131021）圆散布', s10, 0.204, 0.0005);
fprintf('       合 10 夜共 %d 帧（论文 §3.2 写 498 帧）\n', n10);

% ==========================================================================
%  表 2（马欣 2021 对照）—— 只打印，不硬断言
% ==========================================================================
fprintf('\n--- 表 2 参考输出（与马欣 2021 对照；来自 calib_2014 的 shape_cmp，非本包 .mat 字段）\n');
fprintf('       本解天顶像点 (%.1f, %.1f)；r(90°) = %.1f px\n', ...
        S14.calMat.x0, S14.calMat.y0, S14.calMat.r_horizon_px);
fprintf('       论文印的是 (553.5, 536.5) / 468.1 px\n');

% ==========================================================================
fprintf('\n%s\n', repmat('=', 1, 78));
if acc.bad == 0 && acc.miss == 0
    fprintf('  全部 %d 条通过：论文数值可由随包文件逐位复算。\n', acc.n);
else
    fprintf('  %d 条通过，%d 条不符，%d 条缺字段 —— 请核对上表。\n', acc.n - acc.bad, acc.bad, acc.miss);
end
fprintf('%s\n\n', repmat('=', 1, 78));
end

% ==========================================================================
function acc = chk(acc, name, got, want, tol)
acc.n = acc.n + 1;
if abs(got - want) <= tol
    fprintf('  [OK]  %-32s 复算 %-11s 论文 %-11s\n', name, num2str(got, '%.6g'), num2str(want, '%.6g'));
else
    fprintf('  [!!]  %-32s 复算 %-11s 论文 %-11s  ← 差 %.4g（容差 %.4g）\n', ...
            name, num2str(got, '%.6g'), num2str(want, '%.6g'), abs(got - want), tol);
    acc.bad = acc.bad + 1;
end
end

function acc = chklogical(acc, name, got, want)
acc.n = acc.n + 1;
if isequal(got, want)
    fprintf('  [OK]  %-32s %s\n', name, tf(got));
else
    fprintf('  [!!]  %-32s 复算 %s 期望 %s\n', name, tf(got), tf(want));
    acc.bad = acc.bad + 1;
end
end

function [r, s, n] = merge_rot(c14, c13, dropNight)
%MERGE_ROT 把两期的逐夜 rot 合起来算圆均值与圆散布（可选剔除某一夜）。
rot = []; night = {}; n = 0;
for c = {c14, c13}
    cc = c{1};
    if ~isfield(cc, 'rot_per_night'), continue; end
    pn = cc.rot_per_night;
    for i = 1:numel(pn)
        night{end + 1} = pn(i).night;   %#ok<AGROW>
        rot(end + 1)   = pn(i).rot;     %#ok<AGROW>
        n = n + pn(i).n;
    end
end
if ~isempty(dropNight)
    m = strcmp(night, dropNight);
    rot(m) = []; n = n - sum_removed(c14, c13, dropNight);
end
[r, s] = circ(rot);
end

function n = sum_removed(c14, c13, night)
n = 0;
for c = {c14, c13}
    cc = c{1};
    if ~isfield(cc, 'rot_per_night'), continue; end
    pn = cc.rot_per_night;
    for i = 1:numel(pn)
        if strcmp(pn(i).night, night), n = n + pn(i).n; end
    end
end
end

function [m, s] = circ(rotDeg)
v = deg2rad(rotDeg(:)');
R = abs(mean(exp(1i * v)));
R = min(max(R, 1e-12), 1 - 1e-12);
m = mod(atan2(mean(sind(rotDeg(:))), mean(cosd(rotDeg(:)))) * 180 / pi, 360);
s = sqrt(-2 * log(R)) * 180 / pi;
end

function p = relpath(P, f)
p = strrep(f, [P.root filesep], '');
end

function s = tf(v)
if v, s = '是'; else, s = '否'; end
end
