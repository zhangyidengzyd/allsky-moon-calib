function [g, info] = q_read_gray(path)
%Q_READ_GRAY  读一帧灰度 PNG，按**位深**统一到 16 bit 量级（返回 double, 0..65535）。
%
%   g = q_read_gray(path)
%   [g, info] = q_read_gray(path)     % info.bits / info.scale / info.class
%   g = q_read_gray([])               % -> []（调用方自己计数，绝不静默顶替）
%
% 为什么需要它
% ------------
% 本工程的帧有两种位深，而链上**所有阈值都是按 16 bit 写的**：
%
%   二期 2025-2026 : 1024x1024, 16 bit（满量程 65535）；夜天光 ~3400~8000 DN
%   一期 2014      : 1024x1024,  8 bit（满量程   255）；夜天光 ~13~60   灰阶
%
% 典型阈值：饱和判据 65000 DN（moon_calib 的 C.sat_level）、视场边界对比 50 DN
% （zenith_from_rim）、暗场 3441 DN 基座（bubble_night / clarity_rank 要减它）。
% 若直接 imread，8 bit 数据**永远**满足不了 `> 65000` ⇒ moon_calib 的 detect
% 阶段一帧都检不出，而报错信息会把它伪装成"数据里没有月亮"，把人往错方向带。
%
% 处理方式（只做位深对齐，不改内容）
% ---------------------------------
%   uint16 -> double(raw)          原样（16 bit 是基准）
%   uint8  -> double(raw) * 257    线性提升到 16 bit 量级（255*257 = 65535，满量程对齐）
%   double -> 原样；若 max>1 视为已是 DN，否则按 [0,1] 映射到 65535
%   整型且为 3 通道 -> 取第 1 通道
%
% ⚠ 这个函数**只修位深，不修动态范围**。8 bit 一期数据夜天光只有十几到几十灰阶，
%   百分之几的气辉耗竭在那条数据上是被量化吃掉的 —— 位深对齐不会把它变回来。
%   详见 reports\泡团观测_运行手册.md 的"一期数据可用性"一节。
%
% 失败（文件坏 / 不是图像）返回 g = []，不抛错：
%   归档帧有整夜截断损坏的情况（IHDR 完好、IDAT 不全），调用方需要按夜汇总计数。

g = [];
info = struct('bits', NaN, 'scale', 1, 'class', '');
if isempty(path), return; end
if isstring(path), path = char(path); end

try
    im = imread(path);
catch
    g = [];
    return
end
if isempty(im), return; end
if ndims(im) == 3, im = im(:, :, 1); end   % 兼容灰度/彩色两种写盘

info.class = class(im);
if isa(im, 'uint16')
    info.bits = 16;  info.scale = 1;
    g = double(im);
elseif isa(im, 'uint8')
    info.bits = 8;   info.scale = 257;
    g = double(im) * 257;
else
    g = double(im);
    if isempty(g), return; end
    if max(g(:)) > 1
        info.bits = 16;  info.scale = 1;
    else
        info.bits = 1;   info.scale = 65535;
        g = g * 65535;
    end
end
end
