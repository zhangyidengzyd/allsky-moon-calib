%MOON_CALIB_SELFTEST  交付前自检：字段名 / 正反投影往返 / 地平半径。
%
%   moon_calib_selftest          % 在 MATLAB 里直接运行本脚本
%
% 检查三件事（全部只读 data/calib/phase2_native.mat，不需要原始帧）
%   ① calMat 的字段是否齐全（下游 moon_calib_to_pipeline 依赖这些字段名）
%   ② 正向投影 -> 反解的往返误差（应为数值零量级）
%   ③ r(90°) 是否与论文报告值 499.9 px 一致
%
% ★ 本脚本在 src/verify/ 下 ⇒ root 要自己数三层，不要直接调 pkg_root()
%   （pkg_root 是 src/ 下的锚点，从 verify/ 调用会算错一层）。

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));          % <pkg>
addpath(here);                              % src/verify/
addpath(fullfile(root, 'src'));             % 依赖 src/ 下的 moon_calib_project / _unproject
MP = mc_paths();

S = load(MP.calP2Native);
cal = S.calMat;

fprintf('== calMat fields ==\n');
fn = fieldnames(cal);
fprintf('%s ', fn{:}); fprintf('\n');
fprintf('x0=%.3f y0=%.3f f=%.4f a=%.6f rot=%.4f lat=%.3f lon=%.3f\n', ...
    cal.x0, cal.y0, cal.f, cal.a, cal.rot, cal.site_lat, cal.site_lon);

r90 = cal.f * (pi/2 + cal.a * (pi/2)^3);
fprintf('\n== horizon radius ==  r(90deg)=%.2f px   (论文 §3.1 报 499.9 px)\n', r90);

lat = cal.site_lat;
tests = [ 10 30 ; 20 -40 ; 45 75 ];
fprintf('\n== round trip ==\n');
fprintf('  dec      H     z(deg)   az(deg)        x         y  |  z2(deg)  az2(deg)      dz      daz       dr\n');
for k = 1:size(tests,1)
    dec = tests(k,1); Hd = tests(k,2);
    num = -sind(Hd) * cosd(dec);
    den = sind(dec) * cosd(lat) - cosd(dec) * sind(lat) * cosd(Hd);
    az  = mod(rad2deg(atan2(num, den)), 360);
    [x, y, zd, rr] = moon_calib_project(cal, lat, dec, Hd, az);
    [z2, az2, r2]  = moon_calib_unproject(cal, x, y);
    fprintf('%6.1f %6.1f %8.4f %8.4f %10.3f %9.3f | %8.4f %8.4f %8.2e %8.2e %8.2e\n', ...
        dec, Hd, zd, az, x, y, z2, az2, z2-zd, az2-az, r2-rr);
end

fprintf('\n== origin sanity == zenith pixel (%.2f, %.2f) is 0-based; +1 for MATLAB indexing\n', cal.x0, cal.y0);
disp('ROUNDTRIP_OK');
