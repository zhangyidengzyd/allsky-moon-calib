function r = pkg_root()
%PKG_ROOT  包根目录（= src/ 的上一级）。整包的锚点，只有 5 行。
%
%   r = pkg_root()
%
% 任何文件里都不要出现盘符；一律走 pkg_root() / mc_paths()，
% 这样把整个文件夹拷到任何位置都能跑。

here = fileparts(mfilename('fullpath'));
r = fileparts(here);          % src 的上一级 = 包根
end
