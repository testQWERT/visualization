% By lyqmath
% DLUT School of Mathematical Sciences 2008
% BLOG��http://blog.sina.com.cn/lyqmath
function main()
global hfig
hfig = figure;
hold on; box on;
global cc;

aa=imread('blank.png');
cc=imshow(aa);
haxis = gca;
set(cc, 'ButtonDownFcn', @click_ceshi);
function click_ceshi(src, event)
global hfig
global cc
% ��ȡ��ǰ������
xy = get(cc, 'CurrentPoint');
% ��ȡgcf��gca��λ����Ϣ
hpos = get(hfig, 'Position');
apos = get(gca, 'Position');
% ���µ����Ϣ����ȡ�������gca��׼ȷ������Ϣ
x = (xy(1) - apos(1)*hpos(3))/(apos(3)*hpos(3));
y = (xy(2) - apos(2)*hpos(4))/(apos(4)*hpos(4));
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
x = x*(xlim(2) - xlim(1)) + xlim(1);
y = y*(ylim(2) - ylim(1)) + ylim(1);
% ��ע
text(x, y, 'Happy', 'color', rand(3, 1));