clear
clc
data = xlsread('Sum_r.xlsx','Sheet2');
data1 = xlsread('Sum_r.xlsx','Sheet3');
N = data(2:end,1);
Data = [data(2:end,2:end),data1(2:end,7:end)];
temp = Data
marker = {'r.-','ro-','rx-','r+-','r*-','bs-','bd-','kv-','k^-','bs:','bd:','kv:','k^:'};
sD = size(Data);
subplot(3,2,[1 3 5])
for i = 1:sD(2)
    x = N;
    y = Data(:,i)./0.07330;
    loglog(x,y,char(marker(i)))
    hold on
end
legend('SBB-MRT','LIBB-MRT','QIBB-MRT','MR-MRT',...
    'CLI-MRT','PSM-MRT-A','PSM-MRT-B','IBM-MRT-A','IBM-MRT-B',...
    'PSM-SRT-A','PSM-SRT-B','IBM-SRT-A','IBM-SRT-B',...
    'location','NorthEast')
grid on
% axis([0.5 2 0.7 1.3])
subplot(3,2,2)
for i=1:5
    x = N;
    y = Data(:,i)./0.07330;
    loglog(x,y,char(marker(i)))
    hold on
end
grid on
% axis([0.5 2 0.7 1.3])
subplot(3,2,4)
num = [6,7,10,11];
for i=1:numel(num)
    x = N;
    y = Data(:,num(i))./0.07330;
    loglog(x,y,char(marker(num(i))))
    hold on
end
grid on
% axis([0.5 2 0.7 1.3])
subplot(3,2,6)
num = [8,9,12,13];
for i=1:numel(num)
    x = N;
    y = Data(:,num(i))./0.07330;
    loglog(x,y,char(marker(num(i))))
    hold on
end
grid on
% axis([0.5 2 0.7 1.3])