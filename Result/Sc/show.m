clear
clc
data = xlsread('Sum.xlsx');
data1 = xlsread('Sum.xlsx','Sheet2');
tau = data(5:end,1);
Data = [data(5:end,2:end),data1(5:end,7:end)];
marker = {'r.','ro','rx','r+','r*','rs','rd','rv','r^','ks','kd','kv','k^'};
sD = size(Data);
subplot(3,2,[1 3 5])
for i = 1:sD(2)
    x = tau;
    y = Data(:,i)./0.07330;
    plot(x,y,char(marker(i)))
    hold on
end
legend('0-MRT','1-MRT','2-MRT','3-MRT',...
    '4-MRT','-1-MRT','-2-MRT','-3-MRT',...
    '-4-MRT','-1-SRT','-2-SRT','-3-SRT','-4_SRT','location','SouthEast')
grid on
axis([0.5 2 0.7 1.3])
subplot(3,2,2)
for i=1:5
    x = tau;
    y = Data(:,i)./0.07330;
    plot(x,y,char(marker(i)))
    hold on
end
grid on
axis([0.5 2 0.7 1.3])
subplot(3,2,4)
num = [6,9,10,13];
for i=1:numel(num)
    x = tau;
    y = Data(:,num(i))./0.07330;
    plot(x,y,char(marker(num(i))))
    hold on
end
grid on
axis([0.5 2 0.7 1.3])
subplot(3,2,6)
num = [7,8,11,12];
for i=1:numel(num)
    x = tau;
    y = Data(:,num(i))./0.07330;
    plot(x,y,char(marker(num(i))))
    hold on
end
grid on
axis([0.5 2 0.7 1.3])