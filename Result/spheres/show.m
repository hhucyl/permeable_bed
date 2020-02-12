clear
clc
data = xlsread('Sum_p1.xlsx');
data1 = xlsread('Sum_p1.xlsx','Sheet2');
tau = data(2:end,1);
nu = (tau-0.5)./3;
Data = [data(2:end,2:10),data1(2:end,7:end)];
dref1 = mean(data(2:end,end));
dref2 = mean(data(:,end-1));
marker = {'g.-','go-','gx-','g+-','g*-','bs-','bd-','kv-','k^-','bs:','bd:','kv:','k^:'};
markersize = 10;
sD = size(Data);
subplot(3,2,[1 3 5])
nn = [1:5]
for i = 1:numel(nn)
    x = nu;
    y = Data(:,nn(i))./dref1;
    if(i==1) mt = 3*markersize; else mt = markersize; end
    plot(x,y,char(marker(i)),'markersize',mt)
    hold on
    drawnow
end
legend_str = {'SBB-MRT','LIBB-MRT','QIBB-MRT','MR-MRT',...
    'CLI-MRT','PSM-MRT-A','PSM-MRT-B','IBM-MRT-A','IBM-MRT-B',...
    'PSM-SRT-A','PSM-SRT-B','IBM-SRT-A','IBM-SRT-B'
    };
legend(legend_str,'location','southeast');
legend boxoff
% grid on
axis([min(x) max(x) 0.3 1.5])
set(gca,'ytick',0.3:0.1:1.5)
xlabel('Viscosity')
ylabel('\itk^{*}_{\rmsimulated}/k^{*}_{\rmanalytical}')
picture(15)
subplot(3,2,2)
for i=1:5
    x = nu;
    y = Data(:,i)./0.07330;
    if(i==1) mt = 3*markersize; else mt = markersize; end

    h1(i,1) = plot(x,y,char(marker(i)),'markersize',mt);
    
    hold on
end
% grid on
axis([min(x) max(x) 0.7 1.3])
xlabel('Viscosity')
ylabel('\itk^{*}_{\rmsimulated}/k^{*}_{\rmanalytical}')
legend(h1(1:2),legend_str{1:2},'orientation','horizontal','location','north');
legend boxoff
picture(14)
ah = axes('position',get(gca,'position'),'visible','off');
legend(ah,h1(3:5),legend_str{3:5},'orientation','horizontal','location','south');
legend boxoff
picture(14)
% subplot(3,2,4)
% num = [6,7,10,11];
% for i=1:numel(num)
%     x = nu;
%     y = Data(:,num(i))./0.07330;
%     plot(x,y,char(marker(num(i))),'markersize',markersize)
%     hold on
% end
% legend(legend_str{num},'orientation','horizontal','location','south')
% legend boxoff
% % grid on
% axis([min(x) max(x) 0.7 1.3])
% xlabel('Viscosity')
% ylabel('\itk^{*}_{\rmsimulated}/k^{*}_{\rmanalytical}')
% picture(14)
% subplot(3,2,6)
% num = [8,9,12,13];
% for i=1:numel(num)
%     x = nu;
%     y = Data(:,num(i))./0.07330;
%     plot(x,y,char(marker(num(i))),'markersize',markersize)
%     hold on
% end
% % grid on
% axis([min(x) max(x) 0.7 1.3])
% xlabel('Viscosity')
% ylabel('\itk^{*}_{\rmsimulated}/k^{*}_{\rmanalytical}')
% legend(legend_str{num},'orientation','horizontal','location','south')
% legend boxoff
% picture(14)