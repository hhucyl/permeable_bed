clear
clc
file_name = {'Cd_collidetype1_bbtype1_100_0.8.out'};
file = importdata(char(file_name));
loglog(file.data(:,2),file.data(:,3),'g')
hold on 

loglog(file.data(:,2),file.data(:,4),'r')

loglog(file.data(:,2),file.data(:,5),'o')

legend({'$\frac{24}{Re}+\frac{6}{1+\sqrt{Re}}+0.4$','$\frac{24}{9.06^2}(\frac{9.06}{\sqrt{Re}}+1)^2$','Stimulation'},...
    'Interpreter','latex','fontsize',15,'location','SouthWest')
legend boxoff