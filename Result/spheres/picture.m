function picture(front_size,px,py,l)
    if nargin < 4 l=1; end
    
    set(gca,'FontName','Times New Roman','FontSize',front_size);
 
    set(get(gca,'XLabel'),'FontName','Times New Roman','FontSize',front_size);%图上文字为8 point或小5号
    set(get(gca,'YLabel'),'FontName','Times New Roman','FontSize',front_size);
    %set(gcf,'FontName','Times New Roman');
    if l==0
        set(gcf,'unit','centimeters','position',[5 3 px py ]);
        %saves(gcf)
    end
end
