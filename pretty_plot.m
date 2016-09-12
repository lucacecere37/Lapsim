set(gcf,'Units','inches')
set(gcf,'Position',[2 2 8 6])
set(gcf,'PaperPosition',[2 2 8 6])
set(gcf,'color','w')
ax = findall(gcf,'type','axes');

for n = 1:length(ax)
    set(get(ax(n),'Title'),'fontweight','bold')
    set(get(ax(n),'Xlabel'),'fontweight','bold','fontsize',14)
    set(get(ax(n),'Ylabel'),'fontweight','bold','fontsize',14)
    set(get(ax(n),'Zlabel'),'fontweight','bold','fontsize',14)
    set(get(ax(n),'Title'),'fontsize',16)
    
    try
        set(get(ax(n),'Children'),'linewidth',2)
    end
    
    try
        set(get(ax(n),'Children'),'markersize',18)
    end
    
    grid on;
end
