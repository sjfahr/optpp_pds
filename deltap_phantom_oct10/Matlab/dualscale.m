%scatter(x(2:60),rnoise*ones(1,59),'r--','LineWidth', 2)


plot(x(2:60),RODS(2:60),'r','LineWidth', 2)
hold
hline=refline(0,rnoise);
set(hline,'Color','r','LineWidth',2,'LineStyle', :)
set(gca,'Box','off')
ylim([0,.75])
ylabel('RMS Difference (^oC)', 'FontSize', 20)
axes('YAxisLocation','right' ,'YLim',[0,1.77], 'XTick',[],'Color', 'none')
ylabel ('Relative Error', 'FontSize', 20)
