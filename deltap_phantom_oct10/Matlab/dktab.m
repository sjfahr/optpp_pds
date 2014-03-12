figure
for nn=1:max(dakota_tabular(:,5))
    b=0
    x=0
    y=0
    for mm=1:numel(dakota_tabular(:,1))
    if dakota_tabular(mm,5)==nn
        b=b+1
        x(b)=dakota_tabular(mm,2)
        y(b)=dakota_tabular(mm,3)
    end
    end
    scatter(x,y,'o')
    hold all
    k=waitforbuttonpress
    if k==1
    clf
    end
end
    
    