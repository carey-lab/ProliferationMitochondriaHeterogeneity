%% load viability data
WORKDIR  = '~/Develop/ProliferationMitochondriaHeterogeneity/' ;
VIAB = [ WORKDIR 'Data/2020__FastSlowSort__viability_data.xlsx' ] ;

T = readtable( VIAB , 'Sheet' , 'Viability');
T = T( ~regexpcmp(T.Drug,'\+'),:);


T.Speed = categorical(T.Speed);
T.Drug = categorical(T.Drug);
T.CellType = categorical(T.CellType);

U = stack(T,T.Properties.VariableNames(3:end-2)) ;
U.Properties.VariableNames{5} = 'Day';
U.Properties.VariableNames{6} = 'Viab';
U.Day = cellfun(@(X)str2double(regexprep(X,'Day','')),string(U.Day));
%
U = U(U.Day>0,:)

figure;
clrs = get(gca,'ColorOrder');
boxplot(U.Viab, {U.CellType,U.Drug,U.Speed},'Colors',clrs([1 2],:),'FactorGap',[10,10,0])
ylabel('Viability')

%%

%%
figname = '~/Downloads/viability_boxplts_' ;

for drugname = {'Antimycin' 'Oligomycin' 'DMSO'}
    drugname = drugname{:}; 
    for celltype = {'ESC' 'FIB'}
    
    celltype = celltype{:};
    
    Q = U(U.Drug==drugname & U.CellType == celltype , : );
    [~,p] = ttest2(Q.Viab(Q.Speed=='Fast cells'),Q.Viab(Q.Speed=='Slow cells'));

    figure('Position',[100 100 200 200]);
    hold on ;
    clrs = get(gca,'ColorOrder');
    bh = boxplot(Q.Viab,Q.Speed,'widths',0.8);
    ylabel('Viability    (% BrdU^+)');
    title( sprintf('%s +%s p=%0.04f',celltype, drugname , p ) )
    line(xlim,[0 0],'LineStyle','--','Color',[.7 .7 .7])
    
    h = findobj(gca,'Tag','Box');
    patch(get(h(1),'XData'),get(h(1),'YData'),clrs(2,:),'FaceAlpha',1);
    patch(get(h(2),'XData'),get(h(2),'YData'),clrs(1,:),'FaceAlpha',1);
   
    xlabel('CFSE sort')
    ylim([50,100])
    
    print('-dpng',[figname '_box_' drugname '_' celltype '.png'] , '-r600')
    close;
    
    end
end
