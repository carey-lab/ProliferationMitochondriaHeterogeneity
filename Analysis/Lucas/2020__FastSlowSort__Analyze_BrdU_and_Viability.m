%% load data
WORKDIR  = '~/Develop/GrowthPredictionTMRE/' ;
BRDU = [ WORKDIR 'Data/2020__FastSlowSort__Brdu_data.xlsx' ] ;
VIAB = [ WORKDIR 'Data/2020__FastSlowSort__viability_data.xlsx' ] ;

%% look at viability
F1 = readtable( VIAB , 'Sheet' , 'Fibro replicate 1');
F2 = readtable( VIAB , 'Sheet' , 'Fibro replicate 2');
E1 = readtable( VIAB , 'Sheet' , 'ESC replicate 1');
E2 = readtable( VIAB , 'Sheet' , 'ESC replicate 2');

F = F1 ; 
F(:,3:end) = array2table( (table2array(F1(:,3:end)) + table2array(F2(:,3:end)) ) ./ 2  ) ;

E = E1 ; 
E(:,3:end) = array2table( (table2array(E1(:,3:end)) + table2array(E2(:,3:end)) ) ./ 2  ) ;

%% look at BrdU
E = readtable( BRDU , 'Sheet' , 'BrdU_ESC');
F = readtable( BRDU , 'Sheet' , 'BrdU_Fibro');
E = E( ~isnan(E.Replicate1),:);
F = F( ~isnan(F.Replicate1),:);
E.CellsInS = mean([E.Replicate1,E.Replicate2],2) ;
F.CellsInS = mean([F.Replicate1,F.Replicate2],2) ;

E.Type = regexprep(E.Type,' cells','');
F.Type = regexprep(F.Type,' cells','');

F.Speed = categorical( regexprep(F.Type,' .*','') );
E.Speed = categorical( regexprep(E.Type,' .*','') );
F.Drug = categorical( regexprep(F.Type,'.* ','') );
E.Drug = categorical( regexprep(E.Type,'.* ','') );
E.Day = categorical(E.Day);
F.Day = categorical(F.Day);

%% plot averages FIB
figure; 

tiledlayout(2,3)

nexttile;
gscatter(F.O_N(1:4),F.O_N(5:8),F.Drug(1:4),'kbgr',[],30);
xlabel('FAST')
ylabel('SLOW')
line(xlim,xlim)
title('FIB O/N')

nexttile;
gscatter(F.Day2(1:4),F.Day2(5:8),F.Drug(1:4),'kbgr',[],30);
xlabel('FAST')
ylabel('SLOW')
line(xlim,xlim)
title('FIB Day 2')

nexttile;
gscatter(F.Day3(1:4),F.Day3(5:8),F.Drug(1:4),'kbgr',[],30);
xlabel('FAST')
ylabel('SLOW')
line(xlim,xlim)
title('FIB Day 3')

nexttile;
gscatter(F.Day4(1:4),F.Day4(5:8),F.Drug(1:4),'kbgr',[],30);
xlabel('FAST')
ylabel('SLOW')
line(xlim,xlim)
title('FIB Day 4')


nexttile;
gscatter(F.Day5(1:4),F.Day5(5:8),F.Drug(1:4),'kbgr',[],30);
xlabel('FAST')
ylabel('SLOW')
line(xlim,xlim)
title('FIB Day 5')

%% plot average ESC
clrs = 'kbgr';
figure; 

tiledlayout(2,2)

nexttile; hold on ;
gscatter(E.O_N(1:4),E.O_N(5:8),E.Drug(1:4),clrs,[],30);
for I = 1:numel(clrs)
    ph1 = plot(E1.O_N(I),E1.O_N(I+4),'+','Color',clrs(I));
    ph2 = plot(E2.O_N(I),E2.O_N(I+4),'+','Color',clrs(I));
    set(ph1,'HandleVisibility','off')
    set(ph2,'HandleVisibility','off')
end
xlabel('FAST')
ylabel('SLOW')
xlim([40 100]);ylim(xlim);
line(xlim,xlim,'HandleVisibility','off');
title('ESC O/N')

nexttile;hold on ;
gscatter(E.Day2(1:4),E.Day2(5:8),E.Drug(1:4),clrs,[],30);
for I = 1:numel(clrs)
    ph1 = plot(E1.Day2(I),E1.Day2(I+4),'+','Color',clrs(I));
    ph2 = plot(E2.Day2(I),E2.Day2(I+4),'+','Color',clrs(I));
    set(ph1,'HandleVisibility','off')
    set(ph2,'HandleVisibility','off')
end
xlabel('FAST')
ylabel('SLOW')
xlim([40 100]);ylim(xlim);
line(xlim,xlim,'HandleVisibility','off');
title('ESC Day 2')
legend('location','sw')

nexttile;hold on ;
gscatter(E.Day3(1:4),E.Day3(5:8),E.Drug(1:4),clrs,[],30);
for I = 1:numel(clrs)
    ph1 = plot(E1.Day3(I),E1.Day3(I+4),'+','Color',clrs(I));
    ph2 = plot(E2.Day3(I),E2.Day3(I+4),'+','Color',clrs(I));
    set(ph1,'HandleVisibility','off')
    set(ph2,'HandleVisibility','off')
end
xlabel('FAST')
ylabel('SLOW')
xlim([40 100]);ylim(xlim);
line(xlim,xlim,'HandleVisibility','off');
title('ESC Day 3')


%% plot average FIB
clrs = 'kbgr';
figure; 

tiledlayout(2,3)

nexttile; hold on ;
gscatter(F.O_N(1:4),F.O_N(5:8),F.Drug(1:4),clrs,[],30);
for I = 1:numel(clrs)
    ph1 = plot(F1.O_N(I),F1.O_N(I+4),'+','Color',clrs(I));
    ph2 = plot(F2.O_N(I),F2.O_N(I+4),'+','Color',clrs(I));
    set(ph1,'HandleVisibility','off')
    set(ph2,'HandleVisibility','off')
end
xlabel('FAST')
ylabel('SLOW')
xlim([40 100]);ylim(xlim);
line(xlim,xlim,'HandleVisibility','off');
title('FIB O/N')


nexttile;hold on ;
gscatter(F.Day2(1:4),F.Day2(5:8),F.Drug(1:4),clrs,[],30);
for I = 1:numel(clrs)
    ph1 = plot(F1.Day2(I),F1.Day2(I+4),'+','Color',clrs(I));
    ph2 = plot(F2.Day2(I),F2.Day2(I+4),'+','Color',clrs(I));
    set(ph1,'HandleVisibility','off')
    set(ph2,'HandleVisibility','off')
end
xlabel('FAST')
ylabel('SLOW')
xlim([40 100]);ylim(xlim);
line(xlim,xlim,'HandleVisibility','off');
title('FIB Day 2')

nexttile;hold on ;
gscatter(F.Day3(1:4),F.Day3(5:8),F.Drug(1:4),clrs,[],30);
for I = 1:numel(clrs)
    ph1 = plot(F1.Day3(I),F1.Day3(I+4),'+','Color',clrs(I));
    ph2 = plot(F2.Day3(I),F2.Day3(I+4),'+','Color',clrs(I));
    set(ph1,'HandleVisibility','off')
    set(ph2,'HandleVisibility','off')
end
xlabel('FAST')
ylabel('SLOW')
xlim([40 100]);ylim(xlim);
line(xlim,xlim,'HandleVisibility','off');
title('FIB Day 3')



nexttile;hold on ;
gscatter(F.Day4(1:4),F.Day4(5:8),F.Drug(1:4),clrs,[],30);
for I = 1:numel(clrs)
    ph1 = plot(F1.Day4(I),F1.Day4(I+4),'+','Color',clrs(I));
    ph2 = plot(F2.Day4(I),F2.Day4(I+4),'+','Color',clrs(I));
    set(ph1,'HandleVisibility','off')
    set(ph2,'HandleVisibility','off')
end
xlabel('FAST')
ylabel('SLOW')
xlim([40 100]);ylim(xlim);
line(xlim,xlim,'HandleVisibility','off');
title('FIB Day 4')


nexttile;hold on ;
gscatter(F.Day5(1:4),F.Day5(5:8),F.Drug(1:4),clrs,[],30);
for I = 1:numel(clrs)
    ph1 = plot(F1.Day5(I),F1.Day5(I+4),'+','Color',clrs(I));
    ph2 = plot(F2.Day5(I),F2.Day5(I+4),'+','Color',clrs(I));
    set(ph1,'HandleVisibility','off')
    set(ph2,'HandleVisibility','off')
end
xlabel('FAST')
ylabel('SLOW')
xlim([40 100]);ylim(xlim);
line(xlim,xlim,'HandleVisibility','off');
title('FIB Day 5')