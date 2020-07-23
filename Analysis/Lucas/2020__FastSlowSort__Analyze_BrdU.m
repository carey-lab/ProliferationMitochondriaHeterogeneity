%% load data
WORKDIR  = '~/Develop/GrowthPredictionTMRE/' ;
BRDU = [ WORKDIR 'Data/2020__FastSlowSort__Brdu_data.xlsx' ] ;
VIAB = [ WORKDIR 'Data/2020__FastSlowSort__viability_data.xlsx' ] ;

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


% relative to DMSO on that day
F.Rel2DMSO = NaN(height(F),1);
for I = 1:height(F)
    F.Rel2DMSO(I) = log2( ...
        F.CellsInS(I) / F.CellsInS( F.Day==F.Day(I) & F.Drug=='DMSO' & F.Speed==F.Speed(I)) ...
        );
    F.Rel2DMSO_1(I) = log2( ...
        F.Replicate1(I) / F.Replicate1( F.Day==F.Day(I) & F.Drug=='DMSO' & F.Speed==F.Speed(I)) ...
        );
    F.Rel2DMSO_2(I) = log2( ...
        F.Replicate2(I) / F.Replicate2( F.Day==F.Day(I) & F.Drug=='DMSO' & F.Speed==F.Speed(I)) ...
        );
            
end

F.Type=[];
%F.Replicate1=[];
%F.Replicate2=[];


E.Rel2DMSO = NaN(height(E),1);
for I = 1:height(E)
    E.Rel2DMSO(I) = log2( ...
        E.CellsInS(I) / E.CellsInS( E.Day==E.Day(I) & E.Drug=='DMSO' & E.Speed==E.Speed(I)) ...
        );
    E.Rel2DMSO_1(I) = log2( ...
        E.Replicate1(I) / E.Replicate1( E.Day==E.Day(I) & E.Drug=='DMSO' & E.Speed==E.Speed(I)) ...
        );
    E.Rel2DMSO_2(I) = log2( ...
        E.Replicate2(I) / E.Replicate2( E.Day==E.Day(I) & E.Drug=='DMSO' & E.Speed==E.Speed(I)) ...
        );
            
end

E.Type=[];
%E.Replicate1=[];
%E.Replicate2=[];

%%
clrs = 'gr';
figure; 
tiledlayout(2,5)


nexttile; hold on ;
idx=F.Drug=='Antimycin';
gscatter(F.Day(idx),F.Rel2DMSO_1(idx),F.Speed(idx),clrs,'s',5);
gscatter(F.Day(idx),F.Rel2DMSO_2(idx),F.Speed(idx),clrs,'^',5);
gscatter(F.Day(idx),F.Rel2DMSO(idx),F.Speed(idx),clrs,[],20);
title('FIB Antimycin')
ylabel('log2(%S / DMSO)')

nexttile; hold on ;
idx=F.Drug=='Oligomycin';
gscatter(F.Day(idx),F.Rel2DMSO_1(idx),F.Speed(idx),clrs,'s',5);
gscatter(F.Day(idx),F.Rel2DMSO_2(idx),F.Speed(idx),clrs,'^',5);
gscatter(F.Day(idx),F.Rel2DMSO(idx),F.Speed(idx),clrs,[],20);
title('FIB Oligomycin')
ylabel('log2(%S / DMSO)')


nexttile; hold on ;
idx=F.Drug=='DMSO';
gscatter(F.Day(idx),100*F.Replicate1(idx),F.Speed(idx),clrs,'s',5);
gscatter(F.Day(idx),100*F.Replicate2(idx),F.Speed(idx),clrs,'^',5);
gscatter(F.Day(idx),100*F.CellsInS(idx),F.Speed(idx),clrs,[],20);
legend('off')
title('FIB DMSO')
ylabel('% of cells in S')


nexttile; hold on ;
idx=F.Drug=='Antimycin';
gscatter(vertcat(F.Day(idx),F.Day(idx)),100*vertcat(F.Replicate1(idx),F.Replicate2(idx)),vertcat(F.Speed(idx),F.Speed(idx)),clrs,'o',5);
gscatter(F.Day(idx),100*F.CellsInS(idx),F.Speed(idx),clrs,[],20);
legend('off')
title('FIB Antimycin')
ylabel('% of cells in S')

nexttile; hold on ;
idx=F.Drug=='Oligomycin';
gscatter(vertcat(F.Day(idx),F.Day(idx)),100*vertcat(F.Replicate1(idx),F.Replicate2(idx)),vertcat(F.Speed(idx),F.Speed(idx)),clrs,'o',5);
gscatter(F.Day(idx),100*F.CellsInS(idx),F.Speed(idx),clrs,[],20);
legend('off')
title('FIB Oligomycin')
ylabel('% of cells in S')

nexttile; hold on ;
idx=E.Drug=='Antimycin';
gscatter(vertcat(E.Day(idx),E.Day(idx)),vertcat(E.Rel2DMSO_1(idx),E.Rel2DMSO_2(idx)),vertcat(E.Speed(idx),E.Speed(idx)),clrs,'o',5);
gscatter(E.Day(idx),E.Rel2DMSO(idx),E.Speed(idx),clrs,[],20);
title('ESC Antimycin')
ylabel('log2(%S / DMSO)')
legend('off')
xlabel('Day')

nexttile; hold on ;
idx=E.Drug=='Oligomycin';
gscatter(vertcat(E.Day(idx),E.Day(idx)),vertcat(E.Rel2DMSO_1(idx),E.Rel2DMSO_2(idx)),vertcat(E.Speed(idx),E.Speed(idx)),clrs,'o',5);
gscatter(E.Day(idx),E.Rel2DMSO(idx),E.Speed(idx),clrs,[],20);
title('ESC Oligomycin')
ylabel('log2(%S / DMSO)')
legend('off')
xlabel('Day')

nexttile; hold on ;
idx=E.Drug=='DMSO';
gscatter(E.Day(idx),100*E.Replicate1(idx),E.Speed(idx),clrs,'s',5);
gscatter(E.Day(idx),100*E.Replicate2(idx),E.Speed(idx),clrs,'^',5);
gscatter(E.Day(idx),100*E.CellsInS(idx),E.Speed(idx),clrs,[],20);
title('ESC DMSO')
ylabel('% of cells in S')
legend('off')
xlabel('Day')

nexttile; hold on ;
idx=E.Drug=='Antimycin';
gscatter(vertcat(E.Day(idx),E.Day(idx)),100*vertcat(E.Replicate1(idx),E.Replicate2(idx)),vertcat(E.Speed(idx),E.Speed(idx)),clrs,'o',5);
gscatter(E.Day(idx),100*E.CellsInS(idx),E.Speed(idx),clrs,[],20);
title('ESC Antimycin')
ylabel('% of cells in S')
legend('off')
xlabel('Day')

nexttile; hold on ;
idx=E.Drug=='Oligomycin';
gscatter(vertcat(E.Day(idx),E.Day(idx)),100*vertcat(E.Replicate1(idx),E.Replicate2(idx)),vertcat(E.Speed(idx),E.Speed(idx)),clrs,'o',5);
gscatter(E.Day(idx),100*E.CellsInS(idx),E.Speed(idx),clrs,[],20);
title('ESC Oligomycin')
ylabel('% of cells in S')
legend('off')
xlabel('Day')