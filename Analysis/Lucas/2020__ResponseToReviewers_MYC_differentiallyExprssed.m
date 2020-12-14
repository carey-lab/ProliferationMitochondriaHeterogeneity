% 4) The authors report that Myc target genes are more highly expressed in fast proliferating cells, but 
%   is Myc itself differentially expressed? If not, how can this observation then be explained?
% Myc is more highly expressed in fast proliferating ESCs and Fibroblasts across all replicates. 
%
% LBC Response to reviewers

figname = '~/Downloads/MYC_more_highly_expressed_fast_cells__response_to_reviewers' ; 

BASEDIR = '~/Develop/ProliferationMitochondriaHeterogeneity/' ;

FIB16 = 'GeneBehaviorExampleFigure2016/FIBBehavGroups.tsv' ;
ESC16 = 'GeneBehaviorExampleFigure2016/ESCBehavGroups.tsv' ;
FIB18 = 'GeneBehaviorExampleFigure2018/FIB_CFSESep2018.tsv';
ESC18 = 'GeneBehaviorExampleFigure2018/ES_CFSESep2018.tsv' ;


FIB16 = readtable( [ BASEDIR FIB16] , 'FileType','text');
FIB18 = readtable( [ BASEDIR FIB18] , 'FileType','text');
ESC16 = readtable( [ BASEDIR ESC16] , 'FileType','text');
ESC18 = readtable( [ BASEDIR ESC18] , 'FileType','text');

FIB16.Year = repmat( 2016 , height(FIB16) , 1);
ESC16.Year = repmat( 2016 , height(ESC16) , 1);
FIB18.Year = repmat( 2018 , height(FIB18) , 1);
ESC18.Year = repmat( 2018 , height(ESC18) , 1);
FIB16.CellType = categorical( repmat( "FIB" , height(FIB16) , 1) );
ESC16.CellType = categorical( repmat( "ESC" , height(ESC16) , 1) );
FIB18.CellType = categorical( repmat( "FIB" , height(FIB18) , 1) );
ESC18.CellType = categorical( repmat( "ESC" , height(ESC18) , 1) );


GENE_ID = 'ENSMUSG00000022346' ; % MYC
%GENE_ID = 'ENSMUSG00000000303' ; % E_cadherin
%GENE_ID = 'ENSMUSG00000024304' ; % N_cadherin

% rename 18 variables
ESC18.S1_TPM = ESC18.HH_median ; % Slow replicate 1
ESC18.F1_TPM = ESC18.HL_median ; % Fast replicate 1
ESC18.S2_TPM = ESC18.MH_median ; % Slow replicate 2
ESC18.F2_TPM = ESC18.ML_median ; % Fast replicate 2
FIB18.S1_TPM = FIB18.HH_median ; % Slow replicate 1
FIB18.F1_TPM = FIB18.HL_median ; % Fast replicate 1
FIB18.S2_TPM = FIB18.MH_median ; % Slow replicate 2
FIB18.F2_TPM = FIB18.ML_median ; % Fast replicate 2
%
vn16 = ESC16.Properties.VariableNames;
vn18 = ESC18.Properties.VariableNames;
vn16idx = regexpcmp(vn16,'[SF]._TPM$') | strcmp(vn16,'Year') | strcmp(vn16,'CellType') ; 
vn18idx = regexpcmp(vn18,'[SF]._TPM$') | strcmp(vn18,'Year') | strcmp(vn18,'CellType') ; 


EXPR = vertcat( ESC16(regexpcmp(ESC16.Var1,GENE_ID), vn16idx ) , FIB16(regexpcmp(FIB16.Var1,GENE_ID), vn16idx ) ...
    , ESC18(regexpcmp(ESC18.Var1,GENE_ID), vn18idx ) , FIB18(regexpcmp(FIB18.Var1,GENE_ID), vn18idx )) ;

EXPR.log2Ratio_FoverS_1 = log2( EXPR.F1_TPM ./ EXPR.S1_TPM ) ; 
EXPR.log2Ratio_FoverS_2 = log2( EXPR.F2_TPM ./ EXPR.S2_TPM ) ; 


EXPR = stack( EXPR , {'log2Ratio_FoverS_1' 'log2Ratio_FoverS_2'} , 'NewDataVariableName','log2Ratio_FoverS')

%%
figure( 'Position', [99 99 150 200]); 
hold on ;
bh = boxplot(EXPR.log2Ratio_FoverS , {EXPR.CellType})
ylim( [-0.19 max(ylim)])
line(xlim,[0 0],'LineStyle','--','Color',[0.5 0.5 0.5])
ylabel('log_2( FAST / SLOW )     (TPM)')
plot( 1 , EXPR.log2Ratio_FoverS( EXPR.CellType == 'ESC' )  ,'.k' ,'MarkerFaceColor','k')
plot( 2 , EXPR.log2Ratio_FoverS( EXPR.CellType == 'FIB' )  ,'.k','MarkerFaceColor','k')
set(bh(1,1),'Color',[.7 .7 .7])
set(bh(1,2),'Color',[.7 .7 .7])
set(bh(2,1),'Color',[.7 .7 .7])
set(bh(2,2),'Color',[.7 .7 .7])
set(bh(3,1),'Color',[.7 .7 .7])
set(bh(3,2),'Color',[.7 .7 .7])
set(bh(4,1),'Color',[.7 .7 .7])
set(bh(4,2),'Color',[.7 .7 .7])
set(bh(5,1),'Color',[.7 .7 .7])
set(bh(5,2),'Color',[.7 .7 .7])
title('MYC expression')
print( '-dpng' , figname , '-r600') 

% %%
% ES = [ [ ES0.S1_TPM(regexpcmp(ES0.Var1,MYC))  ES0.S2_TPM(regexpcmp(ES0.Var1,MYC)) ] , ...
%        [ ES0.F1_TPM(regexpcmp(ES0.Var1,MYC))  ES0.F2_TPM(regexpcmp(ES0.Var1,MYC)) ] ]
% 
% 
% FIB = [  [ FIB0.S1_TPM(regexpcmp(FIB0.Var1,MYC))  FIB0.S2_TPM(regexpcmp(FIB0.Var1,MYC)) ] , ...
%          [ FIB0.F1_TPM(regexpcmp(FIB0.Var1,MYC))  FIB0.F2_TPM(regexpcmp(FIB0.Var1,MYC)) ] ]
% 
% d = vertcat( [log2(ES(3)/ES(1)) log2(ES(4)/ES(2))]  ,  [log2(FIB(3)/FIB(1)) log2(FIB(4)/FIB(2))] ) ; 
% 
% 
% set(gca,'ColorOrder',[0 0 0 ; 0.4 0.4 0.4])
% bar( mean(d') );
% set(gca,'xticklabel',{'FIB' 'ESC'});
% ylabel('log_2( FAST / SLOW ) (TPM)')
% %legend({'rep 1' 'rep 2'});
% set(gca,'ColorOrder',[0 0 0 ; 0.4 0.4 0.4])
% xlim([0.45 2.55])
% 
% %% showing reps
% E_cadherin = 'ENSMUSG00000000303' ;
% N_cadherin = 'ENSMUSG00000024304' ; 
% gene_names = {'E-cadherin' 'N-cadherin'}
% genes = {E_cadherin , N_cadherin} ; 
% 
% for geneI = 1:numel(genes)
%     gene = genes{geneI};
%     
% ES0(regexpcmp(ES0.Var1,gene), regexpcmp(vn0,'_TPM$') )
% 
% FIB0(regexpcmp(FIB0.Var1,gene), regexpcmp(vn0,'_TPM$') )
% 
% ES = [ [ ES0.S1_TPM(regexpcmp(ES0.Var1,gene))  ES0.S2_TPM(regexpcmp(ES0.Var1,gene)) ] , ...
%        [ ES0.F1_TPM(regexpcmp(ES0.Var1,gene))  ES0.F2_TPM(regexpcmp(ES0.Var1,gene)) ] ]
% 
% 
% FIB = [  [ FIB0.S1_TPM(regexpcmp(FIB0.Var1,gene))  FIB0.S2_TPM(regexpcmp(FIB0.Var1,gene)) ] , ...
%          [ FIB0.F1_TPM(regexpcmp(FIB0.Var1,gene))  FIB0.F2_TPM(regexpcmp(FIB0.Var1,gene)) ] ]
% 
% d = vertcat( [log2(ES(3)/ES(1)) log2(ES(4)/ES(2))]  ,  [log2(FIB(3)/FIB(1)) log2(FIB(4)/FIB(2))] ) ; 
% 
% 
% figure; 
% set(gca,'ColorOrder',[0 0 0 ; 0.4 0.4 0.4])
% bar( d );
% set(gca,'xticklabel',{'FIB' 'ESC'});
% ylabel('log_2( FAST / SLOW ) (TPM)')
% legend({'rep 1' 'rep 2'});
% set(gca,'ColorOrder',[0 0 0 ; 0.4 0.4 0.4])
% xlim([0.45 2.55])
% title( gene_names{geneI} )
% ylim([-1 0.5])
% end
