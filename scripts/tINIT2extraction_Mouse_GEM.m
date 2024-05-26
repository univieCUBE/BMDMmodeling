cd '/Users/rokosango/PhD/MetabModelling/Human-GEM/data/metabolicTasks'

essentialTasks = parseTaskList('metabolicTasks_Essential.txt')

cd /Users/rokosango/Desktop/PhD/RNA-seq/quants/

expression_data = readtable('TPMS.Final.filtered.median.csv') %then delete the first column


data_struct.tissues = expression_data.Properties.VariableNames(2:end);  % sample (tissue) names
data_struct.genes = expression_data.genes;  % gene names
data_struct.levels = table2array(expression_data(:, 2:end)); 
data_struct.threshold = 1;

cd '/Users/rokosango/MetabModelling/cobratoolbox/test/models/mat'
load('Mouse-GEM.mat')
MouseGEM = addBoundaryMets(mouseGEM) % model has to be 'closed' so that tINIT works
%% turn RAVEN like MOUSE GEM into COBRA like MOUSE GEM
cobraMouseModel = ravenCobraWrapper(MouseGEM)
writeCbModel(cobraMouseModel, 'fileName', 'MouseGEM', 'format', 'xlsx')


%It is important to first verify that the reference model (Mouse-GEM) can successfully perform all of the tasks. 
% If the reference model cannot perform a task, then neither can any GEM extracted from that model.

checkTasks(MouseGEM, [], true, false, false, essentialTasks);
% close the  human biomass reaction (MAR13082) because this messes with
% MAR0021 biomass rxn and is unable to carry a flux through this rxn
MouseGEM.rxns(MouseGEM.c == 1) % MAR00021
MouseGEM = setParam(MouseGEM, 'eq', 'MAR13082', 0)

%% use median values for expression data

refModel = MouseGEM;  % the reference model from which the GEM will be extracted
tissue = 'Control';  %IL4,Control or LPS must match the tissue name in data_struct.tissues
celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)
arrayData = data_struct;  % data structure with gene (RNA) abundance information
metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];  % we already loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
params = [];  % additional optimization parameters for the INIT algorithm
paramsFT = [];

IL4_GEM_Median = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

LPS_GEM_Median = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

Control_GEM_Median = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

setRavenSolver('gurobi')
checkTasks(IL4_GEM_Median, [], true, false, false, essentialTasks);
checkTasks(LPS_GEM_Median, [], true, false, false, essentialTasks);
checkTasks(Control_GEM_Median, [], true, false, false, essentialTasks);

%% IL4 GEM - without human biomass reaction 
%Final model statistics:
	%6418 reactions, 2709 genes
	%Mean reaction score: 6.2349
	%Mean gene score: 8.3473
	%Reactions with positive scores: 77.9371%

    IL4_GEM_Median.rxns(IL4_GEM_Median.c == 1) % MAR00021

    solveLP(IL4_GEM_Median)
    Simplified_IL4_GEM_Median = simplifyModel(IL4_GEM_Median)
    solveLP(Simplified_IL4_GEM_Median)

 %% LPS GEM - without human biomass reaction
%Final model statistics:
	%6850 reactions, 2761 genes
	%Mean reaction score: 5.2959
	%Mean gene score: 7.8082
	%Reactions with positive scores: 72.9051%

    LPS_GEM_Median.rxns(LPS_GEM_Median.c == 1) % MAR00021

    solveLP(LPS_GEM_Median)
    Simplified_LPS_GEM_Median = simplifyModel(LPS_GEM_Median)
    solveLP(Simplified_LPS_GEM_Median)

%% Control GEM without human biomass rxn
%Final model statistics:
	%6391 reactions, 2723 genes
	%Mean reaction score: 6.4189
	%Mean gene score: 8.4978
	%Reactions with positive scores: 78.2037%

    solveLP(Control_GEM_Median)
    Simplified_Control_GEM_Median = simplifyModel(Control_GEM_Median)

    constructEquations(Control_GEM_Median, 'MAR00021')
    printRxnFormula(Control_GEM_Median, 'MAR00021') % COBRA function works with RAVEN model

%% Constraining models

solveLP(LPSModel)
solveLP(Simplified_LPS_GEM_Median)

%LPS

[LPSexchangeRxns,LPSexchangeRxnsIndexes]=getExchangeRxns(LPSModel,'both')


iterModel = LPSModel

Z = cell(length(LPSexchangeRxns),1);

for i= 1:length(LPSexchangeRxns)
    iterModel = setParam(iterModel, 'lb', LPSexchangeRxns(i), -1)
    fba = solveLP(iterModel)
        if abs(fba.f) > 0.1
          Z(i) = LPSexchangeRxns(i)
        end
end

writecell(Z, "SinkReactionsLPS.csv")

LPSFluxTable = readtable('SinkReactionsLPS.csv')

LPS_rxnNameList = cellstr(LPSFluxTable.ID)
LPS_consume_values = LPSFluxTable.LB
LPS_produce_values =  LPSFluxTable.UB

LPSModel = setParam(Simplified_LPS_GEM_Median, 'ub', LPS_rxnNameList, LPS_produce_values)
LPSModel = setParam(Simplified_LPS_GEM_Median, 'lb', LPS_rxnNameList, LPS_consume_values)

printConstraints(LPSModel, -1000, 1000, LPS_rxnNameList)

LPS_fba = solveLP(LPSModel,1) %-0.6607
nnz(LPS_fba.x) % 1319


% IL4

[IL4exchangeRxns,IL4exchangeRxnsIndexes]=getExchangeRxns(IL4Model,'both')


iterModel = IL4Model

Z = cell(length(IL4exchangeRxns),1);

for i= 1:length(IL4exchangeRxns)
    iterModel = setParam(iterModel, 'lb', IL4exchangeRxns(i), -1)
    fba = solveLP(iterModel)
        if abs(fba.f) > 0.1
          Z(i) = IL4exchangeRxns(i)
        end
end

writecell(Z, "SinkReactionsIL4.csv")

IL4FluxTable = readtable('SinkReactionsIL4.csv')

IL4_rxnNameList = cellstr(IL4FluxTable.ID)
IL4_consume_values = IL4FluxTable.LB
IL4_produce_values =  IL4FluxTable.UB

IL4Model = setParam(Simplified_IL4_GEM_Median, 'ub', IL4_rxnNameList, IL4_produce_values)
IL4Model = setParam(Simplified_IL4_GEM_Median, 'lb', IL4_rxnNameList, IL4_consume_values)

printConstraints(IL4Model, -1000, 1000, IL4_rxnNameList)
printConstraints(IL4Model, -1000, 1000, 'MAR09034')
h

IL4_fba = solveLP(IL4Model,1) % f: -0.0268
nnz(IL4_fba.x) %873





%  Control

[ControlexchangeRxns,ControlexchangeRxnsIndexes]=getExchangeRxns(ControlModel,'both')


iterModel = ControlModel

Z = cell(length(ControlexchangeRxns),1);

for i= 1:length(ControlexchangeRxns)
    iterModel = setParam(iterModel, 'lb', ControlexchangeRxns(i), -1)
    fba = solveLP(iterModel)
        if abs(fba.f) > 0.1
          Z(i) = ControlexchangeRxns(i)
        end
end


writecell(Z, "SinkReactionsControl.csv")

ControlFluxTable = readtable('SinkReactionsControl.csv')

Control_rxnNameList = cellstr(ControlFluxTable.ID)
Control_consume_values = ControlFluxTable.LB
Control_produce_values =  ControlFluxTable.UB

ControlModel = setParam(Simplified_Control_GEM_Median, 'ub', Control_rxnNameList, Control_produce_values)
ControlModel = setParam(Simplified_Control_GEM_Median, 'lb', Control_rxnNameList, Control_consume_values)

printConstraints(ControlModel, -1000, 1000, Control_rxnNameList)
printConstraints(ControlModel, -1000, 1000, 'MAR09034')


Control_fba = solveLP(ControlModel,1) % f: -0.0267
nnz(Control_fba.x) % 763

%% Exporting to xlsx

  exportToExcelFormat(IL4Model, 'IL4Model.xlsx')
  exportToExcelFormat(LPSModel, 'LPSModel.xlsx')
  exportToExcelFormat(ControlModel, 'ControlModel.xlsx')

  %% Convert to COBRA models
%writeCbModel(CobraIL4Model, 'fileName', 'IL4_GEM_Median', 'format', 'xlsx')
CobraIL4Model = ravenCobraWrapper(IL4Model)
CobraLPSModel = ravenCobraWrapper(LPSModel)
CobraControlModel = ravenCobraWrapper(ControlModel)

solveLP(CobraIL4Model)% f value same as RAVEN model

%% fastFVA - with COBRA models - CPLEX 12.10 & MATLAB2019b needed

% IL4
changeCobraSolver ('ibm_cplex', 'all', 1);

[IL4minFlux,IL4maxFlux] = fastFVA(CobraIL4Model, 90)

FVAIL4 = table(CobraIL4Model.rxns,CobraIL4Model.subSystems, minFlux, maxFlux, 'VariableNames',{'ReactionID','Compartment', 'MinFlux','MaxFlux'})
writetable(FVAIL4, 'FVAIL4.csv')

% LPS

[LPSminFlux,LPSmaxFlux] = fastFVA(CobraLPSModel, 90)

FVALPS = table(CobraLPSModel.rxns, CobraLPSModel.subSystems, LPSminFlux, LPSmaxFlux, 'VariableNames',{'ReactionID','Compartment', 'MinFlux','MaxFlux'})
writetable(FVALPS, 'FVALPS.csv')

% Control

[ControlminFlux, ControlmaxFlux] = fastFVA(CobraControlModel, 90)

FVACONTROL = table(CobraControlModel.rxns, CobraControlModel.subSystems, ControlminFlux, ControlmaxFlux, 'VariableNames',{'ReactionID','Compartment', 'MinFlux','MaxFlux'})
writetable(FVACONTROL, 'FVACONTROL.csv')



%% this is for exporting and visualizing fluxes in R

  % LPS

      LPS_Flux = table(LPSModel.rxns, LPS_fba.x, 'VariableNames',{'ReactionID','FluxLPS'})
      writetable(LPS_Flux,'LPSBiomassFlux.xlsx')  
      % IL4
    
      IL4_Flux = table(IL4Model.rxns, IL4_fba.x, 'VariableNames',{'ReactionID','FluxIL4'})
      writetable(IL4_Flux,'IL4BiomassFlux.xlsx')  
      % Control
      
      Control_Flux = table(ControlModel.rxns, Control_fba.x, 'VariableNames',{'ReactionID','FluxControl'})
      writetable(Control_Flux,'ControlBiomassFlux.xlsx')  



%% TPM > 2 threshold

data_struct.threshold = 2;

Control_GEM_Median_2TPM= getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

%Final model statistics:
	%6375 reactions, 2722 genes
	%Mean reaction score: 6.44
	%Mean gene score: 8.4978
	%Reactions with positive scores: 78.4%


LPS_GEM_Median_2TPM = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);


%Final model statistics:
	%6361 reactions, 2704 genes
	%Mean reaction score: 5.8563
	%Mean gene score: 7.8082
	%Reactions with positive scores: 78.4939%


IL4_GEM_Median_2TPM = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);


%Final model statistics:
	%8713 reactions, 2849 genes
	%Mean reaction score: 4.0622
	%Mean gene score: 8.3473
	%Reactions with positive scores: 57.3855%
  



 %%% from excel to SBML


 SBMLFromExcel('IL4Model.xlsx', 'IL4ModelL.xml','false','false')
 SBMLFromExcel('LPSModel.xlsx', 'LPSModel.xml','false','false')
 SBMLFromExcel('ControlModel.xlsx', 'ControlModel.xml','false','false')

% saving to .mat files
 save('ControlModel.mat','ControlModel')
 save('LPSModel.mat','LPSModel')
 save('IL4Model.mat','IL4Model')





