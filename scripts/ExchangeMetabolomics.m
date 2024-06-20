
% constrain models with ExMet dataset and nutrient dataaset

ExMetLPS = readtable('/Users/rokosango/PhD/Metabolomics/data/ExMetForLPSGEM.xlsx')
AddRxns = readtable('/Users/rokosango/PhD/Metabolomics/data/AddRxnsLPS.xlsx')
ExMetCtrl = readtable('/Users/rokosango/PhD/Metabolomics/data/ExMetForCtrlGEM.xlsx')



%1. add rxns not already in the model
%2. then import nutrient values 
%3. perform fba
solveLP(Simplified_LPS_GEM_Median) % -6.5245
solveLP(Simplified_IL4_GEM_Median) % -2.8955
solveLP(Simplified_Control_GEM_Median) % -3.5865


%1.

Rxns_to_add = string(AddRxns.Rxn_ID)
LPSModel_v2 = addRxnsGenesMets(Simplified_LPS_GEM_Median, mouseGEM, 'MAR09861', true);
LPSModel_v2 = addRxnsGenesMets(LPSModel_v2, mouseGEM, 'MAR09415', true);
LPSModel_v2 = addRxnsGenesMets(LPSModel_v2, mouseGEM, 'MAR01947', true);
LPSModel_v2 = addRxnsGenesMets(LPSModel_v2, mouseGEM, 'MAR11349', true);
LPSModel_v2 = addRxnsGenesMets(LPSModel_v2, mouseGEM, 'MAR09868', true);
LPSModel_v2 = addRxnsGenesMets(LPSModel_v2, mouseGEM, 'MAR09025', true);
LPSModel_v2 = addRxnsGenesMets(LPSModel_v2, mouseGEM, 'MAR09460', true);
LPSModel_v2 = addRxnsGenesMets(LPSModel_v2, mouseGEM, 'MAR09449', true);
LPSModel_v2 = addRxnsGenesMets(LPSModel_v2, mouseGEM, 'MAR09224', true);
LPSModel_v2 = addRxnsGenesMets(LPSModel_v2, mouseGEM, 'MAR11350', true);
LPSModel_v2 = addRxnsGenesMets(LPSModel_v2, mouseGEM, 'MAR10439', true);
LPSModel_v2 = addRxnsGenesMets(LPSModel_v2, mouseGEM, 'MAR11419', true);
LPSModel_v2 = addRxnsGenesMets(LPSModel_v2, mouseGEM, 'MAR10264', true);
LPSModel_v2 = addRxnsGenesMets(LPSModel_v2, mouseGEM, 'MAR10440', true);
LPSModel_v2 = addRxnsGenesMets(LPSModel_v2, mouseGEM, 'MAR11428', true);
LPSModel_v2 = addRxnsGenesMets(LPSModel_v2, mouseGEM, 'MAR09013', true);
LPSModel_v2 = addRxnsGenesMets(LPSModel_v2, mouseGEM, 'MAR09451', true);

constructEquations(LPSModel_v2, 'MAR09259') % already in the model

% conclusion: get an error when trying to add all reactions to model,
% error message is that there are some reactions that are not in source
% model. if I then do it manually, per rxns, then everything works apart
% from MAR09013 and MAR09259 that are already in the model. WTF?????



%2. 

ExchangeRxnsLPS = readtable('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/ExchangeRxns_No_Sinks_LPS.xls')
ExchangeRxnsLPS.LB(642) = -0.5
ExchangeRxnsLPS.ID{642} = 'MAR09451'

LPS_rxnNameList = cellstr(ExchangeRxnsLPS.ID)
LPSModel_v2_consume_values = ExchangeRxnsLPS.LB
LPSModel_v2_produce_values =  ExchangeRxnsLPS.UB

LPSModel_v2_updated = setParam(LPSModel_v2, 'ub', LPS_rxnNameList, LPSModel_v2_produce_values)
LPSModel_v2_updated = setParam(LPSModel_v2_updated, 'lb', LPS_rxnNameList, LPSModel_v2_consume_values)

%print Constraints prints super weird bounds, eg. MAR10265 in
%LPSModel_v2_updated. -1000,1000 instead of 0,0. This is bc this reaction
%is "transport reaction" not "Exchange reaction" in the model. WTF??!!!


solveLP(LPSModel_v2_updated) % -1, 'The problem is infeasible'

%2.5 try now setting all sinks to -0.1

% the following is a sheet with nutrient information, extracellular
% metabolomics for LPS metabolites, and constrains rxns for metabolites
% that are significant in IL4 mets dataset and not in LPS. e.g., L-Arginine
ExchangeRxnsLPS = readtable('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/ExchangeRxns_Sinks_-0.1_LPS.xls')
ExchangeRxnsLPS.LB(642) = -0.5
ExchangeRxnsLPS.ID{642} = 'MAR09451'

LPS_rxnNameList = cellstr(ExchangeRxnsLPS.ID)
LPSModel_v2_consume_values = ExchangeRxnsLPS.LB
LPSModel_v2_produce_values =  ExchangeRxnsLPS.UB

LPSModel_v2_updated = setParam(LPSModel_v2, 'ub', LPS_rxnNameList, LPSModel_v2_produce_values)
LPSModel_v2_updated = setParam(LPSModel_v2_updated, 'lb', LPS_rxnNameList, LPSModel_v2_consume_values)

fba = solveLP(LPSModel_v2_updated)

ExMetLPSFluxTable = table(LPSModel_v2_updated.rxns, fba.x, 'VariableNames',{'ReactionID','Flux'})

writetable(ExMetLPSFluxTable, 'ExMetLPSFluxTable_-0.1.xlsx')


fba_parse = solveLP(LPSModel_v2_updated, 1)

ExMetLPSParsimoniousFluxTable = table(LPSModel_v2_updated.rxns, fba_parse.x, 'VariableNames',{'ReactionID','Flux'})

writetable(ExMetLPSParsimoniousFluxTable, 'ExMetLPSParsimoniousFluxTable_-0.1.xlsx')

% -0.01


ExchangeRxnsLPS = readtable('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/ExchangeRxns_Sinks_-0.01_LPS.xls')
ExchangeRxnsLPS.LB(642) = -0.5
ExchangeRxnsLPS.ID{642} = 'MAR09451'

LPS_rxnNameList = cellstr(ExchangeRxnsLPS.ID)
LPSModel_v2_consume_values = ExchangeRxnsLPS.LB
LPSModel_v2_produce_values =  ExchangeRxnsLPS.UB

LPSModel_v2_updated = setParam(LPSModel_v2, 'ub', LPS_rxnNameList, LPSModel_v2_produce_values)
LPSModel_v2_updated = setParam(LPSModel_v2_updated, 'lb', LPS_rxnNameList, LPSModel_v2_consume_values)

fba = solveLP(LPSModel_v2_updated, 1) % -0.0043

ExMetLPSParsimoniousFluxTable = table(LPSModel_v2_updated.rxns, fba.x, 'VariableNames',{'ReactionID','Flux'})

writetable(ExMetLPSParsimoniousFluxTable, 'ExMetLPSParsimoniousFluxTable_-0.01.xlsx')

% which rxns are needed as sinks
[LPSexchangeRxns,LPSexchangeRxnsIndexes]=getExchangeRxns(LPSModel_v2_updated,'both')


iterModel = LPSModel_v2_updated

Z = cell(length(LPSexchangeRxns),1);

for i= 1:length(LPSexchangeRxns)
    iterModel = setParam(iterModel, 'lb', LPSexchangeRxns(i), -1)
    fba = solveLP(iterModel)
        if abs(fba.f) > 0.1
          Z(i) = LPSexchangeRxns(i)
        end
end


writecell(Z, "SinksToCheckLPS.xlsx")

% with 178 sink reactions, 178 / 651 = 0.2734 (27%) of reactions added at low
% values (-1)

% this is still with one set of reactions at -0.01 LB, let's what happens
% when they're zero
ExchangeRxnsLPS = readtable('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/SinksToCheckLPS.xlsx')
ExchangeRxnsLPS.LB(642) = -0.5
ExchangeRxnsLPS.ID{642} = 'MAR09451'

LPS_rxnNameList = cellstr(ExchangeRxnsLPS.ID)
LPSModel_v2_consume_values = ExchangeRxnsLPS.LB
LPSModel_v2_produce_values =  ExchangeRxnsLPS.UB

LPSModel_v2_updated = setParam(LPSModel_v2, 'ub', LPS_rxnNameList, LPSModel_v2_produce_values)
LPSModel_v2_updated = setParam(LPSModel_v2_updated, 'lb', LPS_rxnNameList, LPSModel_v2_consume_values)

fba = solveLP(LPSModel_v2_updated, 1) % -0.0100

ExMetLPSParsiFluxTableWithCalcSinks = table(LPSModel_v2_updated.rxns, fba_parse.x, 'VariableNames',{'ReactionID','Flux'})

writetable(ExMetLPSParsiFluxTableWithCalcSinks, 'ExMetLPSParsiFluxTableWithCalcSinks.xlsx')

nnz(fba.x) %990 non zero rxns

% -0.01 -> 0 (LB) for those that are not in the of the for loop

ExchangeRxnsLPS = readtable('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/SinksToCheckLPS_v3.xlsx')
ExchangeRxnsLPS.LB(642) = -0.5
ExchangeRxnsLPS.ID{642} = 'MAR09451'

LPS_rxnNameList = cellstr(ExchangeRxnsLPS.ID)
LPSModel_v2_consume_values = ExchangeRxnsLPS.LB
LPSModel_v2_produce_values =  ExchangeRxnsLPS.UB

LPSModel_v2_updated = setParam(LPSModel_v2, 'ub', LPS_rxnNameList, LPSModel_v2_produce_values)
LPSModel_v2_updated = setParam(LPSModel_v2_updated, 'lb', LPS_rxnNameList, LPSModel_v2_consume_values)

fba = solveLP(LPSModel_v2_updated, 1) % error, the problem is unfeasible. also happens when I adjust LB = [-1,-2,-5,-10]
% this means that these set of reactions can't be used alone for biomass
% generation, no matter how much I increase uptake of these exchange into
% the system. also when I switch why reactions are 0 and which -1, still its unfeasible. 
% therefore, it looks like the system relies on the input from
% other reactions that are closed i.e. have zero as an LB. 

% maybe fix for loop.

iterModel2 = LPSModel_v2_updated

Z = cell(length(LPSexchangeRxns),1);

iterModel2 = setParam(iterModel2, 'ub', iterModel2.rxns, 1000)
iterModel2 = setParam(iterModel2, 'lb', iterModel2.rxns, 0)

for i= 1:length(LPSexchangeRxns)
    iterModel2 = setParam(iterModel2, 'lb', LPSexchangeRxns(i), -1)
    fba = solveLP(iterModel2)
        if abs(fba.f) > 0.1
          Z(i) = LPSexchangeRxns(i)
        end
end


% didn't work. maybe because I am changing one reaction at the time, the
% loop opens one reaction (at LB = -1), but keeps all other closed, that is
% why you get no biomass growth because different reactions have to be
% opened at the same time, by various LB values! 

% strategy then -> keep all reactions (for which I have no metabolomics or
% media concentrations) to very low (i.e. LB = -0.001 etc) and do this for
% each model and focus on the differences.


% test8 

ExchangeRxnsLPS = readtable('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/ExchangeRxns_Sinks_-0.01_LPS_test8.xls')
ExchangeRxnsLPS.LB(642) = -0.5
ExchangeRxnsLPS.ID{642} = 'MAR09451'

LPS_rxnNameList = cellstr(ExchangeRxnsLPS.ID)
LPSModel_v2_consume_values = ExchangeRxnsLPS.LB
LPSModel_v2_produce_values =  ExchangeRxnsLPS.UB

LPSModel_v2_updated = setParam(LPSModel_v2, 'ub', LPS_rxnNameList, LPSModel_v2_produce_values)
LPSModel_v2_updated = setParam(LPSModel_v2_updated, 'lb', LPS_rxnNameList, LPSModel_v2_consume_values)

fba = solveLP(LPSModel_v2_updated, 1)

nnz(fba.x) % 633


ExMetLPSParsiFluxTableCurated = table(LPSModel_v2_updated.rxns, fba.x, 'VariableNames',{'ReactionID','Flux'})

writetable(ExMetLPSParsiFluxTableCurated, 'ExMetLPSParsiFluxTableCurated.xlsx')



% optimize LPS for few different things (test 9)

% add NO exchange rxn
LPSModel_v2_updated = addRxnsGenesMets(LPSModel_v2_updated, mouseGEM, 'MAR05009', true);
%LPSModel_v2_updated.rxns(LPSModel_v2_updated.c == 1) % MAR00021
% put biomass rxn flux value as before but now as a LB
% put MAR04358 reaction (ATP production in glycolysis) with 0.1 as LB
% put  MAR04474 reaction (NADPH production in PPP) with 0.1 as LB
% put lactate at 0.1 MAR09135 exchange rxn


ExchangeRxnsLPS = readtable('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/ExchangeRxns_Sinks_-0.01_LPS_test9.xls')
ExchangeRxnsLPS.ID(654) = {'MAR09381'}
ExchangeRxnsLPS.LB(654) = 0.5;
ExchangeRxnsLPS.UB(654) = 1000


LPS_rxnNameList = cellstr(ExchangeRxnsLPS.ID)
LPSModel_v2_consume_values = ExchangeRxnsLPS.LB
LPSModel_v2_produce_values =  ExchangeRxnsLPS.UB
%LPSModel_v2_updated = addRxnsGenesMets(LPSModel_v2_updated, mouseGEM, 'MAR05009', true);
%LPSModel_v2_updated = addRxnsGenesMets(LPSModel_v2_updated, mouseGEM, 'MAR09381', true);
%LPSModel_v2_updated = addRxnsGenesMets(LPSModel_v2_updated, mouseGEM, 'MAR03829', true);



LPSModel_v2_updated = setParam(LPSModel_v2_updated, 'ub', LPS_rxnNameList, LPSModel_v2_produce_values)
LPSModel_v2_updated = setParam(LPSModel_v2_updated, 'lb', LPS_rxnNameList, LPSModel_v2_consume_values)

fba = solveLP(LPSModel_v2_updated, 1) %-5.2056e-04
nnz(fba.x) % 942
ExMetLPSParsiFluxTableCurated = table(LPSModel_v2_updated.rxns, fba.x, 'VariableNames',{'ReactionID','Flux'})
writetable(ExMetLPSParsiFluxTableCurated, 'ExMetLPSParsiFluxTableCurated_test9.xlsx')

%% fill gaps in LPSModel v2 


[newConnectedLPSv2, cannotConnectLPSv2, addedRxnsLPSv2, LPSGapFilledModel, exitFlag]=fillGaps(LPSModel_v2_updated,mouseGEM,'false','true','false')

ExchangeRxnsLPS = readtable('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/ExchangeRxns_Sinks_-0.01_LPS_test9.xls')
ExchangeRxnsLPS.ID(654) = {'MAR09381'}
ExchangeRxnsLPS.LB(654) = 0.5;
ExchangeRxnsLPS.UB(654) = 1000


LPS_rxnNameList = cellstr(ExchangeRxnsLPS.ID)
LPSModel_v2_consume_values = ExchangeRxnsLPS.LB
LPSModel_v2_produce_values =  ExchangeRxnsLPS.UB

LPSModel_v2_updated = LPSGapFilledModel


LPSModel_v2_updated = setParam(LPSModel_v2_updated, 'ub', LPS_rxnNameList, LPSModel_v2_produce_values)
LPSModel_v2_updated = setParam(LPSModel_v2_updated, 'lb', LPS_rxnNameList, LPSModel_v2_consume_values)

%LPSModel_v2_updated = addRxnsGenesMets(LPSModel_v2_updated, mouseGEM, 'MAR05009', true);
%LPSModel_v2_updated = addRxnsGenesMets(LPSModel_v2_updated, mouseGEM, 'MAR09381', true);
%LPSModel_v2_updated = addRxnsGenesMets(LPSModel_v2_updated, mouseGEM, 'MAR03829', true);

fba = solveLP(LPSModel_v2_updated, 1) % -9.5961
nnz(fba.x) %1317


ExMetLPSGapFilledParsiFluxTableCurated = table(LPSModel_v2_updated.rxns, fba.x, 'VariableNames',{'ReactionID','Flux'})

writetable(ExMetLPSGapFilledParsiFluxTableCurated, 'ExMetLPSGapFilledParsiFluxTableCurated.xlsx')


[LPSGapFilledexchangeRxns,LPSGapFilledexchangeRxnsIndexes]=getExchangeRxns(LPSModel_v2_updated,'both')

writecell(LPSGapFilledexchangeRxns, 'ExchangeRxns_Sinks_-0.01_LPSGapFilled.xlsx') % file which contains all rxns and fluxes


% 963 rxns, compared to previous 654. that is 309 more exchange rxns than
% the previous models


% 7950 total rxns in the gap filled model, compared to 6871. that's 1080
% more rxns

% next task. further constrain the flux space by closing these added
% exchanges

printConstraints(LPSModel_v2_updated, -1000, 1000)

LPS_LB_exchanges_gap_filled_model = LPSModel_v2_updated.lb(LPSGapFilledexchangeRxnsIndexes)
writematrix(LB_exchanges_gap_filled_model, 'LPS_LB_exchanges_gap_filled_model.xlsx') %lb for exchange rxns

LPS_UB_exchanges_gap_filled_model = LPSModel_v2_updated.ub(LPSGapFilledexchangeRxnsIndexes)
writematrix(UB_exchanges_gap_filled_model, 'LPS_UB_exchanges_gap_filled_model.xlsx') %ub for exchange rxns
%
% now see how it changes when I delete all reactions
% first close signature reactions and then inspect biomass

% test 10

ExchangeRxnsLPS = readtable('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/ExchangeRxns_Sinks_-0.01_LPSGapFilled.xlsx')

LPS_rxnNameList = cellstr(ExchangeRxnsLPS.ID)
LPSModel_v2_consume_values = ExchangeRxnsLPS.LB
LPSModel_v2_produce_values =  ExchangeRxnsLPS.UB

LPSModel_v2_updated = setParam(LPSModel_v2_updated, 'ub', LPS_rxnNameList, LPSModel_v2_produce_values)
LPSModel_v2_updated = setParam(LPSModel_v2_updated, 'lb', LPS_rxnNameList, LPSModel_v2_consume_values)
LPSModel_v2_updated = addRxnsGenesMets(LPSModel_v2_updated, mouseGEM, 'MAR03816', true);


LPS_fba = solveLP(LPSModel_v2_updated, 1) % -0.0013
nnz(LPS_fba.x) %1283


ExMetLPSGapFilledParsiFluxTableCurated = table(LPSModel_v2_updated.rxns, LPS_fba.x, 'VariableNames',{'ReactionID','Flux'})

writetable(ExMetLPSGapFilledParsiFluxTableCurated, 'LPSFluxTable.xlsx')



%%%% work on replicate LPS replicate models:
%1. 

[newConnectedLPSv2, cannotConnectLPSv2, addedRxnsLPSv2, LPSRep1_GEM_1TPM, exitFlag]=fillGaps(LPSRep1_GEM_1TPM,mouseGEM,'false','true','false')
[newConnectedLPSv2, cannotConnectLPSv2, addedRxnsLPSv2, LPSRep2_GEM_1TPM, exitFlag]=fillGaps(LPSRep2_GEM_1TPM,mouseGEM,'false','true','false')
[newConnectedLPSv2, cannotConnectLPSv2, addedRxnsLPSv2, LPSRep3_GEM_1TPM, exitFlag]=fillGaps(LPSRep3_GEM_1TPM,mouseGEM,'false','true','false')
[newConnectedLPSv2, cannotConnectLPSv2, addedRxnsLPSv2, LPSRep4_GEM_1TPM, exitFlag]=fillGaps(LPSRep4_GEM_1TPM,mouseGEM,'false','true','false')

replicate_models = struct()
replicate_models.LPS1 = LPSRep1_GEM_1TPM
replicate_models.LPS2 = LPSRep2_GEM_1TPM
replicate_models.LPS3 = LPSRep3_GEM_1TPM
replicate_models.LPS4 = LPSRep4_GEM_1TPM


fn = fieldnames(replicate_models);

for k=1:numel(fn)

if ismember(LPS_rxnNameList,replicate_models.(fn{k}).rxns) == 1

    replicate_models.(fn{k}) = setParam(replicate_models.(fn{k}), 'ub', LPS_rxnNameList, LPSModel_v2_produce_values)
    replicate_models.(fn{k}) = setParam(replicate_models.(fn{k}), 'lb', LPS_rxnNameList, LPSModel_v2_consume_values)

end
end




fluxes_replicate_models = struct()

fluxes_replicate_models.LPS1 = struct()
fluxes_replicate_models.LPS2 = struct()
fluxes_replicate_models.LPS3 = struct()
fluxes_replicate_models.LPS4 = struct()


for k=1:numel(fn)
fluxes_replicate_models.(fn{k}) = solveLP(replicate_models.(fn{k}), 1)
fluxes_replicate_models.(fn{k}) = table(fluxes_replicate_models.(fn{k}).x)
fluxes_replicate_models.(fn{k}).Properties.RowNames  = replicate_models.(fn{k}).rxns
fluxes_replicate_models.(fn{k}).Properties.VariableNames = "Flux"
end

writetable(fluxes_replicate_models.LPS1, "Rep1WtModel.csv", 'WriteRowNames',true)
writetable(fluxes_replicate_models.LPS2, "Rep2WtModel.csv", 'WriteRowNames',true)
writetable(fluxes_replicate_models.LPS3, "Rep3WtModel.csv", 'WriteRowNames',true)
writetable(fluxes_replicate_models.LPS4, "Rep4WtModel.csv", 'WriteRowNames',true)

fn = fieldnames(replicate_ko_models);
for k=1:numel(fn)

replicate_ko_models = struct()
replicate_ko_models.LPS1 = struct()
replicate_ko_models.LPS2 = struct()
replicate_ko_models.LPS3 = struct()
replicate_ko_models.LPS4 = struct()

fn = fieldnames(replicate_ko_models);
for k=1:numel(fn)

replicate_ko_models.(fn{k}) = ravenCobraWrapper(replicate_models.(fn{k}))

[grRatio, grRateKO, grRateWT, hasEffect, delRxns, replicate_ko_models.(fn{k})] = singleGeneDeletion(replicate_ko_models.(fn{k}), 'FBA', "Cyp27a1")
replicate_ko_models.(fn{k}) = array2table(replicate_ko_models.(fn{k}))
replicate_ko_models.(fn{k}).Properties.RowNames  = replicate_models.(fn{k}).rxns
replicate_ko_models.(fn{k}).Properties.VariableNames = "Flux"

end

writetable(replicate_ko_models.LPS1, "Rep1KoModel.csv", 'WriteRowNames',true)
writetable(replicate_ko_models.LPS2, "Rep2KoModel.csv", 'WriteRowNames',true)
writetable(replicate_ko_models.LPS3, "Rep3KoModel.csv", 'WriteRowNames',true)
writetable(replicate_ko_models.LPS4, "Rep4KoModel.csv", 'WriteRowNames',true)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% done with LPS, now do the same for IL4 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ExMetIL4 = readtable('/Users/rokosango/PhD/Metabolomics/data/ExMetForIL4GEM.xlsx')


%1. add rxns not already in the model
%2. then import nutrient values 
%3. perform fba
solveLP(Simplified_LPS_GEM_Median) % -6.5245
solveLP(Simplified_IL4_GEM_Median) % -2.8955

%1.

Rxns_to_add = string(ExMetIL4.Rxn_ID)
IL4Model_v2 = addRxnsGenesMets(Simplified_IL4_GEM_Median, mouseGEM, 'MAR11428', true);
IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR10440', true);
IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR10264', true);
IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR11419', true);
IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR10439', true);
IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR11350', true);
IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR09224', true);
IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR09449', true);
IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR09437', true);
%IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR11438', true);
%already in the model
IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR09919', true);
%IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR12200', true);
%already in the model

% LPS reactions that need to be added and then closed since they're only
% for LPS
% then overlay nutrient information on top. if some certain reactions you
% have different values for exMet and nutrient, nutrient takes precedence
IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR09460', true);
IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR09025', true);
IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR09013', true);
IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR09868', true);
IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR11349', true);
IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR01947', true);
IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR09451', true);
IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR09415', true);
IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR09861', true);
IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR08386', true);
IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR11974', true);
IL4Model_v2 = addRxnsGenesMets(IL4Model_v2, mouseGEM, 'MAR09247', true);

constructEquations(IL4Model_v2, 'MAR11438') % already in the model
constructEquations(IL4Model_v2, 'MAR12200') % already in the model

% conclusion: get an error when trying to add all reactions to model,
% error message is that there are some reactions that are not in source
% model. if I then do it manually, per rxns, then everything works apart
% from MAR09013 and MAR09259 that are already in the model. WTF?????

[IL4exchangeRxns,IL4exchangeRxnsIndexes]=getExchangeRxns(IL4Model_v2,'both')

writecell(IL4exchangeRxns, 'ExchangeRxns_Sinks_-0.01_IL4.xlsx') % file which contains all rxns and fluxes


% special test 14

IL4FluxTable = readtable('SinkReactionsIL4_test14.csv')

IL4_rxnNameList = cellstr(IL4FluxTable.ID)
IL4_consume_values = IL4FluxTable.LB
IL4_produce_values =  IL4FluxTable.UB

IL4Model_v2_updated = setParam(IL4Model_v2, 'ub', IL4_rxnNameList, IL4_produce_values)
IL4Model_v2_updated = setParam(IL4Model_v2_updated, 'lb', IL4_rxnNameList, IL4_consume_values)


fba = solveLP(IL4Model_v2_updated, 1) % -2.6300e-05


nnz(fba.x) % 790

%special test 15 - lactate exchange closed, and include urea rxns

IL4FluxTable = readtable('SinkReactionsIL4_test14.csv')

IL4_rxnNameList = cellstr(IL4FluxTable.ID)
IL4_consume_values = IL4FluxTable.LB
IL4_produce_values =  IL4FluxTable.UB

IL4Model_v2_updated = setParam(IL4Model_v2, 'ub', IL4_rxnNameList, IL4_produce_values)
IL4Model_v2_updated = setParam(IL4Model_v2_updated, 'lb', IL4_rxnNameList, IL4_consume_values)

fba = solveLP(IL4Model_v2_updated, 1) % -2.6300e-05


fba = solveLP(IL4Model_v2_updated, 1) % -2.6300e-05


nnz(fba.x) % 778

ExMetIL4ParsiFluxTableCurated = table(IL4Model_v2_updated.rxns, fba.x, 'VariableNames',{'ReactionID','Flux'})

writetable(ExMetIL4ParsiFluxTableCurated, 'ExMetIL4ParsiFluxTableCurated_15.xlsx')


%% try with gap filled IL4 model

[newConnectedIL4v2, cannotConnectIL4v2, addedRxnsIL4v2, IL4GapFilledModel, exitFlag]=fillGaps(IL4Model_v2_updated,mouseGEM,'false','true','false')

IL4FluxTable = readtable('SinkReactionsIL4_test14.csv') 

IL4_rxnNameList = cellstr(IL4FluxTable.ID)
IL4_consume_values = IL4FluxTable.LB
IL4_produce_values =  IL4FluxTable.UB

IL4Model_v2_updated = setParam(IL4GapFilledModel, 'ub', IL4_rxnNameList, IL4_produce_values)
IL4Model_v2_updated = setParam(IL4Model_v2_updated, 'lb', IL4_rxnNameList, IL4_consume_values)



%IL4Model_v2_updated = setParam(IL4Model_v2_updated, 'lb', 'MAR09068', 0.1)
fba = solveLP(IL4Model_v2_updated, 1) % -8.1271

nnz(fba.x) % 1048


[IL4GapFilledexchangeRxns,IL4GapFilledexchangeRxnsIndexes]=getExchangeRxns(IL4Model_v2_updated,'both')

writecell(IL4GapFilledexchangeRxns, 'ExchangeRxns_Sinks_-0.01_IL4GapFilled.xlsx') % file which contains all rxns and fluxes


% 924 rxns, compared to previous 538. that is 386 more exchange rxns than
% the previous models


% 7686 total rxns in the gap filled model, compared to 6441. that's 1245
% more rxns

% next task. further constrain the flux space by closing these added
% exchanges


IL4_LB_exchanges_gap_filled_model = IL4Model_v2_updated.lb(IL4GapFilledexchangeRxnsIndexes)
writematrix(IL4_LB_exchanges_gap_filled_model, 'IL4_LB_exchanges_gap_filled_model.xlsx') %lb for exchange rxns

IL4_UB_exchanges_gap_filled_model = IL4Model_v2_updated.ub(IL4GapFilledexchangeRxnsIndexes)
writematrix(IL4_UB_exchanges_gap_filled_model, 'IL4_UB_exchanges_gap_filled_model.xlsx') % ub for exchange rxns

ExchangeRxnsIL4 = readtable('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/ExchangeRxns_Sinks_-0.01_IL4GapFilled.xlsx')

IL4_rxnNameList = cellstr(ExchangeRxnsIL4.ID)
IL4Model_v2_consume_values = ExchangeRxnsIL4.LB
IL4Model_v2_produce_values =  ExchangeRxnsIL4.UB

IL4Model_v2_updated = setParam(IL4Model_v2_updated, 'ub', IL4_rxnNameList, IL4Model_v2_produce_values)
IL4Model_v2_updated = setParam(IL4Model_v2_updated, 'lb', IL4_rxnNameList, IL4Model_v2_consume_values)

IL4_fba = solveLP(IL4Model_v2_updated, 1) % -0.0018, -0.0015 for LPS
nnz(IL4_fba.x) %1103

% -9.2085e-04 is biomass when I raise urea + ornithine reaction from 0.1 to
% 0.2


ExMetIL4GapFilledParsiFluxTableCurated = table(IL4Model_v2_updated.rxns, fba.x, 'VariableNames',{'ReactionID','Flux'})

writetable(ExMetIL4GapFilledParsiFluxTableCurated, 'IL4FluxTable.xlsx')




%choose special test 14
% conclusion - can close 457 reactions and still have growth.
% 558 (total) - 457 = 101 / 558 = 18 %
% in another words, I need 18 % of total reactions being sinks at low
% values i.e. 0.01 as a LB.

%ExMetIL4ParsiFluxTableCurated = table(IL4Model_v2_updated.rxns, fba.x, 'VariableNames',{'ReactionID','Flux'})

%writetable(ExMetIL4ParsiFluxTableCurated, 'ExMetIL4ParsiFluxTableCurated_15.xlsx')


%% Control sink constraint
Rxns_to_add = string(ExMetCtrl.Rxn_ID)

%%% important note:
% if you supply a list of reactions to this function, if ANY of the
% reactions are already in the model, the function crashes. meaning, if
% within those reactions there are novel ones, they won't be added and will
% be missed by this tool. Hence, I need to manually check membership of
% each reaction. another WTF!!!!!

% LPS reactions that need to be added and then closed since they're only
% for LPS
% then overlay nutrient information on top. if some certain reactions you
% have different values for exMet and nutrient, nutrient takes precedence

% now go to Combined LPS + IL4 ExMet, look for metabolites that are
% significant for both and specific for condition, and close those rxns for control
% as I don't want to provide control model metabolites for which no
% evidence is given


[ControlExchangeRxns,ControlExchangeRxnsIndexes]=getExchangeRxns(ControlModel_v2,'both')

writecell(ControlExchangeRxns, "ControlExchangeRxns.xlsx")



% test 1 - lactate exchange closed


CtrlFluxTable = readtable('SinksReactionsControl_test1.csv')

Ctrl_rxnNameList = cellstr(CtrlFluxTable.ID)
Ctrl_consume_values = CtrlFluxTable.LB
Ctrl_produce_values =  CtrlFluxTable.UB

ControlModel_v2_updated = setParam(ControlModel_v2, 'ub', Ctrl_rxnNameList, Ctrl_produce_values)
ControlModel_v2_updated = setParam(ControlModel_v2_updated, 'lb', Ctrl_rxnNameList, Ctrl_consume_values)


setRavenSolver('gurobi')
fba = solveLP(ControlModel_v2_updated, 1) %-2.6300e-05

nnz(fba.x) % 959


%%%% conclusion: 115 / 567 = 20% reactions have to be added as sinks at low
%%%% concs (-0.01). and still have producing biomass


ExMetCtrlParsiFluxTableCurated = table(ControlModel_v2_updated.rxns, fba.x, 'VariableNames',{'ReactionID','Flux'})

writetable(ExMetCtrlParsiFluxTableCurated, 'ExMetCtrlParsiFluxTableCurated.xlsx')




%% Gap fill for Control GEM

[newConnectedControlv2, cannotConnectControlv2, addedRxnsControlv2, ControlGapFilledModel, exitFlag]=fillGaps(ControlModel_v2_updated,mouseGEM,'false','true','false')


[ControlGapFilledexchangeRxns,ControlGapFilledexchangeRxnsIndexes]=getExchangeRxns(ControlGapFilledModel,'both')

writecell(ControlGapFilledexchangeRxns, 'ExchangeRxns_Sinks_-0.01_ControlGapFilled.xlsx') % file which contains all rxns and fluxes


Control_LB_exchanges_gap_filled_model = ControlGapFilledModel.lb(ControlGapFilledexchangeRxnsIndexes)
writematrix(Control_LB_exchanges_gap_filled_model, 'Control_LB_exchanges_gap_filled_model.xlsx') %lb for exchange rxns

Control_UB_exchanges_gap_filled_model = ControlGapFilledModel.ub(ControlGapFilledexchangeRxnsIndexes)
writematrix(Control_UB_exchanges_gap_filled_model, 'Control_UB_exchanges_gap_filled_model.xlsx') % ub for exchange rxns



% 914 exchange reactions in gap filled model compared to 567 in control v2
% = 347 more exchanges

% 6405 rxns in COntrol v2 vs 7668 in Control gap filled = 1263  more rxns


CtrlFluxTable = readtable('ExchangeRxns_Sinks_-0.01_ControlGapFilled.xlsx')

Ctrl_rxnNameList = cellstr(CtrlFluxTable.ID)
Ctrl_consume_values = CtrlFluxTable.LB
Ctrl_produce_values =  CtrlFluxTable.UB

ControlGapFilledModel = setParam(ControlGapFilledModel, 'ub', Ctrl_rxnNameList, Ctrl_produce_values)
ControlGapFilledModel = setParam(ControlGapFilledModel, 'lb', Ctrl_rxnNameList, Ctrl_consume_values)

ControlGapFilledModel = addRxnsGenesMets(ControlGapFilledModel, mouseGEM, 'MAR03816', true);

Control_fba = solveLP(ControlGapFilledModel, 1) % -0.0011 for Control, -0.0013 for IL4 -0.0015 for LPS

nnz(Control_fba.x) %947


ExMetCtrlGapFilledParsiFluxTableCurated = table(ControlGapFilledModel.rxns, Control_fba.x, 'VariableNames',{'ReactionID','Flux'})

writetable(ExMetCtrlGapFilledParsiFluxTableCurated, 'ControlFluxTable.xlsx')



%check for essential genes with CheckTasksGenes

TasksPath = '/Users/rokosango/PhD/Modelling/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.xlsx'
[LPSTaskReport, LPSessentialGenes, taskStructure] = checkTasksGenes(LPSModel_v2_updated, TasksPath, [], [], true)
[IL4TaskReport, IL4essentialGenes, taskStructure] = checkTasksGenes(IL4Model_v2_updated, TasksPath, [], [], true)
[CtrlTaskReport, CtrlessentialGenes, taskStructure] = checkTasksGenes(ControlModel_v2_updated, TasksPath, [], [], true)

writematrix(LPSessentialGenes, 'LPSessentialGenes.csv')
writematrix(IL4essentialGenes, 'IL4essentialGenes.csv')
writematrix(CtrlessentialGenes, 'CtrlessentialGenes.csv')

writecell(LPSModel_v2_updated.genes, 'LPSGapFilledModelGenes.csv')
writecell(IL4Model_v2_updated.genes, 'IL4GapFilledModelGenes.csv')
writecell(ControlModel_v2_updated.genes, 'ControlGapFilledModelGenes.csv')


%check for gene essentiality with singleGeneDeletion.m 
%LPS

LPSModel_v2_updated = addRxnsGenesMets(LPSModel_v2_updated, mouseGEM, 'MAR03816', true);


CobraLPSGapFilledModel = ravenCobraWrapper(LPSModel_v2_updated)

[LPSgrRatio, LPSgrRateKO, LPSgrRateWT, LPShasEffect, LPSdelRxns, LPSfluxSolution] = singleGeneDeletion(CobraLPSGapFilledModel, ...
    'FBA', CobraLPSGapFilledModel.genes)

LPSGapFilledEssGenes = table(LPSgrRatio, LPSgrRateKO, LPShasEffect, CobraLPSGapFilledModel.genes)
writetable(LPSGapFilledEssGenes, 'LPSGapFilledEssGenes.csv')


%IL4

CobraIL4GapFilledModel = ravenCobraWrapper(IL4Model_v2_updated)

[IL4grRatio, IL4grRateKO, IL4grRateWT, IL4hasEffect, IL4delRxns, IL4fluxSolution] = singleGeneDeletion(CobraIL4GapFilledModel, ...
    'FBA', CobraIL4GapFilledModel.genes)

IL4GapFilledEssGenes = table(IL4grRatio, IL4grRateKO, IL4hasEffect, CobraIL4GapFilledModel.genes)
writetable(IL4GapFilledEssGenes, 'IL4GapFilledEssGenes.csv')


%Control
CobraCtrlGapFilledModel = ravenCobraWrapper(ControlGapFilledModel)

[CtrlgrRatio, CtrlgrRateKO, CtrlgrRateWT, CtrlhasEffect, CtrldelRxns, CtrlfluxSolution] = singleGeneDeletion(CobraCtrlGapFilledModel, ...
    'FBA', CobraCtrlGapFilledModel.genes)

CtrlGapFilledEssGenes = table(CtrlgrRatio, CtrlgrRateKO, CtrlhasEffect, CobraCtrlGapFilledModel.genes)
writetable(CtrlGapFilledEssGenes, 'CtrlGapFilledEssGenes.csv')



% fastFVA - MATLAB 2019b & CPlex optimization studio 1210 is needed ! also
% note long waiting time for the calculation. 1-2 hrs for each fastFVA run
% with individual models

[LPSminFlux,LPSmaxFlux] = fastFVA(CobraLPSGapFilledModel, 90)

[min,max] = fastFVA(model, 90)

[LPSminFlux,LPSmaxFlux] = fastFVA(CobraLPSGapFilledModel, 90)

[IL4minFlux,IL4maxFlux] = fastFVA(CobraIL4GapFilledModel, 90)

[ControlminFlux, ControlmaxFlux] = fastFVA(CobraCtrlGapFilledModel, 90)

FVALPS = table(CobraLPSGapFilledModel.rxns, CobraLPSGapFilledModel.subSystems, LPSminFlux, LPSmaxFlux, 'VariableNames',{'ReactionID','Compartment', 'MinFlux','MaxFlux'})
writetable(FVALPS, 'FVALPS.csv')

FVAIL4 = table(CobraIL4GapFilledModel.rxns,CobraIL4GapFilledModel.subSystems, IL4minFlux, IL4maxFlux, 'VariableNames',{'ReactionID','Compartment', 'MinFlux','MaxFlux'})
writetable(FVAIL4, 'FVAIL4.csv')

FVACONTROL = table(CobraCtrlGapFilledModel.rxns, CobraCtrlGapFilledModel.subSystems, ControlminFlux, ControlmaxFlux, 'VariableNames',{'ReactionID','Compartment', 'MinFlux','MaxFlux'})
writetable(FVACONTROL, 'FVACONTROL.csv')

FVALPS = table(CobraLPSGapFilledModel.rxns, CobraLPSGapFilledModel.subSystems, LPSminFlux, LPSmaxFlux, 'VariableNames',{'ReactionID','Compartment', 'MinFlux','MaxFlux'})
writetable(FVALPS, 'FVALPS.csv')




% sampling

changeCobraSolver('gurobi')

[LPSModelSampling, LPSsamples] = sampleCbModel(CobraLPSGapFilledModel, 'LPSModelFile', 'ACHR')

[IL4ModelSampling, IL4samples] = sampleCbModel(CobraIL4GapFilledModel, 'IL4ModelFile', 'ACHR')

[ControlModelSampling, Controlsamples] = sampleCbModel(CobraCtrlGapFilledModel, 'ControlModelFile', 'ACHR')

%try also with LPS GAPDH Knockout Model

[Gapdh_KO_ModelSampling, Gapdh_KO_Modelsamples] = sampleCbModel(Gapdh_KO_Model, 'Gapdh_KO_ModelFile', 'ACHR')

Gapdh_KO_Model

% try also with LPS Cyp27a1 Knockout Flux


LPS_KO_Model = CobraLPSGapFilledModel

reactions_to_knockout = CobraLPSGapFilledModel.rxns(string(CobraLPSGapFilledModel.subSystems) == 'Bile acid biosynthesis')

LPS_KO_Model = changeRxnBounds(LPS_KO_Model, reactions_to_knockout, 0, 'b')



optimizeCbModel(LPS_KO_Model)

[KO_LPSModelSampling, KO_LPSsamples] = sampleCbModel(LPS_KO_Model, 'KO_LPSModelFile', 'ACHR')

% try also with IL-4 Cyp27a1 Knockout Flux

IL4_KO_Model = CobraIL4GapFilledModel
reactions_to_knockout = CobraIL4GapFilledModel.rxns(string(CobraIL4GapFilledModel.subSystems) == 'Bile acid biosynthesis')
IL4_KO_Model = changeRxnBounds(IL4_KO_Model, reactions_to_knockout, 0, 'b')
optimizeCbModel(IL4_KO_Model)

[KO_IL4ModelSampling, KO_IL4samples] = sampleCbModel(IL4_KO_Model, 'KO_IL4ModelFile', 'ACHR')


%try also with Ctrl Cyp27a1 Knockout

Ctrl_KO_Model = CobraCtrlGapFilledModel
reactions_to_knockout = CobraCtrlGapFilledModel.rxns(string(CobraCtrlGapFilledModel.subSystems) == 'Bile acid biosynthesis')
Ctrl_KO_Model = changeRxnBounds(Ctrl_KO_Model, reactions_to_knockout, 0, 'b')

% IL-4 KO FVA

[KO_IL4minFlux,KO_IL4maxFlux] = fastFVA(IL4_KO_Model, 90)


FVA_KO_LPS = table(LPS_KO_Model.rxns, LPS_KO_Model.subSystems, KO_LPSminFlux, KO_LPSmaxFlux, 'VariableNames',{'ReactionID','Compartment', 'MinFlux','MaxFlux'})
writetable(FVA_KO_LPS, 'FVA_KO_LPS.csv')

[Gapdh_KO_minFlux, Gapdh_KO_maxFlux] = fastFVA(Gapdh_KO_Model, 90)

FVA_Gapdh_KO_LPS = table(Gapdh_KO_Model.rxns, Gapdh_KO_Model.subSystems, Gapdh_KO_minFlux, Gapdh_KO_maxFlux, 'VariableNames',{'ReactionID','Compartment', 'MinFlux','MaxFlux'})
writetable(FVA_Gapdh_KO_LPS, 'FVA_Gapdh_KO_LPS.csv')

FVA_Gapdh_KO_LPS(FVA_Gapdh_KO_LPS.ReactionID == "MAR04373",:)

%get ID of reactions from each sampled Model, and then extract as
%spreadsheet each `points` dataset


writecell(LPSModelSampling.rxns, 'LPSModelSamplingRxns')
writecell(IL4ModelSampling.rxns, 'IL4ModelSamplingRxns')
writecell(KO_LPSModelSampling.rxns, 'KO_LPSModelSamplingRxns')
writecell(KO_IL4ModelSampling.rxns, 'KO_IL4ModelSamplingRxns')


% next come 'point' - these are saved as .mat files - 10 files in total
%LPS datasets - save as xlsx one by one

csvwrite('LPSModelFile1.csv', points)
csvwrite('LPSModelFile2.csv', points)
csvwrite('LPSModelFile3.csv', points)
csvwrite('LPSModelFile4.csv', points)
csvwrite('LPSModelFile5.csv', points)
csvwrite('LPSModelFile6.csv', points)
csvwrite('LPSModelFile7.csv', points)
csvwrite('LPSModelFile8.csv', points)
csvwrite('LPSModelFile9.csv', points)
csvwrite('LPSModelFile10.csv', points)

% next come 'point' - these are saved as .mat files - 10 files in total
%IL4 datasets - save as xlsx one by one

csvwrite('IL4ModelFile1.csv', points)
csvwrite('IL4ModelFile2.csv', points)
csvwrite('IL4ModelFile3.csv', points)
csvwrite('IL4ModelFile4.csv', points)
csvwrite('IL4ModelFile5.csv', points)
csvwrite('IL4ModelFile6.csv', points)
csvwrite('IL4ModelFile7.csv', points)
csvwrite('IL4ModelFile8.csv', points)
csvwrite('IL4ModelFile9.csv', points)
csvwrite('IL4ModelFile10.csv', points)


% next come 'point' - these are saved as .mat files - 10 files in total
%LPS KO datasets - save as xlsx one by one
%example: load('KO_LPSModelFile_1.mat')

csvwrite('KO_LPSModelFile1.csv', points)
csvwrite('KO_LPSModelFile2.csv', points)
csvwrite('KO_LPSModelFile3.csv', points)
csvwrite('KO_LPSModelFile4.csv', points)
csvwrite('KO_LPSModelFile5.csv', points)
csvwrite('KO_LPSModelFile6.csv', points)
csvwrite('KO_LPSModelFile7.csv', points)
csvwrite('KO_LPSModelFile8.csv', points)
csvwrite('KO_LPSModelFile9.csv', points)
csvwrite('KO_LPSModelFile10.csv', points)


% next come 'point' - these are saved as .mat files - 10 files in total
%LPS KO datasets - save as xlsx one by one
%example: load('KO_LPSModelFile_1.mat')

csvwrite('KO_IL4ModelFile1.csv', points)
csvwrite('KO_IL4ModelFile2.csv', points)
csvwrite('KO_IL4ModelFile3.csv', points)
csvwrite('KO_IL4ModelFile4.csv', points)
csvwrite('KO_IL4ModelFile5.csv', points)
csvwrite('KO_IL4ModelFile6.csv', points)
csvwrite('KO_IL4ModelFile7.csv', points)
csvwrite('KO_IL4ModelFile8.csv', points)
csvwrite('KO_IL4ModelFile9.csv', points)
csvwrite('KO_IL4ModelFile10.csv', points)



% Example for the standard Gaussian distribution
%numVar = size(CobraLPSGapFilledModel.S, 2);
%options.nWorkers = 1
%CobraLPSGapFilledModel.vMean = zeros(numVar, 1); 
%CobraLPSGapFilledModel.vCov = ones(numVar, 1);
%[LPSModelSamplingRHMC, LPSsamplesRHMC] = sampleCbModel(CobraLPSGapFilledModel, 'LPSModelFileRHMC', 'RHMC', options);



% get flux distribution of the knock out models LPS/IL4/Control

geneList = ["Cyp27a1" "Phgdh" "Alox15" "Setdb2", "Sphk1", "Sphk2", "Acly"]

%Gene Deletion with pFBA instead of FBA - inserted a function in
%singleGeneDeletion.m

Gapdh_KO_Model = setParam(CobraLPSGapFilledModel, 'eq', 'MAR04373', 0) 
%"Gapdhs" "Gapdh-ps15" gene knockouts don't work for some reason

Gapdh_fba = solveLP(Gapdh_KO_Model, 1)


Arg1_KO_Model = setParam(CobraLPSGapFilledModel, 'eq', ['MAR08426', 'MAR03816'], 0)

Arg1_fba = solveLP(Arg1_KO_Model, 1)


[grRatio, grRateKO, grRateWT, hasEffect, delRxns, LPSGeneDelFluxSolution] = singleGeneDeletion(CobraLPSGapFilledModel, 'FBA', geneList)
LPSGeneDelFluxSolution = array2table(LPSGeneDelFluxSolution)

A = LPSGeneDelFluxSolution
B = table(Gapdh_fba.x, Arg1_fba.x)

LPSGeneDelFluxSolution = [A B]
LPSGeneDelFluxSolution.Properties.VariableNames = [geneList, "Gapdh" "Arg1"]
LPSGeneDelFluxSolution.Properties.RowNames = CobraLPSGapFilledModel.rxns


%same genelist in IL4 setting


Gapdh_KO_Model = setParam(CobraIL4GapFilledModel, 'eq', 'MAR04373', 0) %"Gapdhs" "Gapdh-ps15" gene knockouts don't work for some reason

Gapdh_fba = solveLP(Gapdh_KO_Model, 1)

Arg1_KO_Model = setParam(CobraIL4GapFilledModel, 'eq', ['MAR08426', 'MAR03816'], 0)

Arg1_fba = solveLP(Arg1_KO_Model, 1)

[grRatio, grRateKO, grRateWT, hasEffect, delRxns, IL4GeneDelFluxSolution] = singleGeneDeletion(CobraIL4GapFilledModel, 'FBA', geneList)
IL4GeneDelFluxSolution = array2table(IL4GeneDelFluxSolution)


A = IL4GeneDelFluxSolution
B = table(Gapdh_fba.x, Arg1_fba.x)
IL4GeneDelFluxSolution = [A B]
IL4GeneDelFluxSolution.Properties.VariableNames = [geneList, "Gapdh" "Arg1"]
IL4GeneDelFluxSolution.Properties.RowNames = CobraIL4GapFilledModel.rxns


%same genelist in Control setting except Arg1

%add reaction MAR03816 to the Control Model (ARGN1)

%CobraCtrlGapFilledModel = ravenCobraWrapper(ControlGapFilledModel)

%ControlGapFilledModel = addRxnsGenesMets(ControlGapFilledModel, mouseGEM, 'MAR03816', true);

%CobraCtrlGapFilledModel = ravenCobraWrapper(ControlGapFilledModel)


Gapdh_KO_Model = setParam(CobraCtrlGapFilledModel, 'eq', 'MAR04373', 0) %"Gapdhs" "Gapdh-ps15" gene knockouts don't work for some reason

Gapdh_fba = solveLP(Gapdh_KO_Model, 1)

Arg1_KO_Model = setParam(CobraCtrlGapFilledModel, 'eq', ['MAR08426', 'MAR03816'], 0)

Arg1_fba = solveLP(Arg1_KO_Model, 1)


[grRatio, grRateKO, grRateWT, hasEffect, delRxns, CtrlGeneDelFluxSolution] = singleGeneDeletion(CobraCtrlGapFilledModel, 'FBA', geneList)
CtrlGeneDelFluxSolution = array2table(CtrlGeneDelFluxSolution)

A = CtrlGeneDelFluxSolution
B = table(Gapdh_fba.x, Arg1_fba.x)
CtrlGeneDelFluxSolution = [A B]

CtrlGeneDelFluxSolution.Properties.VariableNames = [geneList, "Gapdh" "Arg1"]
CtrlGeneDelFluxSolution.Properties.RowNames = CobraCtrlGapFilledModel.rxns


writetable(LPSGeneDelFluxSolution, "LPSGeneDelFluxSolution.xlsx", 'WriteRowNames',true)
writetable(IL4GeneDelFluxSolution, "IL4GeneDelFluxSolution.xlsx", 'WriteRowNames',true)
writetable(CtrlGeneDelFluxSolution, "CtrlGeneDelFluxSolution.xlsx", 'WriteRowNames',true)



% checking for futile cycles in LPS WT vs GAPDH mutant

% * alternatively, you fix growth rate to the optimum and 
% PGM to its maximum (1000) and then do an FVA on all other 
% reactions (or just the one that you are interested in). 
% You will find that some reactions that were previously 
% not fixed (when PGM was variable) are now fixed too. Those
% reactions that became fixed, are part of the cycle. Repeat with PGM fixed to -1000

Gapdh_fba = solveLP(Gapdh_KO_Model, 1) %-0.0011
solveLP(CobraLPSGapFilledModel, 1) %-0.0013

Gapdh_KO_Model.rxns(Gapdh_KO_Model.c == 1) %MAR00021

Gapdh_KO_Model_futile_cycle = setParam(Gapdh_KO_Model, 'lb', 'MAR00021', 0.0011)
Gapdh_KO_Model_futile_cycle = setParam(Gapdh_KO_Model_futile_cycle, 'lb', 'MAR04365', 1000)

WT_Model_futile_cycle = setParam(CobraLPSGapFilledModel, 'lb', 'MAR00021', 0.0013)
WT_Model_futile_cycle = setParam(WT_Model_futile_cycle, 'lb', 'MAR04365', 1000)



[Gapdh_KO_Model_futile_cycle_minFlux,Gapdh_KO_Model_futile_cycle_maxFlux] = fastFVA(Gapdh_KO_Model_futile_cycle, 90)

[WT_Model_futile_cycle_minFlux,WT_Model_futile_cycle_maxFlux] = fastFVA(WT_Model_futile_cycle, 90)

FVA_Gapdh_KO_LPS_futile_cycle = table(Gapdh_KO_Model_futile_cycle.rxns, ...
                        Gapdh_KO_Model_futile_cycle.subSystems, ...
                        Gapdh_KO_Model_futile_cycle_minFlux, ...
                        Gapdh_KO_Model_futile_cycle_maxFlux, ...
                        'VariableNames',{'ReactionID','Compartment', 'MinFlux','MaxFlux'})
                    
writetable(FVA_Gapdh_KO_LPS_futile_cycle, 'FVA_Gapdh_KO_LPS_futile_cycle.csv')

FVA_WT_Model_futile_cycle = table(WT_Model_futile_cycle.rxns, ...
                        WT_Model_futile_cycle.subSystems, ...
                        WT_Model_futile_cycle_minFlux, ...
                        WT_Model_futile_cycle_maxFlux, ...
                        'VariableNames',{'ReactionID','Compartment', 'MinFlux','MaxFlux'})
       
writetable(FVA_WT_Model_futile_cycle, 'FVA_WT_Model_futile_cycle.csv')




% you remove your carbon source and (say in the mutant) maximise for 
% PGM—the resulting flux distribution will be a cycle.


Gapdh_KO_Model_futile_cycle = setParam(Gapdh_KO_Model, 'eq', 'MAR09034', 0)
Gapdh_KO_Model_futile_cycle = setParam(Gapdh_KO_Model_futile_cycle, 'obj', 'MAR04365', 1)

fba_futile_gapdh_ko = solveLP(Gapdh_KO_Model_futile_cycle, 1)

FluxTable_futile_gapdh = table(Gapdh_KO_Model_futile_cycle.rxns, ...
                        fba_futile_gapdh_ko.x, ...
                        'VariableNames',{'ReactionID', 'Flux'})


writetable(FluxTable_futile_gapdh, 'FluxTable_futile_gapdh.xlsx')


%%%% write Cobra models to mat

save('IL4Model.mat', "IL4Model")
save('ControlModel.mat', "ControlModel")
save('LPSModel.mat', "LPSModel")

