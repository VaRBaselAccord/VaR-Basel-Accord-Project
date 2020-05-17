%%  VaR Basel Accord Project
%   Math Master Project 1
%   Christopher Burhan (19290751)
%   Master of Science (Actuarial and Financial Science)


% close all
% clear all
% clc

% S&P500 analysis
FirmClosePrice                                                              = readtable('EikonData.xlsx','Sheet','Sheet2','PreserveVariableNames',true);
FirmClosePriceadj                                                           = rmmissing(FirmClosePrice,2);
Dates                                                                       = datenum(FirmClosePrice.Date(2:end,:));

[number,TEXT,everything]                                                    = xlsread('EikonData.xlsx','Sheet3');
Ticker                                                                      = string(TEXT(2:506,8));
companylist                                                                 = TEXT(2:506,7);

WE                                                                          = 500;                                                                      % Estimation window

FirmNames                                                                   = FirmClosePriceadj.Properties.VariableNames(1:end);
FirmPrice                                                                   = FirmClosePriceadj(WE+3:end,FirmNames).Variables;
FirmPrice                                                                   = str2double(FirmPrice);
FirmReturn                                                                  = diff(log(FirmPrice));
FirmReturnTable                                                             = array2table(FirmReturn(1:end,:));
FirmReturnTable.Properties.VariableNames                                    = (FirmNames);
FirmRisk                                                                    = std(cell2mat(table2cell(FirmReturnTable)));
FirmRetn                                                                    = mean(cell2mat(table2cell(FirmReturnTable)));
FirmSharpeRatio                                                             = FirmRetn./FirmRisk;

% Sort Data
results                                                                     = table(FirmNames',(FirmRisk)',(FirmRetn)',(FirmSharpeRatio)','VariableNames',{'Firm','Risk','Return','SharpeRatio'});
results                                                                     = sortrows(results,{'SharpeRatio'},{'ascend'});

% Fitting distribution to Sharpe Ratio
pd                                                                          = fitdist(results.SharpeRatio,'normal');
h                                                                           = chi2gof(results.SharpeRatio,'CDF',pd);

if h == 0
    fprintf('Fail to reject null hypothesis\n\n')
else
    fprintf('Reject null hypothesis\n\n')
end

figure
histfit(results.SharpeRatio,[],'Normal')

% Risk treshold
HighRisk                                                                    = pd.mu + norminv(0.75) * pd.sigma;
LowRisk                                                                     = pd.mu + norminv(0.25) * pd.sigma;

xl1                                                                         = xline(LowRisk,'-.k','25%');
xl1.LineWidth                                                               = 2;
xl1.LabelVerticalAlignment                                                  = 'top';
xl1.LabelHorizontalAlignment                                                = 'left';

xl2                                                                         = xline(HighRisk,'-.k','75%');
xl2.LineWidth                                                               = 2;
xl2.LabelVerticalAlignment                                                  = 'top';
xl2.LabelHorizontalAlignment                                                = 'left';

xl3                                                                         = xline(pd.mu,'-.k','Average');
xl3.LineWidth                                                               = 2;
xl3.LabelVerticalAlignment                                                  = 'top';
xl3.LabelHorizontalAlignment                                                = 'left';


% Scale adjustment (daily to Annually)
start_date                                                                  = str2double(datestr(Dates(WE+2),'yyyy'));
end_date                                                                    = str2double(datestr(Dates(end),'yyyy')); 
scale                                                                       = 252*((end_date - start_date));

%%
% Extract firm moment & names
FirmCovar                                                                   = cov(cell2mat(table2cell(FirmReturnTable)));
FirmRskDaily                                                                = sqrt(diag(FirmCovar));
FirmRetnDaily                                                               = FirmRetn;

FirmRiskAnnual                                                              = sqrt(scale) * FirmRskDaily;
FirmReturnAnnual                                                            = scale * FirmRetnDaily;

figure;
scatter(FirmRiskAnnual,FirmReturnAnnual, 25, 'b', 'Filled');
hold on
for k = 1:length(FirmNames)
    text(FirmRiskAnnual(k) + 0.01, FirmReturnAnnual(k), FirmNames{k},'FontWeight','bold','FontSize', 11, 'Interpreter','Latex');
end
hold off
grid on
title('\textbf{S\&P500 constituents}','Interpreter','Latex','FontSize',12);
xlabel('Standard Deviation of Portfolio Returns (Annualized)','Interpreter','Latex','FontSize',12);
ylabel('Expected Returns (Annualized)','Interpreter','Latex','FontSize',12);
set(gca,'TickLabelInterpreter','latex','FontSize',12)

%% Risk averse Portfolio

RiskAverseLogic                                                             = results.SharpeRatio >= HighRisk;
FirmRiskAverse                                                              = datasample(results(RiskAverseLogic,:),5,'replace',false);
FirmRiskAverseNames                                                         = table2cell(FirmRiskAverse(:,1))';
FirmRiskAverseRsk                                                           = cell2mat(table2cell(FirmRiskAverse(:,2)));
FirmRiskAverseRetn                                                          = cell2mat(table2cell(FirmRiskAverse(:,3)));
FirmRiskAverseReturn                                                        = FirmReturnTable(:,FirmRiskAverseNames);
FirmRiskAverseCovar                                                         = cov(cell2mat(table2cell(FirmRiskAverseReturn)));

pRiskAverse                                                                 = Portfolio('AssetList',FirmRiskAverseNames,'RiskFreeRate',0);
pRiskAverse                                                                 = setAssetMoments(pRiskAverse,FirmRiskAverseRetn,FirmRiskAverseCovar); 

% Set weight constraints
AssetBoundsRiskAverse                                                       = [0, 1];
pRiskAverse                                                                 = setBounds(pRiskAverse, AssetBoundsRiskAverse(1), AssetBoundsRiskAverse(2));

% Set budget constraints (if the fund permits leverage up to 10% then the porfolio can be invested from 100% to 110% in risky assets otherwise 100% invested in risky assets).
pRiskAverse                                                                 = setBudget(pRiskAverse, 1, 1);

% Solver settings
pRiskAverse                                                                 = setSolver(pRiskAverse,'quadprog','Display','off','ConstraintTolerance',1.0e-8,'OptimalityTolerance',1.0e-8,'StepTolerance',1.0e-8,'MaxIterations',10000);

% Basis Comparison (i.e. equal weight portfolio)
pRiskAverse                                                                 = setInitPort(pRiskAverse,1/pRiskAverse.NumAssets);                                    % Set equal weight for initial portfolio
[eRiskAversersk,eRiskAverseret]                                             = estimatePortMoments(pRiskAverse,pRiskAverse.InitPort);

% Estimate the frontier line
pRiskAversewgt                                                              = estimateFrontier(pRiskAverse,100);
[pRiskAversersk,pRiskAverseret]                                             = estimatePortMoments(pRiskAverse,pRiskAversewgt);

% Set up optimisation problem which maximises sharpe ratio 
pRiskAverse                                                                 = setInitPort(pRiskAverse, 0);
sRiskAversewgt                                                              = estimateMaxSharpeRatio(pRiskAverse,'method','direct');
[sRiskAversersk,sRiskAverseret]                                             = estimatePortMoments(pRiskAverse,sRiskAversewgt);

% Extract asset moment & names
RiskAverseRskDaily                                                          = sqrt(diag(pRiskAverse.AssetCovar));
RiskAverseRetnDaily                                                         = pRiskAverse.AssetMean;
    
% Conversion from daily to annualy
RiskAverseRskAnnual                                                         = sqrt(scale) * RiskAverseRskDaily;
RiskAverseRetnAnnual                                                        = scale * RiskAverseRetnDaily;
portRiskAverseRskAnnual                                                     = sqrt(scale) * pRiskAversersk;
portRiskAverseRetnAnnual                                                    = scale * pRiskAverseret;
equalRiskAverseRskAnnual                                                    = sqrt(scale) * eRiskAversersk;
equalRiskAverseRetnAnnual                                                   = scale * eRiskAverseret;
sharpRiskAverseRskAnnual                                                    = sqrt(scale) * sRiskAversersk;
sharpRiskAverseRetnAnnual                                                   = scale * sRiskAverseret;

% Plot efficient frontier of risk averse portfolio that attains maximum Sharpe ratio 
figure
scatter(RiskAverseRskAnnual(sRiskAversewgt>0.001),RiskAverseRetnAnnual(sRiskAversewgt>0.001), 25, 'r', 'Filled');
hold on
scatter(RiskAverseRskAnnual(sRiskAversewgt<0.001),RiskAverseRetnAnnual(sRiskAversewgt<0.001), 25, 'k', 'Filled');
for k = 1:length(FirmRiskAverseNames)
    text(RiskAverseRskAnnual(k) + 0.01, RiskAverseRetnAnnual(k), FirmRiskAverseNames{k},'FontWeight','bold','FontSize', 11, 'Interpreter','Latex');
end
line(portRiskAverseRskAnnual,portRiskAverseRetnAnnual,'LineWidth',1.5,'color','b')
scatter(equalRiskAverseRskAnnual,equalRiskAverseRetnAnnual, 25, 'm', 'Filled')
text(equalRiskAverseRskAnnual + 0.01, equalRiskAverseRetnAnnual, "\textbf{equal weight}", 'Fontsize', 11,'Interpreter','Latex')
scatter(sharpRiskAverseRskAnnual,sharpRiskAverseRetnAnnual, 25, 'g', 'Filled')
text(sharpRiskAverseRskAnnual + 0.01, sharpRiskAverseRetnAnnual - 0.0025 , "\textbf{Max Sharpe ratio}", 'Fontsize', 11,'Interpreter','Latex')
hold off;
title('\textbf{Efficient Frontier of risk averse portfolio that attains Maximum Sharpe Ratio}','Interpreter','Latex','FontSize',12);
grid on
xlabel('Standard Deviation of Portfolio Returns (Annualized)','Interpreter','Latex','FontSize',12);
ylabel('Mean of Portfolio Returns (Annualized)','Interpreter','Latex','FontSize',12);
legend('Firm with weight more than 0','Firm with weight less than 0','Location','southeast')
set(legend,'Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex','FontSize',12)

% Display the weight
RiskAverseFirm                                                              = strings([505,1]);
RiskAverseFirm(1:length(FirmRiskAverseNames ))                              = FirmRiskAverseNames;

tf1                                                                         = cell(length(FirmRiskAverseNames),1);
for j1                                                                      = 1:length(FirmRiskAverseNames)
    tf1{j1}                                                                 = (string(FirmRiskAverseNames{j1}) == Ticker);
end

RiskAverseCompany                                                           = cell(length(tf1),1);
for l1 = 1:length(tf1)
    RiskAverseCompany{l1}                                                   = string(companylist(tf1{l1}));
end

RiskAversePort                                                              = table(FirmRiskAverseNames',RiskAverseCompany,sRiskAversewgt*100,'VariableNames',{'Ticker','Company Name','Weight (%)'});
disp(RiskAversePort);

% Display the return and risk for portfolio that maximises Sharpe Ratio)
fprintf ('Return for Max. Sharpe Ratio portfolio (High Risk) is %0.2f%%\n', sRiskAverseret*scale*100)
fprintf ('Risk for Max. Sharpe Ratio portfolio (High Risk) is %0.2f%%\n', sRiskAversersk*sqrt(scale)*100)

pRiskAversesratio                                                           = (portRiskAverseRetnAnnual - pRiskAverse.RiskFreeRate) ./ portRiskAverseRskAnnual;                                 
sRiskAversesratio                                                           = (sharpRiskAverseRetnAnnual - pRiskAverse.RiskFreeRate) / sharpRiskAverseRskAnnual; 

% clf;
figure
subplot(2,1,1);
scatter(RiskAverseRskAnnual(sRiskAversewgt>0.001),RiskAverseRetnAnnual(sRiskAversewgt>0.001), 25, 'r', 'Filled');
hold on
scatter(RiskAverseRskAnnual(sRiskAversewgt<0.001),RiskAverseRetnAnnual(sRiskAversewgt<0.001), 25, 'k', 'Filled');
for k = 1:length(FirmRiskAverseNames)
    text(RiskAverseRskAnnual(k) + 0.01, RiskAverseRetnAnnual(k), FirmRiskAverseNames{k},'FontWeight','bold','FontSize', 11, 'Interpreter','Latex');
end
line(portRiskAverseRskAnnual,portRiskAverseRetnAnnual,'LineWidth',1.5,'color','b')
scatter(equalRiskAverseRskAnnual,equalRiskAverseRetnAnnual, 25, 'm', 'Filled')
text(equalRiskAverseRskAnnual + 0.01, equalRiskAverseRetnAnnual, "\textbf{equal weight}", 'Fontsize', 11,'Interpreter','Latex')
scatter(sharpRiskAverseRskAnnual,sharpRiskAverseRetnAnnual, 50, 'g', 'Filled')
text(sharpRiskAverseRskAnnual + 0.01, sharpRiskAverseRetnAnnual - 0.0025 , "\textbf{Max Sharpe ratio}", 'Fontsize', 11,'Interpreter','Latex')
hold off;
title('\textbf{Efficient Frontier of risk averse portfolio that attains Maximum Sharpe Ratio}','Interpreter','Latex','FontSize',12);
grid on
xlabel('Standard Deviation of Portfolio Returns (Annualized)','Interpreter','Latex','FontSize',12);
ylabel('Mean of Portfolio Returns (Annualized)','Interpreter','Latex','FontSize',12);
legend('Firm with weight more than 0','Firm with weight less than 0','Location','southeast')
set(legend,'Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex','FontSize',12)
% xlim([min(portRiskAverseRskAnnual)-0.2 max(RiskAverseRskAnnual)+0.1])
% xlim([0.5 max(RiskAverseRskAnnual)+0.1])
subplot(2,1,2);
plot(portRiskAverseRskAnnual, pRiskAversesratio, 'b-','MarkerFaceColor','b', 'LineWidth', 2.0,'MarkerSize',3);
hold on
scatter(sharpRiskAverseRskAnnual, sRiskAversesratio,25,'g','Filled');
title('\textbf{Sharpe Ratio}','Interpreter','Latex','FontSize',12);
xlabel('Standard Deviation of Portfolio Return (Annualized)','Interpreter','Latex','FontSize',12);
% ylabel('Sharpe Ratio','Interpreter','Latex','FontSize',15);
set(gca,'TickLabelInterpreter','latex','FontSize',12)
% xlim([min(portRiskAverseRskAnnual)-0.2 max(RiskAverseRskAnnual)+0.1])
xlim([0.7 1.4])
grid on
hold off

% Display the Sharpe Ratio 
fprintf ('Sharpe Ratio for the opitmised portfolio (High Risk) is %0.3f\n', sRiskAversesratio)

%% Risk neutral Portfolio

RiskNeutralLogic                                                            = results.SharpeRatio >= LowRisk & results.SharpeRatio <= HighRisk;
FirmRiskNeutral                                                             = datasample(results(RiskNeutralLogic,:),5,'replace',false);   
FirmRiskNeutralNames                                                        = table2cell(FirmRiskNeutral(:,1))';
FirmRiskNeutralRsk                                                          = cell2mat(table2cell(FirmRiskNeutral(:,2)));
FirmRiskNeutralRetn                                                         = cell2mat(table2cell(FirmRiskNeutral(:,3)));
FirmRiskNeutralReturn                                                       = FirmReturnTable(:,FirmRiskNeutralNames);
FirmRiskNeutralCovar                                                        = cov(cell2mat(table2cell(FirmRiskNeutralReturn)));

pRiskNeutral                                                                = Portfolio('AssetList',FirmRiskNeutralNames,'RiskFreeRate',0);
pRiskNeutral                                                                = setAssetMoments(pRiskNeutral,FirmRiskNeutralRetn,FirmRiskNeutralCovar); 

% Set weight constraints
AssetBoundsRiskNeutral                                                      = [0, 1];
pRiskNeutral                                                                = setBounds(pRiskNeutral, AssetBoundsRiskNeutral(1), AssetBoundsRiskNeutral(2));

% Set budget constraints (if the fund permits leverage up to 10% then the porfolio can be invested from 100% to 110% in risky assets otherwise 100% invested in risky assets).
pRiskNeutral                                                                = setBudget(pRiskNeutral, 1, 1);

% Solver settings
pRiskNeutral                                                                = setSolver(pRiskNeutral,'quadprog','Display','off','ConstraintTolerance',1.0e-8,'OptimalityTolerance',1.0e-8,'StepTolerance',1.0e-8,'MaxIterations',10000);

% Basis Comparison (i.e. equal weight portfolio)
pRiskNeutral                                                                = setInitPort(pRiskNeutral,1/pRiskNeutral.NumAssets);                                    % Set equal weight for initial portfolio
[eRiskNeutralrsk,eRiskNeutralret]                                           = estimatePortMoments(pRiskNeutral,pRiskNeutral.InitPort);

% Estimate the frontier line
pRiskNeutralwgt                                                             = estimateFrontier(pRiskNeutral,100);
[pRiskNeutralrsk,pRiskNeutralret]                                           = estimatePortMoments(pRiskNeutral,pRiskNeutralwgt);

% Set up optimisation problem which maximises sharpe ratio 
pRiskNeutral                                                                = setInitPort(pRiskNeutral, 0);
sRiskNeutralwgt                                                             = estimateMaxSharpeRatio(pRiskNeutral,'method','direct');
[sRiskNeutralrsk,sRiskNeutralret]                                           = estimatePortMoments(pRiskNeutral,sRiskNeutralwgt);

% Extract asset moment & names
RiskNeutralRskDaily                                                         = sqrt(diag(pRiskNeutral.AssetCovar));
RiskNeutralRetnDaily                                                        = pRiskNeutral.AssetMean;

% Conversion from daily to annualy
RiskNeutralRskAnnual                                                        = sqrt(scale) * RiskNeutralRskDaily;
RiskNeutralRetnAnnual                                                       = scale * RiskNeutralRetnDaily;
portRiskNeutralRskAnnual                                                    = sqrt(scale) * pRiskNeutralrsk;
portRiskNeutralRetnAnnual                                                   = scale * pRiskNeutralret;
equalRiskNeutralRskAnnual                                                   = sqrt(scale) * eRiskNeutralrsk;
equalRiskNeutralRetnAnnual                                                  = scale * eRiskNeutralret;
sharpRiskNeutralRskAnnual                                                   = sqrt(scale) * sRiskNeutralrsk;
sharpRiskNeutralRetnAnnual                                                  = scale * sRiskNeutralret;
    
% Plot efficient frontier with portfolio that attains maximum Sharpe ratio 
figure
scatter(RiskNeutralRskAnnual(sRiskNeutralwgt>0.001),RiskNeutralRetnAnnual(sRiskNeutralwgt>0.001), 25, 'r', 'Filled');
hold on
scatter(RiskNeutralRskAnnual(sRiskNeutralwgt<0.001),RiskNeutralRetnAnnual(sRiskNeutralwgt<0.001), 25, 'k', 'Filled');
for k = 1:length(FirmRiskNeutralNames)
    text(RiskNeutralRskAnnual(k) + 0.01, RiskNeutralRetnAnnual(k), FirmRiskNeutralNames{k},'FontWeight','bold','FontSize', 11, 'Interpreter','Latex');
end
line(portRiskNeutralRskAnnual,portRiskNeutralRetnAnnual,'LineWidth',1.5,'color','b')
scatter(equalRiskNeutralRskAnnual,equalRiskNeutralRetnAnnual, 25, 'm', 'Filled')
text(equalRiskNeutralRskAnnual + 0.01, equalRiskNeutralRetnAnnual, "\textbf{equal weight}", 'Fontsize', 11,'Interpreter','Latex')
scatter(sharpRiskNeutralRskAnnual,sharpRiskNeutralRetnAnnual, 25, 'g', 'Filled')
text(sharpRiskNeutralRskAnnual + 0.01, sharpRiskNeutralRetnAnnual - 0.0025 , "\textbf{Max Sharpe ratio}", 'Fontsize', 11,'Interpreter','Latex')
hold off;
title('\textbf{Efficient Frontier risk neutral portfolio that attains Maximum Sharpe Ratio}','Interpreter','Latex','FontSize',12);
grid on
xlabel('Standard Deviation of Portfolio Returns (Annualized)','Interpreter','Latex','FontSize',12);
ylabel('Mean of Portfolio Returns (Annualized)','Interpreter','Latex','FontSize',12);
legend('Firm with weight more than 0','Firm with weight less than 0','Location','southeast')
set(legend,'Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex','FontSize',12)

% Display the weight
RiskNeutralFirm                                                             = strings([505,1]);
RiskNeutralFirm(1:length(FirmRiskNeutralNames))                             = FirmRiskNeutralNames;

tf2                                                                         = cell(length(FirmRiskNeutralNames),1);
for j2                                                                      = 1:length(FirmRiskNeutralNames)
    tf2{j2}                                                                 = (string(FirmRiskNeutralNames{j2}) == Ticker);
end

RiskNeutralCompany                                                          = cell(length(tf2),1);
for l2                                                                      = 1:length(tf2)
    RiskNeutralCompany{l2}                                                  = string(companylist(tf2{l2}));
end

RiskNeutralPort                                                             = table(FirmRiskNeutralNames',RiskNeutralCompany,sRiskNeutralwgt*100,'VariableNames',{'Ticker','Company Name','Weight (%)'});
disp(RiskNeutralPort);

% Display the return and risk for portfolio that maximises Sharpe Ratio)
fprintf ('Return for Max. Sharpe Ratio portfolio (Middle Risk) is %0.2f%%\n', sRiskNeutralret*scale*100)
fprintf ('Risk for Max. Sharpe Ratio portfolio (Middle Risk) is %0.2f%%\n', sRiskNeutralrsk*sqrt(scale)*100)

pRiskNeutralsratio                                                          = (portRiskNeutralRetnAnnual - pRiskNeutral.RiskFreeRate) ./ portRiskNeutralRskAnnual;                                 
sRiskNeutralsratio                                                          = (sharpRiskNeutralRetnAnnual - pRiskNeutral.RiskFreeRate) / sharpRiskNeutralRskAnnual; 

% clf;
figure
subplot(2,1,1);
scatter(RiskNeutralRskAnnual(sRiskNeutralwgt>0.001),RiskNeutralRetnAnnual(sRiskNeutralwgt>0.001), 25, 'r', 'Filled');
hold on
scatter(RiskNeutralRskAnnual(sRiskNeutralwgt<0.001),RiskNeutralRetnAnnual(sRiskNeutralwgt<0.001), 25, 'k', 'Filled');
for k = 1:length(FirmRiskNeutralNames)
    text(RiskNeutralRskAnnual(k) + 0.01, RiskNeutralRetnAnnual(k), FirmRiskNeutralNames{k},'FontWeight','bold','FontSize', 11, 'Interpreter','Latex');
end
line(portRiskNeutralRskAnnual,portRiskNeutralRetnAnnual,'LineWidth',1.5,'color','b')
scatter(equalRiskNeutralRskAnnual,equalRiskNeutralRetnAnnual, 25, 'm', 'Filled')
text(equalRiskNeutralRskAnnual + 0.01, equalRiskNeutralRetnAnnual, "\textbf{equal weight}", 'Fontsize', 11,'Interpreter','Latex')
scatter(sharpRiskNeutralRskAnnual,sharpRiskNeutralRetnAnnual, 25, 'g', 'Filled')
text(sharpRiskNeutralRskAnnual + 0.01, sharpRiskNeutralRetnAnnual - 0.0025 , "\textbf{Max Sharpe ratio}", 'Fontsize', 11,'Interpreter','Latex')
hold off;
title('\textbf{Efficient Frontier of risk neutral portfolio that attains Maximum Sharpe Ratio}','Interpreter','Latex','FontSize',12);
grid on
xlabel('Standard Deviation of Portfolio Returns (Annualized)','Interpreter','Latex','FontSize',12);
ylabel('Mean of Portfolio Return (Annualized)','Interpreter','Latex','FontSize',12);

subplot(2,1,2);
plot(portRiskNeutralRskAnnual, pRiskNeutralsratio, 'b-','MarkerFaceColor','b', 'LineWidth', 2.0,'MarkerSize',3);
hold on
scatter(sharpRiskNeutralRskAnnual, sRiskNeutralsratio,25,'g','Filled');
title('\textbf{Sharpe Ratio}','Interpreter','Latex','FontSize',12);
xlabel('Standard Deviation of Portfolio Return (Annualized)','Interpreter','Latex','FontSize',12);
ylabel('Sharpe Ratio','Interpreter','Latex','FontSize',12);
set(gca,'TickLabelInterpreter','latex','FontSize',12)
grid on
hold off

% Display the Sharpe Ratio 
fprintf ('Sharpe Ratio for the opitmised portfolio (Middle Risk) is %0.3f\n', sRiskNeutralsratio)

%% Risk Seeking Portfolio

RiskSeekingLogic                                                            = results.SharpeRatio <= LowRisk;
FirmRiskSeeking                                                             = datasample(results(RiskSeekingLogic,:),5,'replace',false);
FirmRiskSeekingNames                                                        = table2cell(FirmRiskSeeking(:,1))';
FirmRiskSeekingRsk                                                          = cell2mat(table2cell(FirmRiskSeeking(:,2)));
FirmRiskSeekingRetn                                                         = cell2mat(table2cell(FirmRiskSeeking(:,3)));
FirmRiskSeekingReturn                                                       = FirmReturnTable(:,FirmRiskSeekingNames);
FirmRiskSeekingCovar                                                        = cov(cell2mat(table2cell(FirmRiskSeekingReturn)));

pRiskSeeking                                                                = Portfolio('AssetList',FirmRiskSeekingNames,'RiskFreeRate',0);
pRiskSeeking                                                                = setAssetMoments(pRiskSeeking,FirmRiskSeekingRetn,FirmRiskSeekingCovar); 

% Set weight constraints
AssetBoundsRiskSeeking                                                      = [0, 1];
pRiskSeeking                                                                = setBounds(pRiskSeeking, AssetBoundsRiskSeeking(1), AssetBoundsRiskSeeking(2));

% Set budget constraints (if the fund permits leverage up to 10% then the porfolio can be invested from 100% to 110% in risky assets otherwise 100% invested in risky assets).
pRiskSeeking                                                                = setBudget(pRiskSeeking, 1, 1);

% Solver settings
pRiskSeeking                                                                = setSolver(pRiskSeeking,'quadprog','Display','off','ConstraintTolerance',1.0e-8,'OptimalityTolerance',1.0e-8,'StepTolerance',1.0e-8,'MaxIterations',10000);

% Basis Comparison (i.e. equal weight portfolio)
pRiskSeeking                                                                = setInitPort(pRiskSeeking,1/pRiskSeeking.NumAssets);                                    % Set equal weight for initial portfolio
[eRiskSeekingrsk,eRiskSeekingret]                                           = estimatePortMoments(pRiskSeeking,pRiskSeeking.InitPort);

% Estimate the frontier line
pRiskSeekingwgt                                                             = estimateFrontier(pRiskSeeking,100);
[pRiskSeekingrsk,pRiskSeekingret]                                           = estimatePortMoments(pRiskSeeking,pRiskSeekingwgt);

% Set up optimisation problem which maximises sharpe ratio 
pRiskSeeking                                                                = setInitPort(pRiskSeeking, 0);
sRiskSeekingwgt                                                             = estimateMaxSharpeRatio(pRiskSeeking,'method','direct');
[sRiskSeekingrsk,sRiskSeekingret]                                           = estimatePortMoments(pRiskSeeking,sRiskSeekingwgt);

% Extract asset moment & names
RiskSeekingRskDaily                                                         = sqrt(diag(pRiskSeeking.AssetCovar));
RiskSeekingRetnDaily                                                        = pRiskSeeking.AssetMean;

% Conversion from daily to annualy
RiskSeekingRskAnnual                                                        = sqrt(scale) * RiskSeekingRskDaily;
RiskSeekingRetnAnnual                                                       = scale * RiskSeekingRetnDaily;
portRiskSeekingRskAnnual                                                    = sqrt(scale) * pRiskSeekingrsk;
portRiskSeekingRetnAnnual                                                   = scale * pRiskSeekingret;
equalRiskSeekingRskAnnual                                                   = sqrt(scale) * eRiskSeekingrsk;
equalRiskSeekingRetnAnnual                                                  = scale * eRiskSeekingret;
sharpRiskSeekingRskAnnual                                                   = sqrt(scale) * sRiskSeekingrsk;
sharpRiskSeekingRetnAnnual                                                  = scale * sRiskSeekingret;

% Plot efficient frontier with portfolio that attains maximum Sharpe ratio 
figure
scatter(RiskSeekingRskAnnual(sRiskSeekingwgt>0.001),RiskSeekingRetnAnnual(sRiskSeekingwgt>0.001), 25, 'r', 'Filled');
hold on
scatter(RiskSeekingRskAnnual(sRiskSeekingwgt<0.001),RiskSeekingRetnAnnual(sRiskSeekingwgt<0.001), 25, 'k', 'Filled');
for k = 1:length(FirmRiskSeekingNames)
    text(RiskSeekingRskAnnual(k) + 0.01, RiskSeekingRetnAnnual(k), FirmRiskSeekingNames{k},'FontWeight','bold','FontSize', 11, 'Interpreter','Latex');
end
line(portRiskSeekingRskAnnual,portRiskSeekingRetnAnnual,'LineWidth',1.5,'color','b')
scatter(equalRiskSeekingRskAnnual,equalRiskSeekingRetnAnnual, 25, 'r', 'Filled')
text(equalRiskSeekingRskAnnual + 0.01, equalRiskSeekingRetnAnnual, "\textbf{equal weight}", 'Fontsize', 11,'Interpreter','Latex')
scatter(sharpRiskSeekingRskAnnual,sharpRiskSeekingRetnAnnual, 25, 'g', 'Filled')
text(sharpRiskSeekingRskAnnual + 0.01, sharpRiskSeekingRetnAnnual - 0.08 , "\textbf{Max Sharpe ratio}", 'Fontsize', 11,'Interpreter','Latex')
hold off;
title('\textbf{Efficient Frontier of risk seeking portfolio that attains Maximum Sharpe Ratio}','Interpreter','Latex','FontSize',12);
grid on
xlabel('Standard Deviation of Portfolio Returns (Annualized)','Interpreter','Latex','FontSize',12);
ylabel('Mean of Portfolio Returns (Annualized)','Interpreter','Latex','FontSize',12);
legend('Firm with weight more than 0','Firm with weight less than 0','Location','southwest')
set(legend,'Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex','FontSize',12)

% Display the weight
RiskSeekingFirm                                                             = strings([505,1]);
RiskSeekingFirm(1:length(FirmRiskSeekingNames))                             = FirmRiskSeekingNames;

tf3                                                                         = cell(length(FirmRiskSeekingNames),1);
for j3                                                                      = 1:length(FirmRiskSeekingNames)
    tf3{j3}                                                                 = (string(FirmRiskSeekingNames{j3}) == Ticker);
end

RiskSeekingCompany                                                          = cell(length(tf3),1);
for l3 = 1:length(tf3)
    RiskSeekingCompany{l3}                                                  = string(companylist(tf3{l3}));
end

RiskSeekingPort                                                             = table(FirmRiskSeekingNames',RiskSeekingCompany,sRiskSeekingwgt*100,'VariableNames',{'Ticker','Company Name','Weight (%)'});
disp(RiskSeekingPort);

% Display the return and risk for portfolio that maximises Sharpe Ratio)
fprintf ('Return for Max. Sharpe Ratio portfolio (Low Risk) is %0.2f%%\n', sRiskSeekingret*scale*100)
fprintf ('Risk for Max. Sharpe Ratio portfolio (Low Risk) is %0.2f%%\n', sRiskSeekingrsk*sqrt(scale)*100)

pRiskSeekingsratio                                                          = (portRiskSeekingRetnAnnual - pRiskSeeking.RiskFreeRate) ./ portRiskSeekingRskAnnual;                                 
sRiskSeekingsratio                                                          = (sharpRiskSeekingRetnAnnual - pRiskSeeking.RiskFreeRate) / sharpRiskSeekingRskAnnual; 

% clf;
figure
subplot(2,1,1);
scatter(RiskSeekingRskAnnual(sRiskSeekingwgt>0.001),RiskSeekingRetnAnnual(sRiskSeekingwgt>0.001), 25, 'r', 'Filled');
hold on
scatter(RiskSeekingRskAnnual(sRiskSeekingwgt<0.001),RiskSeekingRetnAnnual(sRiskSeekingwgt<0.001), 25, 'k', 'Filled');
for k = 1:length(FirmRiskSeekingNames)
    text(RiskSeekingRskAnnual(k) + 0.01, RiskSeekingRetnAnnual(k), FirmRiskSeekingNames{k},'FontWeight','bold','FontSize', 11, 'Interpreter','Latex');
end
line(portRiskSeekingRskAnnual,portRiskSeekingRetnAnnual,'LineWidth',1.5,'color','b')
scatter(equalRiskSeekingRskAnnual,equalRiskSeekingRetnAnnual, 25, 'r', 'Filled')
text(equalRiskSeekingRskAnnual + 0.01, equalRiskSeekingRetnAnnual, "\textbf{equal weight}", 'Fontsize', 11,'Interpreter','Latex')
scatter(sharpRiskSeekingRskAnnual,sharpRiskSeekingRetnAnnual, 25, 'g', 'Filled')
text(sharpRiskSeekingRskAnnual + 0.01, sharpRiskSeekingRetnAnnual - 0.08 , "\textbf{Max Sharpe ratio}", 'Fontsize', 11,'Interpreter','Latex')
hold off;
title('\textbf{Efficient Frontier of risk seeking portfolio that attains Maximum Sharpe Ratio}','Interpreter','Latex','FontSize',12);
grid on
xlabel('Standard Deviation of Portfolio Returns (Annualized)','Interpreter','Latex','FontSize',12);
ylabel('Mean of Portfolio Returns (Annualized)','Interpreter','Latex','FontSize',12);
legend('Firm with weight more than 0','Firm with weight less than 0','Location','southwest')
set(legend,'Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex','FontSize',12)
% xlim([(min(portRiskSeekingRskAnnual))-0.2 (max(RiskSeekingRskAnnual))+0.1])

subplot(2,1,2);
plot(portRiskSeekingRskAnnual, pRiskSeekingsratio, 'b-','MarkerFaceColor','b', 'LineWidth', 2.0,'MarkerSize',3);
hold on
scatter(sharpRiskSeekingRskAnnual, sRiskSeekingsratio,25,'g','Filled');
title('\textbf{Sharpe Ratio}','Interpreter','Latex','FontSize',12);
xlabel('Standard Deviation of Portfolio Return (Annualized)','Interpreter','Latex','FontSize',12);
% ylabel('Sharpe Ratio','Interpreter','Latex','FontSize',15);
% set(legend,'Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex','FontSize',12)
xlim([0.6 2.2])
grid on
hold off

% Display the Sharpe Ratio 
fprintf ('Sharpe Ratio for the opitmised portfolio (Low Risk) is %0.3f\n', sRiskSeekingsratio)

%% Summary of optimum portfolio weight
clc

% Display the weight of all portfolios
fprintf('\n\nRisk averse portfolio\n\n');
disp(RiskAversePort);
fprintf('\n\nRisk neutral portfolio\n\n');
disp(RiskNeutralPort);
fprintf('\n\nRisk seeking portfolio\n\n');
disp(RiskSeekingPort);

%% Running Bootstrap and Determining VaR @ 1%
n                                                                           = 1000;                                                                           % Number of bootstrap replication
PortfolioFirmPrice                                                          = FirmClosePriceadj(2:end,FirmNames).Variables;
PortfolioFirmPrice                                                          = str2double(PortfolioFirmPrice);
PortfolioFirmReturn                                                         = diff(log(PortfolioFirmPrice));
T                                                                           = length(PortfolioFirmReturn);                                                    % Number of observation for stock return 
prob                                                                        = 0.01;                                                                           % Probability
op                                                                          = WE*prob;                                                                        % Location of VaR ((WE*p)th value) 
expectedviolation                                                           = prob*(T-WE);

% Expected violations logic
if expectedviolation - floor(expectedviolation) >= 0.5
    exp = ceil(expectedviolation);
else
    exp = floor(expectedviolation);
end

% Extract portfolio return 
PortfolioReturn                                                             = PortfolioFirmReturn';
PortfolioReturnTable                                                        = table(FirmNames',PortfolioReturn,FirmRisk',FirmSharpeRatio','VariableNames',{'Firm','Return','Risk','SharpeRatio'});
PortfolioReturnTable                                                        = sortrows(PortfolioReturnTable,{'SharpeRatio'},{'ascend'});
PortfolioReturnFirm                                                         = string(PortfolioReturnTable.Firm);

% High Risk Portfolio
tf4                                                                         = cell(length(FirmRiskAverseNames),1);
for j4                                                                      = 1:length(FirmRiskAverseNames)
    tf4{j4}                                                                 = string(FirmRiskAverseNames{j4}) == PortfolioReturnFirm;
end

for l4                                                                      = 1:length(tf4)
    RiskAversePortfolioReturnTable(l4,:)                                    = PortfolioReturnTable(tf4{l4},:);
end
RiskAversePortReturn                                                        = sRiskAversewgt .* cell2mat(table2cell(RiskAversePortfolioReturnTable(:,2)));

% Middle Risk Portfolio
tf5                                                                         = cell(length(FirmRiskNeutralNames),1);
for j5                                                                      = 1:length(FirmRiskNeutralNames)
    tf5{j5}                                                                 = string(FirmRiskNeutralNames{j5}) == PortfolioReturnFirm;
end

for l5 = 1:length(tf5)
    RiskNeutralPortfolioReturnTable(l5,:)                                   = PortfolioReturnTable(tf5{l5},:);
end
RiskNeutralPortReturn                                                       = sRiskNeutralwgt .* cell2mat(table2cell(RiskNeutralPortfolioReturnTable(:,2)));

% Low Risk Portfolio
tf6                                                                         = cell(length(FirmRiskSeekingNames),1);
for j6                                                                      = 1:length(FirmRiskSeekingNames)
    tf6{j6}                                                                 = string(FirmRiskSeekingNames{j6}) == PortfolioReturnFirm;
end

for l6 = 1:length(tf6)
    RiskSeekingPortfolioReturnTable(l6,:)                                   = PortfolioReturnTable(tf6{l6},:);
end
RiskSeekingPortReturn                                                       = sRiskSeekingwgt .* cell2mat(table2cell(RiskSeekingPortfolioReturnTable(:,2)));

% Portfolio Return
PortReturn                                                                  = cell(1,3);
PortReturn{1}                                                               = (sum(RiskAversePortReturn))';
PortReturn{2}                                                               = (sum(RiskNeutralPortReturn))';
PortReturn{3}                                                               = (sum(RiskSeekingPortReturn))';

%% Risk averse Portfolio
tic
parfor  t                                                                   = WE + 1:T                                              
        t1                                                                  = t - WE;                                                                                  % Start of data window
        t2                                                                  = t - 1;                                                                                   % End of data window
        window                                                              = PortReturn{1}(t1:t2);                                                                    % Data for estimation     
        rtcell                                                              = cell(n,1);                                                    % create empty cell for rt 
         VaRcell                                                            = cell(n,1);                                                    % create empty cell for VaR forecast   
            for i                                                           = 1:n
                y                                                           = sort(datasample(window,length(window)));                      % Sort the uniformly pick data sample from lowest to highest
                rtcell{i}                                                   = y                                                             % Put them in different rows
                VaRcell{i}                                                  = rtcell{i}(op,:);                                              % Find the VaR from each row (pick the 1% value)
                VaRforecastRiskAverse{t}                                    = mean(cell2mat(VaRcell))                                       % Mean of the VaR i.e. final VaR value
            end      
end
toc

% Backtest 
% Plotting the stock returns and VaR forecast                                       
cRiskAverse                                                                 = VaRforecastRiskAverse;                                       
cRiskAverse(:,1:WE)                                                         = [];  
bRiskAverse                                                                 = {(WE+1+1:T+1)' (cRiskAverse)'};
bRiskAverse{:,1}                                                            = (WE+1+1:T+1)';
xax2RiskAverse                                                              = Dates(bRiskAverse{:,1});                                                   
bRiskAverse{:,2}                                                            = cell2mat(bRiskAverse{2});

VaRRiskAverse250                                                            = bRiskAverse{:,2};
RiskAversereturnsadjusted                                                   = PortReturn{1}; 

% VaR violation throughout VaR forecast window
RiskAversenov250                                                            = length(find(RiskAversereturnsadjusted(WE+1:T)<VaRRiskAverse250));                         % Total Number of Violation throughout VaR forecast window
RiskAverseloc_violation2_250                                                = find(RiskAversereturnsadjusted(WE+1:T)<VaRRiskAverse250);                                 % Locations where violation occur
RiskAverseVR250                                                             = length(find(RiskAversereturnsadjusted(WE+1:T)<VaRRiskAverse250))/exp;                     % Violation Ratios = Observed no of siolations / Expected no of violations 

figure;
set(groot,'DefaultAxesTickLabelInterpreter','Tex');
hold on
plot(xax2RiskAverse,RiskAversereturnsadjusted(WE+1:T));
plot(xax2RiskAverse,VaRRiskAverse250,'LineWidth',1.25);
yline(0,'-.r','LineWidth',1.5);

% Plot configuration
ytickformat('%.2f');
ylim([-0.2 0.2]);
ylabel('Return','Interpreter','latex');
legend('Portfolio returns','VaR forecast','Location','northwest');
title('Backtest risk averse','Interpreter','Latex','FontSize',12);
datetick('x','mmm-yyyy','keeplimits');
set(legend,'Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex','FontSize',12)

% Number of VaR violation
RiskAverseComparison                                                        = cell2mat({PortReturn{1}(WE+1:T) cell2mat(cRiskAverse)'});
ViolationRiskAverse                                                         = cell(1,T);
for     x                                                                   = WE+1+250:T
        x1                                                                  = x-(WE+250);
        x2                                                                  = x-WE;
        ViolationRiskAverse{x}                                              = sum(RiskAverseComparison(x1:x2,1)<RiskAverseComparison(x1:x2,2));
        ViolationRiskAverse(:,WE+1:WE+250)                                  = {0};
end

% Plot of cumulative violations over the previous 250 days
figure;
plot(Dates(WE+1:T),cell2mat(ViolationRiskAverse)','LineWidth',1.5)
ytickformat('%.2f');
ylabel('Cumulative violations over previous 250 days','Interpreter','latex');
title('Risk averse portfolio','Interpreter','Latex','FontSize',12);
datetick('x','mmm-yyyy');
set(gca,'TickLabelInterpreter','latex','FontSize',12)

% Determining multiplication factor "k"
kRiskAverse                                                                 = cell2mat(cell(1,T));
for     t                                                                   = WE+1:T
    if ViolationRiskAverse{t} <= 4
        kRiskAverse(t)  = 0;

    else
        if ViolationRiskAverse{t} == 5
            kRiskAverse(t)   = 0.40;

        else
            if ViolationRiskAverse{t} == 6
                kRiskAverse(t)   = 0.50;

            else
                if ViolationRiskAverse{t} == 7
                    kRiskAverse(t) = 0.65;

                else
                    if ViolationRiskAverse{t} == 8
                        kRiskAverse(t)   = 0.75;

                    else
                        if ViolationRiskAverse{t} == 9
                            kRiskAverse(t)   = 0.85;

                        else
                            if ViolationRiskAverse{t} >= 10
                                kRiskAverse(t)   = 1.00;  
                            end

                        end

                    end

                end

            end

        end

    end
end

% Average VaR over the last 60 days
AverageVaRRiskAverse                                                        = cell(1,T);
for     y                                                                   = WE+1+59:T
        y1                                                                  = y-(WE+59);
        y2                                                                  = y-WE;
        AverageVaRRiskAverse{y}                                             = mean(RiskAverseComparison(y1:y2,2));
        AverageVaRRiskAverse(:,WE+1:WE+59)                                  = {0};
        
end

AverageVaRRiskAverse60                                                      = cellfun(@times, AverageVaRRiskAverse, num2cell(3+kRiskAverse), 'UniformOutput', false);

% Capital Charge 
previousdayVaRforecastRiskAverse                                            = cell2mat(VaRforecastRiskAverse');
previousdayVaRforecastRiskAverse(1,:)                                       = [];
previousdayVaRforecastRiskAverse(2370,1)                                    = 0;
CapitalChargeRiskAverse                                                     = mean(abs(min(cell2mat({cell2mat(AverageVaRRiskAverse60') previousdayVaRforecastRiskAverse}),[],2)));

% Proportion out of green zone
RiskAverseoutofGreenZone                                                    = (sum(kRiskAverse(WE+1:T) > 0.4)/length(kRiskAverse(WE+1:T)))*100;

% Print result 
formatSpec1                                                                 = 'Average Capital Charge is %5.4f%%, %.2f%% out of green zone.\n';
fprintf(formatSpec1,CapitalChargeRiskAverse*100,RiskAverseoutofGreenZone);

%% Risk neutral Portfolio
tic
parfor  t                                                                   = WE + 1:T                                              
        t1                                                                  = t - WE;                                                                                 % Start of data window
        t2                                                                  = t - 1;                                                                                  % End of data window
        window                                                              = PortReturn{2}(t1:t2);                                                                   % Data for estimation     
        rtcell                                                              = cell(n,1);                                                 % create empty cell for rt 
         VaRcell                                                            = cell(n,1);                                                 % create empty cell for VaR forecast   
            for i                                                           = 1:n
                y                                                           = sort(datasample(window,length(window)));                   % Sort the uniformly pick data sample from lowest to highest
                rtcell{i}                                                   = y                                                          % Put them in different rows
                VaRcell{i}                                                  = rtcell{i}(op,:);                                           % Find the VaR from each row (pick the 1% value)
                VaRforecastRiskNeutral{t}                                   = mean(cell2mat(VaRcell))                                    % Mean of the VaR i.e. final VaR value
            end
              
end
toc

% Backtest 
                                    
cRiskNeutral                                                                = VaRforecastRiskNeutral;                                       
cRiskNeutral(:,1:WE)                                                        = [];  
bRiskNeutral                                                                = {(WE+1+1:T+1)' (cRiskNeutral)'};
bRiskNeutral{:,1}                                                           = (WE+1+1:T+1)';
xax2RiskNeutral                                                             = Dates(bRiskNeutral{:,1});                                                   
bRiskNeutral{:,2}                                                           = cell2mat(bRiskNeutral{2});

VaRRiskNeutral250                                                           = bRiskNeutral{:,2};
RiskNeutralreturnsadjusted                                                  = PortReturn{2}; 

% VaR violation throughout VaR forecast window
RiskNeutralnov250                                                           = length(find(RiskNeutralreturnsadjusted(WE+1:T)<VaRRiskNeutral250));                         % Total number of VaR Violation throughout VaR forecast window
RiskNeutralloc_violation2_250                                               = find(RiskNeutralreturnsadjusted(WE+1:T)<VaRRiskNeutral250);                                 % Locations where violation occur
RiskNeutralVR250                                                            = length(find(RiskNeutralreturnsadjusted(WE+1:T)<VaRRiskNeutral250))/exp;                     % Violation Ratios = Observed no of siolations / Expected no of violations

figure;
set(groot,'DefaultAxesTickLabelInterpreter','Tex');
hold on
plot(xax2RiskNeutral,RiskNeutralreturnsadjusted(WE+1:T));
plot(xax2RiskNeutral,VaRRiskNeutral250,'LineWidth',1.25);
yline(0,'-.r','LineWidth',1.5);

% Plot configuration
ytickformat('%.2f');
ylim([-0.2 0.2]);
ylabel('Return','Interpreter','latex');
legend('Portfolio returns','VaR forecast','Location','northwest');
title('Backtest risk neutral portfolio','Interpreter','Latex','FontSize',12);
datetick('x','mmm-yyyy','keeplimits');
set(legend,'Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex','FontSize',12)

% Number of VaR violation
RiskNeutralComparison                                                       = cell2mat({PortReturn{2}(WE+1:T) cell2mat(cRiskNeutral)'});
ViolationRiskNeutral                                                        = cell(1,T);
for     x                                                                   = WE+1+250:T
        x1                                                                  = x-(WE+250);
        x2                                                                  = x-WE;
        ViolationRiskNeutral{x}                                             = sum(RiskNeutralComparison(x1:x2,1)<RiskNeutralComparison(x1:x2,2));
        ViolationRiskNeutral(:,WE+1:WE+250)                                 = {0};
end

% Plot of cumulative violations over the previous 250 days
figure;
plot(Dates(WE+1:T),cell2mat(ViolationRiskNeutral)','LineWidth',1.5)
ytickformat('%.2f');
ylabel('Cumulative violations over previous 250 days','Interpreter','latex');
title('Risk neutral portfolio','Interpreter','Latex','FontSize',12);
datetick('x','mmm-yyyy');
set(gca,'TickLabelInterpreter','latex','FontSize',12)

% Determining multiplication factor "k"
kRiskNeutral                                                                = cell2mat(cell(1,T));
for     t                                                                   = WE+1:T
    if ViolationRiskNeutral{t} <= 4
        kRiskNeutral(t)  = 0;

    else
        if ViolationRiskNeutral{t} == 5
            kRiskNeutral(t)   = 0.40;

        else
            if ViolationRiskNeutral{t} == 6
                kRiskNeutral(t)   = 0.50;

            else
                if ViolationRiskNeutral{t} == 7
                    kRiskNeutral(t) = 0.65;

                else
                    if ViolationRiskNeutral{t} == 8
                        kRiskNeutral(t)   = 0.75;

                    else
                        if ViolationRiskNeutral{t} == 9
                            kRiskNeutral(t)   = 0.85;

                        else
                            if ViolationRiskNeutral{t} >= 10
                                kRiskNeutral(t)   = 1.00;  
                            end

                        end

                    end

                end

            end

        end

    end
end

% Average VaR over the last 60 days
AverageVaRRiskNeutral                                                       = cell(1,T);
for     y                                                                   = WE+1+59:T
        y1                                                                  = y-(WE+59);
        y2                                                                  = y-WE;
        AverageVaRRiskNeutral{y}                                            = mean(RiskNeutralComparison(y1:y2,2));
        AverageVaRRiskNeutral(:,WE+1:WE+59)                                 = {0};
        
end

AverageVaRRiskNeutral60                                                     = cellfun(@times, AverageVaRRiskNeutral, num2cell(3+kRiskNeutral), 'UniformOutput', false);

% Capital Charge 
previousdayVaRforecastRiskNeutral                                           = cell2mat(VaRforecastRiskNeutral');
previousdayVaRforecastRiskNeutral(1,:)                                      = [];
previousdayVaRforecastRiskNeutral(2370,1)                                   = 0;
CapitalChargeRiskNeutral                                                    = mean(abs(min(cell2mat({cell2mat(AverageVaRRiskNeutral60') previousdayVaRforecastRiskNeutral}),[],2)));

% Proportion out of green zone
RiskNeutraloutofGreenZone                                                   = (sum(kRiskNeutral(WE+1:T) > 0.4)/length(kRiskNeutral(WE+1:T)))*100;

% Print result 
formatSpec1                                                                 = 'Average Capital Charge is %5.4f%%, %.2f%% out of green zone.\n';
fprintf(formatSpec1,CapitalChargeRiskNeutral*100,RiskNeutraloutofGreenZone);

%% Risk seeking Portfolio
tic
parfor  t                                                                   = WE + 1:T                                              
        t1                                                                  = t - WE;                                                                                 % Start of data window
        t2                                                                  = t - 1;                                                                                  % End of data window
        window                                                              = PortReturn{3}(t1:t2);                                                                   % Data for estimation     
        rtcell                                                              = cell(n,1);                                                    % create empty cell for rt 
         VaRcell                                                            = cell(n,1);                                                    % create empty cell for VaR forecast   
            for i                                                           = 1:n
                y                                                           = sort(datasample(window,length(window)));                      % Sort the uniformly pick data sample from lowest to highest
                rtcell{i}                                                   = y                                                             % Put them in different rows
                VaRcell{i}                                                  = rtcell{i}(op,:);                                              % Find the VaR from each row (pick the 1% value)
                VaRforecastRiskSeeking{t}                                   = mean(cell2mat(VaRcell))                                       % Mean of the VaR i.e. final VaR value
            end      
end
toc

% Backtest after the first 250 business days
                                    
cRiskSeeking                                                                = VaRforecastRiskSeeking;                                       
cRiskSeeking(:,1:WE)                                                        = [];  
bRiskSeeking                                                                = {(WE+1+1:T+1)' (cRiskSeeking)'};
bRiskSeeking{:,1}                                                           = (WE+1+1:T+1)';
xax2RiskSeeking                                                             = Dates(bRiskSeeking{:,1});                                                           
bRiskSeeking{:,2}                                                           = cell2mat(bRiskSeeking{2});

VaRRiskSeeking250                                                           = bRiskSeeking{:,2};
RiskSeekingreturnsadjusted                                                  = PortReturn{3}; 

% VaR violation throughout VaR forecast window
RiskSeekingnov250                                                           = length(find(RiskSeekingreturnsadjusted(WE+1:T)<VaRRiskSeeking250));                         % Total number of Violation throughout VaR forecast window
RiskSeekingloc_violation2_250                                               = find(RiskSeekingreturnsadjusted(WE+1:T)<VaRRiskSeeking250);                                 % Locations where violation occur
RiskSeekingVR250                                                            = length(find(RiskSeekingreturnsadjusted(WE+1:T)<VaRRiskSeeking250))/exp;                     % Violation Ratios = Observed no of siolations / Expected no of violations

figure;
set(groot,'DefaultAxesTickLabelInterpreter','Tex');
hold on
plot(xax2RiskSeeking,RiskSeekingreturnsadjusted(WE+1:T));
plot(xax2RiskSeeking,VaRRiskSeeking250,'LineWidth',1.25);
yline(0,'-.r','LineWidth',1.5);

% Plot configuration
ytickformat('%.2f');
ylim([-0.2 0.2]);
ylabel('Return','Interpreter','latex');
legend('Portfolio returns','VaR forecast','Location','northwest');
title('Backtest risk seeking portfolio','Interpreter','Latex','FontSize',12);
datetick('x','mmm-yyyy','keeplimits');
set(legend,'Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex','FontSize',12)

% Number of VaR violation
RiskSeekingComparison                                                       = cell2mat({PortReturn{3}(WE+1:T) cell2mat(cRiskSeeking)'});
ViolationRiskSeeking                                                        = cell(1,T);
for     x                                                                   = WE+1+250:T
        x1                                                                  = x-(WE+250);
        x2                                                                  = x-WE;
        ViolationRiskSeeking{x}                                             = sum(RiskSeekingComparison(x1:x2,1)<RiskSeekingComparison(x1:x2,2));
        ViolationRiskSeeking(:,WE+1:WE+250)                                 = {0};
end

% Plot of cumulative violations over the previous 250 days
figure;
plot(Dates(WE+1:T),cell2mat(ViolationRiskSeeking)','LineWidth',1.5)
xlim([730820 734137])
ytickformat('%.2f');
ylabel('Cumulative violations over previous 250 days','Interpreter','latex');
title('Risk seeking portfolio','Interpreter','Latex','FontSize',12);
datetick('x','mmm-yyyy');
set(gca,'TickLabelInterpreter','latex','FontSize',12)

% Determining multiplication factor "k"
kRiskSeeking                                                                = cell2mat(cell(1,T));
for     t                                                                   = WE+1:T
    if ViolationRiskSeeking{t} <= 4
        kRiskSeeking(t)  = 0;

    else
        if ViolationRiskSeeking{t} == 5
            kRiskSeeking(t)   = 0.40;

        else
            if ViolationRiskSeeking{t} == 6
                kRiskSeeking(t)   = 0.50;

            else
                if ViolationRiskSeeking{t} == 7
                    kRiskSeeking(t) = 0.65;

                else
                    if ViolationRiskSeeking{t} == 8
                        kRiskSeeking(t)   = 0.75;

                    else
                        if ViolationRiskSeeking{t} == 9
                            kRiskSeeking(t)   = 0.85;

                        else
                            if ViolationRiskSeeking{t} >= 10
                                kRiskSeeking(t)   = 1.00;  
                            end

                        end

                    end

                end

            end

        end

    end
end

% Average VaR over the last 60 days
AverageVaRRiskSeeking                                                       = cell(1,T);
for     y                                                                   = WE+1+59:T
        y1                                                                  = y-(WE+59);
        y2                                                                  = y-WE;
        AverageVaRRiskSeeking{y}                                            = mean(RiskSeekingComparison(y1:y2,2));
        AverageVaRRiskSeeking(:,WE+1:WE+59)                                 = {0};
        
end

AverageVaRRiskSeeking60                                                     = cellfun(@times, AverageVaRRiskSeeking, num2cell(3+kRiskSeeking), 'UniformOutput', false);

% Capital Charge 
previousdayVaRforecastRiskSeeking                                           = cell2mat(VaRforecastRiskSeeking');
previousdayVaRforecastRiskSeeking(1,:)                                      = [];
previousdayVaRforecastRiskSeeking(2370,1)                                   = 0;
CapitalChargeRiskSeeking                                                    = mean(abs(min(cell2mat({cell2mat(AverageVaRRiskSeeking60') previousdayVaRforecastRiskSeeking}),[],2)));

% Proportion out of green zone
RiskSeekingoutofGreenZone                                                   = (sum(kRiskSeeking(WE+1:T) > 0.4)/length(kRiskSeeking(WE+1:T)))*100;

% Print result 
formatSpec1                                                                 = 'Average Capital Charge is %5.4f%%, %.2f%% out of green zone.\n';
fprintf(formatSpec1,CapitalChargeRiskSeeking*100,RiskSeekingoutofGreenZone);

%% Summary of result
clc

summaryresult                                                               = table([RiskAversenov250;RiskNeutralnov250;RiskSeekingnov250],...,
                                                                                    [exp;exp;exp],...,
                                                                                    [RiskAverseoutofGreenZone;RiskNeutraloutofGreenZone;RiskSeekingoutofGreenZone],...,
                                                                                    [100*CapitalChargeRiskAverse;100*CapitalChargeRiskNeutral;100*CapitalChargeRiskSeeking],...,
                                                                                    'RowNames',{'Risk Averse','Risk Neutral','Risk Seeking'},...,
                                                                                    'VariableNames',{'No. of VaR violation','Expected no. of VaR violation','Proportion out of green zone (%)','Average Capital Charge (%)'});
                                       
                                       
                                       
disp(summaryresult);
%% Scenario testing 

adjustment                                                                  = 0.01;                                                                             % Shift up VaR
adjustment1                                                                 = -0.01;                                                                            % Shift down VaR

% Risk averse portfolio
% Backtest after the first 250 business days
% Plotting the stock returns and VaR forecast over the last 250 business days
cRiskAverseadjusted                                                         = cell2mat(VaRforecastRiskAverse)+adjustment;                                       
bRiskAverseadjusted                                                         = {(WE+1+1:T+1)' (cRiskAverseadjusted)'};
bRiskAverseadjusted{:,1}                                                    = (WE+1+1:T+1)';
xax2RiskAverseadjusted                                                      = Dates(bRiskAverseadjusted{:,1});                                                    
bRiskAverseadjusted{:,2}                                                    = (bRiskAverseadjusted{2});

cRiskAverseadjusted1                                                        = cell2mat(VaRforecastRiskAverse)+adjustment1;                                       
bRiskAverseadjusted1                                                        = {(WE+1+1:T+1)' (cRiskAverseadjusted1)'};
bRiskAverseadjusted1{:,1}                                                   = (WE+1+1:T+1)';
xax2RiskAverseadjusted1                                                     = Dates(bRiskAverseadjusted1{:,1});                                                  
bRiskAverseadjusted1{:,2}                                                   = (bRiskAverseadjusted1{2});

VaRRiskAverse250adjusted                                                    = bRiskAverseadjusted{:,2};
VaRRiskAverse250adjusted1                                                   = bRiskAverseadjusted1{:,2};

% Total number of VaR violation 
RiskAversenov250adjusted                                                    = length(find(RiskAversereturnsadjusted(WE+1:T)<VaRRiskAverse250adjusted));                         % Total number of VaR violation throughout VaR forecast window
RiskAverseloc_violation2_250adjusted                                        = find(RiskAversereturnsadjusted(WE+1:T)<VaRRiskAverse250adjusted);                                 % Locations where violation occur
RiskAverseVR250adjusted                                                     = length(find(RiskAversereturnsadjusted(WE+1:T)<VaRRiskAverse250adjusted))/exp;                     % Violation Ratios = Observed no of siolations / Expected no of violations

RiskAversenov250adjusted1                                                   = length(find(RiskAversereturnsadjusted(WE+1:T)<VaRRiskAverse250adjusted1));                        % Total number of VaR violation throughout VaR forecast window
RiskAverseloc_violation2_250adjusted1                                       = find(RiskAversereturnsadjusted(WE+1:T)<VaRRiskAverse250adjusted1);                                % Locations where violation occur
RiskAverseVR250adjusted1                                                    = length(find(RiskAversereturnsadjusted(WE+1:T)<VaRRiskAverse250adjusted1))/exp;                    % Violation Ratios = Observed no of siolations / Expected no of violations

figure;
hold on
plot(xax2RiskAverse,RiskAversereturnsadjusted(WE+1:T));
plot(xax2RiskAverse,VaRRiskAverse250,'LineWidth',1.25);
plot(xax2RiskAverseadjusted,VaRRiskAverse250adjusted,'LineWidth',1.25);
plot(xax2RiskAverseadjusted,VaRRiskAverse250adjusted1,'LineWidth',1.25,'color','c');
yline(0,'-.r','LineWidth',1.5);
set(groot,'DefaultAxesTickLabelInterpreter','Tex');

% Plot configuration
ytickformat('%.2f');
ylim([-0.2 0.2]);
ylabel('Return','Interpreter','latex');
legend('Portfolio returns','VaR forecast (original)','VaR forecast (shift up)','VaR forecast (shift down)','Location','northwest');
% legend('Portfolio returns','VaR forecast (original)','VaR forecast (shift up)','Location','northwest');
title('Backtest risk averse portfolio','Interpreter','Latex','FontSize',12);
datetick('x','mmm-yyyy','keeplimits');
set(legend,'Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex','FontSize',12)

% Number of VaR violation over the previous 250 days
RiskAverseadjustedComparison                                                = cell2mat({PortReturn{1}(WE+1:T) (cRiskAverseadjusted)'});
ViolationRiskAverseadjusted                                                 = cell(1,T);
for     x                                                                   = WE+1+250:T
        x1                                                                  = x-(WE+250);
        x2                                                                  = x-WE;
        ViolationRiskAverseadjusted{x}                                      = sum(RiskAverseadjustedComparison(x1:x2,1)<RiskAverseadjustedComparison(x1:x2,2));
        ViolationRiskAverseadjusted(:,WE+1:WE+250)                          = {0};
end

% Determining multiplication factor "k"
kRiskAverseadjusted                                                         = cell2mat(cell(1,T));
for     t                                                                   = WE+1:T
    if ViolationRiskAverseadjusted{t} <= 4
        kRiskAverseadjusted(t)  = 0;

    else
        if ViolationRiskAverseadjusted{t} == 5
            kRiskAverseadjusted(t)   = 0.40;

        else
            if ViolationRiskAverseadjusted{t} == 6
                kRiskAverseadjusted(t)   = 0.50;

            else
                if ViolationRiskAverseadjusted{t} == 7
                    kRiskAverseadjusted(t) = 0.65;

                else
                    if ViolationRiskAverseadjusted{t} == 8
                        kRiskAverseadjusted(t)   = 0.75;

                    else
                        if ViolationRiskAverseadjusted{t} == 9
                            kRiskAverseadjusted(t)   = 0.85;

                        else
                            if ViolationRiskAverseadjusted{t} >= 10
                                kRiskAverseadjusted(t)   = 1.00;  
                            end

                        end

                    end

                end

            end

        end

    end
end

% Average VaR over the last 60 days
AverageVaRRiskAverseadjusted                                                = cell(1,T);
for     y                                                                   = WE+1+59:T
        y1                                                                  = y-(WE+59);
        y2                                                                  = y-WE;
        AverageVaRRiskAverseadjusted{y}                                     = mean(RiskAverseadjustedComparison(y1:y2,2));
        AverageVaRRiskAverseadjusted(:,WE+1:WE+59)                          = {0};
        
end

AverageVaRRiskAverseadjusted60                                              = cellfun(@times, AverageVaRRiskAverseadjusted, num2cell(3+kRiskAverseadjusted), 'UniformOutput', false);

% Capital Charge 
previousdayVaRforecastRiskAverseadjusted                                    = cRiskAverseadjusted';
previousdayVaRforecastRiskAverseadjusted(1,:)                               = [];
previousdayVaRforecastRiskAverseadjusted(2370,1)                            = 0;
CapitalChargeRiskAverseadjusted                                             = mean(abs(min(cell2mat({cell2mat(AverageVaRRiskAverseadjusted60') previousdayVaRforecastRiskAverseadjusted}),[],2)));

% Proportion out of green zone
RiskAverseadjustedoutofGreenZone                                            = (sum(kRiskAverseadjusted(WE+1:T) > 0.4)/length(kRiskAverseadjusted(WE+1:T)))*100;

% Print result 
formatSpec1                                                                 = 'Average Capital Charge is %5.4f%%, %.2f%% out of green zone.\n';
fprintf(formatSpec1,CapitalChargeRiskAverseadjusted*100,RiskAverseadjustedoutofGreenZone);

% Number of VaR violation
RiskAverseadjusted1Comparison                                               = cell2mat({PortReturn{1}(WE+1:T) (cRiskAverseadjusted1)'});
ViolationRiskAverseadjusted1                                                = cell(1,T);
for     x                                                                   = WE+1+250:T
        x1                                                                  = x-(WE+250);
        x2                                                                  = x-WE;
        ViolationRiskAverseadjusted1{x}                                     = sum(RiskAverseadjusted1Comparison(x1:x2,1)<RiskAverseadjusted1Comparison(x1:x2,2));
        ViolationRiskAverseadjusted1(:,WE+1:WE+250)                         = {0};
end

% Determining multiplication factor "k"
kRiskAverseadjusted1                                                        = cell2mat(cell(1,T));
for     t                                                                   = WE+1:T
    if ViolationRiskAverseadjusted1{t} <= 4
        kRiskAverseadjusted1(t)  = 0;

    else
        if ViolationRiskAverseadjusted1{t} == 5
            kRiskAverseadjusted1(t)   = 0.40;

        else
            if ViolationRiskAverseadjusted1{t} == 6
                kRiskAverseadjusted1(t)   = 0.50;

            else
                if ViolationRiskAverseadjusted1{t} == 7
                    kRiskAverseadjusted1(t) = 0.65;

                else
                    if ViolationRiskAverseadjusted1{t} == 8
                        kRiskAverseadjusted1(t)   = 0.75;

                    else
                        if ViolationRiskAverseadjusted1{t} == 9
                            kRiskAverseadjusted1(t)   = 0.85;

                        else
                            if ViolationRiskAverseadjusted1{t} >= 10
                                kRiskAverseadjusted1(t)   = 1.00;  
                            end

                        end

                    end

                end

            end

        end

    end
end

% Average VaR over the last 60 days
AverageVaRRiskAverseadjusted1                                               = cell(1,T);
for     y                                                                   = WE+1+59:T
        y1                                                                  = y-(WE+59);
        y2                                                                  = y-WE;
        AverageVaRRiskAverseadjusted1{y}                                    = mean(RiskAverseadjusted1Comparison(y1:y2,2));
        AverageVaRRiskAverseadjusted1(:,WE+1:WE+59)                         = {0};
        
end

AverageVaRRiskAverseadjusted160                                             = cellfun(@times, AverageVaRRiskAverseadjusted1, num2cell(3+kRiskAverseadjusted1), 'UniformOutput', false);

% Capital Charge 
previousdayVaRforecastRiskAverseadjusted1                                   = cRiskAverseadjusted1';
previousdayVaRforecastRiskAverseadjusted1(1,:)                              = [];
previousdayVaRforecastRiskAverseadjusted1(2370,1)                           = 0;
CapitalChargeRiskAverseadjusted1                                            = mean(abs(min(cell2mat({cell2mat(AverageVaRRiskAverseadjusted160') previousdayVaRforecastRiskAverseadjusted1}),[],2)));

% Proportion out of green zone
RiskAverseadjusted1outofGreenZone                                           = (sum(kRiskAverseadjusted1(WE+1:T) > 0.4)/length(kRiskAverseadjusted1(WE+1:T)))*100;

% Print result 
formatSpec1                                                                 = 'Average Capital Charge is %5.4f%%, %.2f%% out of green zone.\n';
fprintf(formatSpec1,CapitalChargeRiskAverseadjusted1*100,RiskAverseadjusted1outofGreenZone);


% Risk neutral portfolio
% Backtest after the first 250 business days
% Plotting the stock returns and VaR forecast over the last 250 business days
cRiskNeutraladjusted                                                        = cell2mat(VaRforecastRiskNeutral)+adjustment;                                       
bRiskNeutraladjusted                                                        = {(WE+1+1:T+1)' (cRiskNeutraladjusted)'};
bRiskNeutraladjusted{:,1}                                                   = (WE+1+1:T+1)';
xax2RiskNeutraladjusted                                                     = Dates(bRiskNeutraladjusted{:,1});                                                   
bRiskNeutraladjusted{:,2}                                                   = (bRiskNeutraladjusted{2});

cRiskNeutraladjusted1                                                       = cell2mat(VaRforecastRiskNeutral)+adjustment1;                                       
bRiskNeutraladjusted1                                                       = {(WE+1+1:T+1)' (cRiskNeutraladjusted1)'};
bRiskNeutraladjusted1{:,1}                                                  = (WE+1+1:T+1)';
xax2RiskNeutraladjusted1                                                    = Dates(bRiskNeutraladjusted1{:,1});                                                  
bRiskNeutraladjusted1{:,2}                                                  = (bRiskNeutraladjusted1{2});

VaRRiskNeutral250adjusted                                                   = bRiskNeutraladjusted{:,2};
VaRRiskNeutral250adjusted1                                                  = bRiskNeutraladjusted1{:,2};

% Total number of VaR violation
RiskNeutralnov250adjusted                                                   = length(find(RiskNeutralreturnsadjusted(WE+1:T)<VaRRiskNeutral250adjusted));                         % Total number of VaR violation throughout VaR forecast window
RiskNeutralloc_violation2_250adjusted                                       = find(RiskNeutralreturnsadjusted(WE+1:T)<VaRRiskNeutral250adjusted);                                 % Locations where violation occur
RiskNeutralVR250adjusted                                                    = length(find(RiskNeutralreturnsadjusted(WE+1:T)<VaRRiskNeutral250adjusted))/exp;                     % Violation Ratios = Observed no of siolations / Expected no of violations
        
RiskNeutralnov250adjusted1                                                  = length(find(RiskNeutralreturnsadjusted(WE+1:T)<VaRRiskNeutral250adjusted1));                        % Total number of VaR violation throughout VaR forecast window
RiskNeutralloc_violation2_250adjusted1                                      = find(RiskNeutralreturnsadjusted(WE+1:T)<VaRRiskNeutral250adjusted1);                                % Locations where violation occur
RiskNeutralVR250adjusted1                                                   = length(find(RiskNeutralreturnsadjusted(WE+1:T)<VaRRiskNeutral250adjusted1))/exp;                    % Violation Ratios = Observed no of siolations / Expected no of violations

figure;
hold on
plot(xax2RiskNeutral,RiskNeutralreturnsadjusted(WE+1:T));
plot(xax2RiskNeutral,VaRRiskNeutral250,'LineWidth',1.25);
plot(xax2RiskNeutraladjusted,VaRRiskNeutral250adjusted,'LineWidth',1.25);
plot(xax2RiskNeutraladjusted,VaRRiskNeutral250adjusted1,'LineWidth',1.25,'color','c');
yline(0,'-.r','LineWidth',1.5);
set(groot,'DefaultAxesTickLabelInterpreter','Tex');

% Plot configuration
ytickformat('%.2f');
ylim([-0.2 0.2]);
ylabel('Return','Interpreter','latex');
legend('Portfolio returns','VaR forecast (original)','VaR forecast (shift up)','VaR forecast (shift down)','Location','northwest');
% legend('Portfolio returns','VaR forecast (original)','VaR forecast (shift up)','Location','northwest');
title('Backtest risk neutral portfolio','Interpreter','Latex','FontSize',12);
datetick('x','mmm-yyyy','keeplimits');
set(legend,'Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex','FontSize',12)

% Number of VaR violation over the previous 250 days
RiskNeutraladjustedComparison                                               = cell2mat({PortReturn{2}(WE+1:T) (cRiskNeutraladjusted)'});
    ViolationRiskNeutraladjusted                                            = cell(1,T);
for     x                                                                   = WE+1+250:T
        x1                                                                  = x-(WE+250);
        x2                                                                  = x-WE;
        ViolationRiskNeutraladjusted{x}                                     = sum(RiskNeutraladjustedComparison(x1:x2,1)<RiskNeutraladjustedComparison(x1:x2,2));
        ViolationRiskNeutraladjusted(:,WE+1:WE+250)                         = {0};
end

% Determining multiplication factor "k"
kRiskNeutraladjusted                                                        = cell2mat(cell(1,T));
for     t                                                                   = WE+1:T
    if ViolationRiskNeutraladjusted{t} <= 4
        kRiskNeutraladjusted(t)  = 0;

    else
        if ViolationRiskNeutraladjusted{t} == 5
            kRiskNeutraladjusted(t)   = 0.40;

        else
            if ViolationRiskNeutraladjusted{t} == 6
                kRiskNeutraladjusted(t)   = 0.50;

            else
                if ViolationRiskNeutraladjusted{t} == 7
                    kRiskNeutraladjusted(t) = 0.65;

                else
                    if ViolationRiskNeutraladjusted{t} == 8
                        kRiskNeutraladjusted(t)   = 0.75;

                    else
                        if ViolationRiskNeutraladjusted{t} == 9
                            kRiskNeutraladjusted(t)   = 0.85;

                        else
                            if ViolationRiskNeutraladjusted{t} >= 10
                                kRiskNeutraladjusted(t)   = 1.00;  
                            end

                        end

                    end

                end

            end

        end

    end
end

% Average VaR over the last 60 days
AverageVaRRiskNeutraladjusted = cell(1,T);
for     y                                                                   = WE+1+59:T
        y1                                                                  = y-(WE+59);
        y2                                                                  = y-WE;
        AverageVaRRiskNeutraladjusted{y}                                    = mean(RiskNeutraladjustedComparison(y1:y2,2));
        AverageVaRRiskNeutraladjusted(:,WE+1:WE+59)                         = {0};
        
end

AverageVaRRiskNeutraladjusted60                                             = cellfun(@times, AverageVaRRiskNeutraladjusted, num2cell(3+kRiskNeutraladjusted), 'UniformOutput', false);

% Capital Charge 
previousdayVaRforecastRiskNeutraladjusted                                   = cRiskNeutraladjusted';
previousdayVaRforecastRiskNeutraladjusted(1,:)                              = [];
previousdayVaRforecastRiskNeutraladjusted(2370,1)                           = 0;
CapitalChargeRiskNeutraladjusted                                            = mean(abs(min(cell2mat({cell2mat(AverageVaRRiskNeutraladjusted60') previousdayVaRforecastRiskNeutraladjusted}),[],2)));

% Proportion out of green zone
RiskNeutraladjustedoutofGreenZone                                           = (sum(kRiskNeutraladjusted(WE+1:T) > 0.4)/length(kRiskNeutraladjusted(WE+1:T)))*100;

% Print result 
formatSpec1                                                                 = 'Average Capital Charge is %5.4f%%, %.2f%% out of green zone.\n';
fprintf(formatSpec1,CapitalChargeRiskNeutraladjusted*100,RiskNeutraladjustedoutofGreenZone);

% Number of VaR violation
RiskNeutraladjusted1Comparison                                              = cell2mat({PortReturn{2}(WE+1:T) (cRiskNeutraladjusted1)'});
ViolationRiskNeutraladjusted1                                               = cell(1,T);
for     x                                                                   = WE+1+250:T
        x1                                                                  = x-(WE+250);
        x2                                                                  = x-WE;
        ViolationRiskNeutraladjusted1{x}                                    = sum(RiskNeutraladjusted1Comparison(x1:x2,1)<RiskNeutraladjusted1Comparison(x1:x2,2));
        ViolationRiskNeutraladjusted1(:,WE+1:WE+250)                        = {0};
end

% Determining multiplication factor "k"
kRiskNeutraladjusted1                                                       = cell2mat(cell(1,T));
for     t                                                                   = WE+1:T
    if ViolationRiskNeutraladjusted1{t} <= 4
        kRiskNeutraladjusted1(t)  = 0;

    else
        if ViolationRiskNeutraladjusted1{t} == 5
            kRiskNeutraladjusted1(t)   = 0.40;

        else
            if ViolationRiskNeutraladjusted1{t} == 6
                kRiskNeutraladjusted1(t)   = 0.50;

            else
                if ViolationRiskNeutraladjusted1{t} == 7
                    kRiskNeutraladjusted1(t) = 0.65;

                else
                    if ViolationRiskNeutraladjusted1{t} == 8
                        kRiskNeutraladjusted1(t)   = 0.75;

                    else
                        if ViolationRiskNeutraladjusted1{t} == 9
                            kRiskNeutraladjusted1(t)   = 0.85;

                        else
                            if ViolationRiskNeutraladjusted1{t} >= 10
                                kRiskNeutraladjusted1(t)   = 1.00;  
                            end

                        end

                    end

                end

            end

        end

    end
end

% Average VaR over the last 60 days
AverageVaRRiskNeutraladjusted1                                              = cell(1,T);
for     y                                                                   = WE+1+59:T
        y1                                                                  = y-(WE+59);
        y2                                                                  = y-WE;
        AverageVaRRiskNeutraladjusted1{y}                                   = mean(RiskNeutraladjusted1Comparison(y1:y2,2));
        AverageVaRRiskNeutraladjusted1(:,WE+1:WE+59)                        = {0};
        
end

AverageVaRRiskNeutraladjusted160                                            = cellfun(@times, AverageVaRRiskNeutraladjusted1, num2cell(3+kRiskNeutraladjusted1), 'UniformOutput', false);

% Capital Charge 
previousdayVaRforecastRiskNeutraladjusted1                                  = cRiskNeutraladjusted1';
previousdayVaRforecastRiskNeutraladjusted1(1,:)                             = [];
previousdayVaRforecastRiskNeutraladjusted1(2370,1)                          = 0;
CapitalChargeRiskNeutraladjusted1                                           = mean(abs(min(cell2mat({cell2mat(AverageVaRRiskNeutraladjusted160') previousdayVaRforecastRiskNeutraladjusted1}),[],2)));

% Proportion out of green zone 
RiskNeutraladjusted1outofGreenZone                                          = (sum(kRiskNeutraladjusted1(WE+1:T) > 0.4)/length(kRiskNeutraladjusted1(WE+1:T)))*100;

% Print result 
formatSpec1                                                                 = 'Average Capital Charge is %5.4f%%, %.2f%% out of green zone.\n';
fprintf(formatSpec1,CapitalChargeRiskNeutraladjusted1*100,RiskNeutraladjusted1outofGreenZone);


% Risk seeking portfolio
% Backtest after the first 250 business days
% Plotting the stock returns and VaR forecast over the last 250 business days
cRiskSeekingadjusted                                                        = cell2mat(VaRforecastRiskSeeking)+adjustment;                                       
bRiskSeekingadjusted                                                        = {(WE+1+1:T+1)' (cRiskSeekingadjusted)'};
bRiskSeekingadjusted{:,1}                                                   = (WE+1+1:T+1)';
xax2RiskSeekingadjusted                                                     = Dates(bRiskSeekingadjusted{:,1});                                                   
bRiskSeekingadjusted{:,2}                                                   = (bRiskSeekingadjusted{2});

cRiskSeekingadjusted1                                                       = cell2mat(VaRforecastRiskSeeking)+adjustment1;                                       
bRiskSeekingadjusted1                                                       = {(WE+1+1:T+1)' (cRiskSeekingadjusted1)'};
bRiskSeekingadjusted1{:,1}                                                  = (WE+1+1:T+1)';
xax2RiskSeekingadjusted1                                                    = Dates(bRiskSeekingadjusted1{:,1});                                                   
bRiskSeekingadjusted1{:,2}                                                  = (bRiskSeekingadjusted1{2});

VaRRiskSeeking250adjusted                                                   = bRiskSeekingadjusted{:,2};
VaRRiskSeeking250adjusted1                                                  = bRiskSeekingadjusted1{:,2};

% Total number of VaR violation
RiskSeekingnov250adjusted                                                   = length(find(RiskSeekingreturnsadjusted(WE+1:T)<VaRRiskSeeking250adjusted));                         % Number of Violation
RiskSeekingloc_violation2_250adjusted                                       = find(RiskSeekingreturnsadjusted(WE+1:T)<VaRRiskSeeking250adjusted);                                 % Locations where violation occur
RiskSeekingVR250adjusted                                                    = length(find(RiskSeekingreturnsadjusted(WE+1:T)<VaRRiskSeeking250adjusted))/exp;       % Violation Ratios = Observed no of siolations / Expected no of violations

RiskSeekingnov250adjusted1                                                  = length(find(RiskSeekingreturnsadjusted(WE+1:T)<VaRRiskSeeking250adjusted1));                        % Number of Violation
RiskSeekingloc_violation2_250adjusted1                                      = find(RiskSeekingreturnsadjusted(WE+1:T)<VaRRiskSeeking250adjusted1);                                % Locations where violation occur
RiskSeekingVR250adjusted1                                                   = length(find(RiskSeekingreturnsadjusted(WE+1:T)<VaRRiskSeeking250adjusted1))/exp;      % Violation Ratios = Observed no of siolations / Expected no of violations

figure;
hold on
plot(xax2RiskSeeking,RiskSeekingreturnsadjusted(WE+1:T));
plot(xax2RiskSeeking,VaRRiskSeeking250,'LineWidth',1.25);
plot(xax2RiskSeekingadjusted,VaRRiskSeeking250adjusted,'LineWidth',1.25);
plot(xax2RiskSeekingadjusted,VaRRiskSeeking250adjusted1,'LineWidth',1.25,'color','c');
yline(0,'-.r','LineWidth',1.5);
set(groot,'DefaultAxesTickLabelInterpreter','Tex');

% Plot configuration
ytickformat('%.2f');
ylim([-0.2 0.2]);
ylabel('Return','Interpreter','latex');
legend('Portfolio returns','VaR forecast (original)','VaR forecast (shift up)','VaR forecast (shift down)','Location','northwest');legend('Portfolio returns','VaR forecast (original)','VaR forecast (shift up)','VaR forecast (shift down)','Location','northwest');
% legend('Portfolio returns','VaR forecast (original)','VaR forecast (shift up)','Location','northwest');
title('Backtest risk seeking portfolio','Interpreter','Latex','FontSize',12);
datetick('x','mmm-yyyy','keeplimits');
set(legend,'Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex','FontSize',12)

% Number of VaR violation over the previous 250 days
RiskSeekingadjustedComparison                                               = cell2mat({PortReturn{3}(WE+1:T) (cRiskSeekingadjusted)'});
ViolationRiskSeekingadjusted                                                = cell(1,T);
for     x                                                                   = WE+1+250:T
        x1                                                                  = x-(WE+250);
        x2                                                                  = x-WE;
        ViolationRiskSeekingadjusted{x}                                     = sum(RiskSeekingadjustedComparison(x1:x2,1)<RiskSeekingadjustedComparison(x1:x2,2));
        ViolationRiskSeekingadjusted(:,WE+1:WE+250)                         = {0};
end

% Determining multiplication factor "k"
kRiskSeekingadjusted                                                        = cell2mat(cell(1,T));
for     t                                                                   = WE+1:T
    if ViolationRiskSeekingadjusted{t} <= 4 
        kRiskSeekingadjusted(t)  = 0;

    else
        if ViolationRiskSeekingadjusted{t} == 5
            kRiskSeekingadjusted(t)   = 0.40;

        else
            if ViolationRiskSeekingadjusted{t} == 6
                kRiskSeekingadjusted(t)   = 0.50;

            else
                if ViolationRiskSeekingadjusted{t} == 7
                    kRiskSeekingadjusted(t) = 0.65;

                else
                    if ViolationRiskSeekingadjusted{t} == 8
                        kRiskSeekingadjusted(t)   = 0.75;

                    else
                        if ViolationRiskSeekingadjusted{t} == 9
                            kRiskSeekingadjusted(t)   = 0.85;

                        else
                            if ViolationRiskSeekingadjusted{t} >= 10
                                kRiskSeekingadjusted(t)   = 1.00;  
                            end

                        end

                    end

                end

            end

        end

    end
end

% Average VaR over the last 60 days
AverageVaRRiskSeekingadjusted                                               = cell(1,T);
for     y                                                                   = WE+1+59:T
        y1                                                                  = y-(WE+59);
        y2                                                                  = y-WE;
        AverageVaRRiskSeekingadjusted{y}                                    = mean(RiskSeekingadjustedComparison(y1:y2,2));
        AverageVaRRiskSeekingadjusted(:,WE+1:WE+59)                         = {0};
        
end

AverageVaRRiskSeekingadjusted60                                             = cellfun(@times, AverageVaRRiskSeekingadjusted, num2cell(3+kRiskSeekingadjusted), 'UniformOutput', false);

% Capital Charge 
previousdayVaRforecastRiskSeekingadjusted                                   = cRiskSeekingadjusted';
previousdayVaRforecastRiskSeekingadjusted(1,:)                              = [];
previousdayVaRforecastRiskSeekingadjusted(2370,1)                           = 0;
CapitalChargeRiskSeekingadjusted                                            = mean(abs(min(cell2mat({cell2mat(AverageVaRRiskSeekingadjusted60') previousdayVaRforecastRiskSeekingadjusted}),[],2)));

% Proportion out of green zone
RiskSeekingadjustedoutofGreenZone                                           = (sum(kRiskSeekingadjusted(WE+1:T) > 0.4)/length(kRiskSeekingadjusted(WE+1:T)))*100;

% Print result 
formatSpec1                                                                 = 'Average Capital Charge is %5.4f%%, %.2f%% out of green zone.\n';
fprintf(formatSpec1,CapitalChargeRiskSeekingadjusted*100,RiskSeekingadjustedoutofGreenZone);

% Number of VaR violation over the previous 250 days
RiskSeekingadjusted1Comparison                                              = cell2mat({PortReturn{3}(WE+1:T) (cRiskSeekingadjusted1)'});
ViolationRiskSeekingadjusted1                                               = cell(1,T);
for     x                                                                   = WE+1+250:T
        x1                                                                  = x-(WE+250);
        x2                                                                  = x-WE;
        ViolationRiskSeekingadjusted1{x}                                    = sum(RiskSeekingadjusted1Comparison(x1:x2,1)<RiskSeekingadjusted1Comparison(x1:x2,2));
        ViolationRiskSeekingadjusted1(:,WE+1:WE+250)                        = {0};
end

% Determining multiplication factor "k"
kRiskSeekingadjusted1                                                       = cell2mat(cell(1,T));
for     t                                                                   = WE+1:T
    if ViolationRiskSeekingadjusted1{t} <= 4
        kRiskSeekingadjusted1(t)  = 0;

    else
        if ViolationRiskSeekingadjusted1{t} == 5
            kRiskSeekingadjusted1(t)   = 0.40;

        else
            if ViolationRiskSeekingadjusted1{t} == 6
                kRiskSeekingadjusted1(t)   = 0.50;

            else
                if ViolationRiskSeekingadjusted1{t} == 7
                    kRiskSeekingadjusted1(t) = 0.65;

                else
                    if ViolationRiskSeekingadjusted1{t} == 8
                        kRiskSeekingadjusted1(t)   = 0.75;

                    else
                        if ViolationRiskSeekingadjusted1{t} == 9
                            kRiskSeekingadjusted1(t)   = 0.85;

                        else
                            if ViolationRiskSeekingadjusted1{t} >= 10
                                kRiskSeekingadjusted1(t)   = 1.00;  
                            end

                        end

                    end

                end

            end

        end

    end
end

% Average VaR over the last 60 days
AverageVaRRiskSeekingadjusted1                                              = cell(1,T);
for     y                                                                   = WE+1+59:T
        y1                                                                  = y-(WE+59);
        y2                                                                  = y-WE;
        AverageVaRRiskSeekingadjusted1{y}                                   = mean(RiskSeekingadjusted1Comparison(y1:y2,2));
        AverageVaRRiskSeekingadjusted1(:,WE+1:WE+59)                        = {0};       
end

AverageVaRRiskSeekingadjusted160                                            = cellfun(@times, AverageVaRRiskSeekingadjusted1, num2cell(3+kRiskSeekingadjusted1), 'UniformOutput', false);

% Capital Charge 
previousdayVaRforecastRiskSeekingadjusted1                                  = cRiskSeekingadjusted1';
previousdayVaRforecastRiskSeekingadjusted1(1,:)                             = [];
previousdayVaRforecastRiskSeekingadjusted1(2370,1)                          = 0;
CapitalChargeRiskSeekingadjusted1                                           = mean(abs(min(cell2mat({cell2mat(AverageVaRRiskSeekingadjusted160') previousdayVaRforecastRiskSeekingadjusted1}),[],2)));

% Proportion out of green zone
RiskSeekingadjusted1outofGreenZone                                          = (sum(kRiskSeekingadjusted1(WE+1:T) > 0.4)/length(kRiskSeekingadjusted1(WE+1:T)))*100;

% Print result 
formatSpec1                                                                 = 'Average Capital Charge is %5.4f%%, %.2f%% out of green zone.\n';
fprintf(formatSpec1,CapitalChargeRiskSeekingadjusted1*100,RiskSeekingadjusted1outofGreenZone);



summaryresultadjusted1                                                      = table([RiskAversenov250adjusted;RiskAversenov250;RiskAversenov250adjusted1;RiskNeutralnov250adjusted;RiskNeutralnov250;RiskNeutralnov250adjusted1;RiskSeekingnov250adjusted;RiskSeekingnov250;RiskSeekingnov250adjusted1],...,
                                                                                    [exp;exp;exp;exp;exp;exp;exp;exp;exp],...,
                                                                                    [RiskAverseadjustedoutofGreenZone;RiskAverseoutofGreenZone;RiskAverseadjusted1outofGreenZone;RiskNeutraladjustedoutofGreenZone;RiskNeutraloutofGreenZone;RiskNeutraladjusted1outofGreenZone;RiskSeekingadjustedoutofGreenZone;RiskSeekingoutofGreenZone;RiskSeekingadjusted1outofGreenZone],...,
                                                                                    [100*CapitalChargeRiskAverseadjusted;100*CapitalChargeRiskAverse;100*CapitalChargeRiskAverseadjusted1;100*CapitalChargeRiskNeutraladjusted;100*CapitalChargeRiskNeutral;100*CapitalChargeRiskNeutraladjusted1;100*CapitalChargeRiskSeekingadjusted;100*CapitalChargeRiskSeeking;100*CapitalChargeRiskSeekingadjusted1],...,
                                                                                    'RowNames',{'Risk averse (shift up)','Risk averse (original)','Risk averse (shift down)','Risk neutral (shift up)','Risk neutral (original)','Risk neutral (shift down)','Risk seeking (shift up)','Risk seeking (original)','Risk seeking (shift down)'},...,
                                                                                    'VariableNames',{'Observed no. of violations','Expected no. of violations','Proportion out of green zone (%)','Capital Charge (%)'});

summaryresultadjusted2                                                      = table([RiskAversenov250adjusted;RiskAversenov250;RiskNeutralnov250adjusted;RiskNeutralnov250;RiskSeekingnov250adjusted;RiskSeekingnov250],...,
                                                                                   [exp;exp;exp;exp;exp;exp],...,
                                                                                   [RiskAverseadjustedoutofGreenZone;RiskAverseoutofGreenZone;RiskNeutraladjustedoutofGreenZone;RiskNeutraloutofGreenZone;RiskSeekingadjustedoutofGreenZone;RiskSeekingoutofGreenZone],...,
                                                                                   [100*CapitalChargeRiskAverseadjusted;100*CapitalChargeRiskAverse;100*CapitalChargeRiskNeutraladjusted;100*CapitalChargeRiskNeutral;100*CapitalChargeRiskSeekingadjusted;100*CapitalChargeRiskSeeking],...,
                                                                                   'RowNames',{'Risk averse (shift up)','Risk averse (original)','Risk neutral (shift up)','Risk neutral (original)','Risk seeking (shift up)','Risk seeking (original)'},...,
                                                                                   'VariableNames',{'Observed no. of violations','Expected no. of violations','Proportion out of green zone (%)','Capital Charge (%)'});
                                            
                                            
% print -depsc epsFig change epsFig to the desired name

%% Summary of Portfolio weight and Back test result
clc

% Display the weight of all portfolios
fprintf('\n\nPortfolio with High Risk Firms\n\n');
disp(RiskAversePort);
fprintf('\n\nPortfolio with Middle Risk Firms\n\n');
disp(RiskNeutralPort);
fprintf('\n\nPortfolio with Low Risk Firms\n\n');
disp(RiskSeekingPort);
  
% Display the back test result of all portfolios (adjusted and non adjusted)
fprintf('\n\nSummary of back test result\n\n');
disp(summaryresultadjusted1);