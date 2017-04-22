function [netChoice,hyperparams]= chooseNetwork( X, t, pathtohyperparameters)

global 	CHOSENHYPERPARAMS
CHOSENHYPERPARAMS =struct();

% trainingFunctions={'trainoss'};%,'traingdx','traingda','traingdm'};%, 'trainlm','traincgb','trainbfg','trainscg'}; %%Levenberg-Marquardt backpropagation (trainlm) works best
% errorFunctions={'crossentropy'};%,'sae','mae','mse','sse'}; %%Sum absolute error (sae) performance function works best
% transferFunctions={'radbas'}; %,'tansig','logsig'
% lr = maxlinlr(X');%%Maximum learning rate

%%Read hyperparameters
[~,~,hp]=xlsread(pathtohyperparameters); hp{1,1}='Parameter';
hp=cell2table(hp, 'VariableNames',hp(1,:),'Rownames',hp(:,1)); hp(1,:)=[]; hp(:,'Parameter')=[];

nNeurons=getRandomHyperWithinRange('Hiddenlayersize',hp);

netChoice =deepNet(X,t,nNeurons,hp);
netChoice.trainParam.sigma=getRandomHyperWithinRange('LearninRatesigma',hp);
netChoice.trainParam.lambda=getRandomHyperWithinRange('LearninRatelambda',hp);
netChoice.inputs{1}.processFcns=getRandomHyperOption('Inputpreprocessing',hp);%%Preprocessing functions

%%Adapt network
% netChoice.layerWeights{1,1}.learnParam.lr=lr;
% netChoice.layerWeights{1,1}.learnParam.mc=mc;
disp(['Training network with --> Training function: ',netChoice.trainFcn,'; nNeurons: ', num2str(netChoice.layers{1}.size),'; Performance function: ', netChoice.performFcn,'; Transfer function: ', netChoice.layers{1}.transferFcn]);
disp(['Adapt function: ',netChoice.adaptFcn])
netChoice= train(netChoice,X',t');

if length(netChoice.layers)>2 %%If this is a deep network make further specs to learning functions
    netChoice.layerWeights{2,1}.learnFcn='learngdm';
    netChoice.layerWeights{3,2}.learnFcn='learngdm';
    netChoice.adaptFcn='adaptwb';
end
disp('Choosen network characteristics:')
disp(['Training function: ',netChoice.trainFcn,'; nNeurons: ', num2str(netChoice.layers{1}.size),'; Performance function: ', netChoice.performFcn,'; Transfer function: ', netChoice.layers{1}.transferFcn]);
disp(['Adapt function: ',netChoice.adaptFcn])

getRandomHyperWithinRange('Numberofadaptations',hp);
hyperparams=CHOSENHYPERPARAMS;
end


function deepnet =deepNet(X,t, neurons,hp)

%%get hyperparameters to be set:
l2wr=getRandomHyperWithinRange('L2WeightRegularization',hp);
sr=getRandomHyperWithinRange('SparsityRegularization',hp);
sp=getRandomHyperWithinRange('SparsityProportion',hp);
me=getRandomHyperWithinRange('MaxEpochs',hp);

autoenc1 = trainAutoencoder(X',neurons, ...
    'MaxEpochs',me, ...
    'L2WeightRegularization',l2wr, ...
    'SparsityRegularization',sr, ...
    'SparsityProportion',sp, ...
    'ScaleData', false);

feat1 = encode(autoenc1,X');

autoenc2 = trainAutoencoder(feat1,ceil(neurons/2), ...
    'MaxEpochs',round(me/4), ...
    'L2WeightRegularization',round(l2wr/2), ...
    'SparsityRegularization',sr, ...
    'SparsityProportion',sp, ...
    'ScaleData', false);

feat2 = encode(autoenc2,feat1);

softnet = trainSoftmaxLayer(feat2,t','MaxEpochs',me);

deepnet = stack(autoenc1,autoenc2,softnet);


end

function x=getRandomHyperWithinRange(hyperName,hp)
%%helper function for random access of a given hyperparameter within given
%%range
global 	CHOSENHYPERPARAMS

x=table2array(hp{hyperName,'Default'}); %% Default value
if  table2array(hp{hyperName,'Set'}) %%set only if this hyperparameter is to be set
    maxV=table2array(hp{hyperName,'max'});  minV=table2array(hp{hyperName,'min'});
    range=maxV-minV;
    magnitude=max(1,1000/range);
    x=randi(round(magnitude* range),1,1)/magnitude+minV;
    if mod(maxV,1)==0 && mod(minV,1)==0
        x=round(x);
    end
end
CHOSENHYPERPARAMS=setfield(CHOSENHYPERPARAMS,char(hyperName),x);
end

function x=getRandomHyperOption(hyperName,hp)
%%helper function for random choice of a given hyperparameter from given
%%list of options
global 	CHOSENHYPERPARAMS

x=(hp{hyperName,'Default'}); %% Default value
if  strcmp(x{:},'NaN')
    x={};
end
if  table2array(hp{hyperName,'Set'}) %%set only if this hyperparameter is to be set
    options=strsplit(char(table2cell(hp(hyperName,'Options'))),', ');
    x=options(randi(length(options),1,1));
end
CHOSENHYPERPARAMS=setfield(CHOSENHYPERPARAMS,char(hyperName),x);
end
