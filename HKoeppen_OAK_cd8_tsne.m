% Code to create tSNE plots, train, and apply SVM classification on
% readouts from immunophenotyping algorithms for CD8 density in CK pos and
% CK neg regions

%%  import POPLAR data to train, manual phenotype 1 = desert, 2 = excluded, 3 = desert

man_pheno = Poplar_data(:,2);
meas = cell2mat(Poplar_data(:,3:22));

for ii = 1 : height(Poplar_data)
    man_label(ii,1) = cell2mat(man_pheno(ii));
end





%%
%create TSNE plot

rng('default') % for reproducibility
Y = tsne(meas,'Algorithm','exact','Distance','mahalanobis');
subplot(2,2,1)
gscatter(Y(:,1),Y(:,2),man_label)
title('Mahalanobis')

rng('default') % for fair comparison
Y = tsne(meas,'Algorithm','exact','Distance','cosine');
subplot(2,2,2)
gscatter(Y(:,1),Y(:,2),man_label)
title('Cosine')

rng('default') % for fair comparison
Y = tsne(meas,'Algorithm','exact','Distance','chebychev');
subplot(2,2,3)
gscatter(Y(:,1),Y(:,2),man_label)
title('Chebychev')

rng('default') % for fair comparison
Y = tsne(meas,'Algorithm','exact','Distance','euclidean');
%subplot(2,2,4)
h = gscatter(Y(:,1),Y(:,2),man_label , 'kmy') ;
title('Oak: train')
igroup = h(3);
igroup.Color = '#D95319';

%  SVM clustering based on training
Mdl = fitcecoc(meas,man_label);

rng default
Mdl = fitcecoc(meas, man_label,'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus'))

%   save('Pop_train_whole_tissue_SVM.m', 'Mdl')

% predict on train
labels_train = predict(Mdl,meas);

%%  import test data OAK 10 bins

man_pheno = OAKfinalSVMdata(:,1);
meas_test = cell2mat(OAKfinalSVMdata(:,2:21));

for ii = 1 : height(OAKfinalSVMdata)
    man_label(ii,1) = man_pheno{ii};
end

% predict test data labels

labels = predict(Mdl,meas_test);

% Plot TSNE of manual on test
rng('default') % for fair comparison
Y = tsne(meas_test,'Algorithm','exact','Distance','euclidean');
%subplot(2,2,4)
h = gscatter(Y(:,1),Y(:,2),man_label_test , 'kmy') ;
title('Oak: test')
igroup = h(3);
igroup.Color = '#D95319';

% Plot TSNE of predicted on test
rng('default') % for fair comparison
Y = tsne(meas_test,'Algorithm','exact','Distance','euclidean');
%subplot(2,2,4)
h = gscatter(Y(:,1),Y(:,2),labels , 'kmy') ;
title('Oak: test')
igroup = h(3);
igroup.Color = '#D95319';


%%  import test data OAK, 20 bins

meas_test = cell2mat(oak_import(:,1:20));

% predict test data labels
labels = predict(Mdl,meas_test);

% Plot TSNE of manual on test
rng('default') % for fair comparison
Y = tsne(meas_test,'Algorithm','exact','Distance','euclidean');
%subplot(2,2,4)
h = gscatter(Y(:,1),Y(:,2),man_label , 'kmy') ;
title('Oak: test')
igroup = h(3);
igroup.Color = '#D95319';

% Plot TSNE of predicted on test
rng('default') % for fair comparison
Y = tsne(meas_test,'Algorithm','exact','Distance','euclidean');
%subplot(2,2,4)
h = gscatter(Y(:,1),Y(:,2),labels , 'kmy') ;
title('Oak: test')
igroup = h(3);
igroup.Color = '#D95319';

%%  import test data IMP130

%meas_test = cell2mat(IMP13010020tilessummary(:,1:20));

meas_test = IMP130svmdata;

% predict test data labels
labels = predict(Mdl,meas_test);

% Plot TSNE of manual on test
rng('default') % for fair comparison
Y = tsne(meas_test,'Algorithm','exact','Distance','euclidean');
%subplot(2,2,4)
h = gscatter(Y(:,1),Y(:,2),man_label , 'kmy') ;
title('Oak: test')
igroup = h(3);
igroup.Color = '#D95319';

% Plot TSNE of predicted on test
rng('default') % for fair comparison
Y = tsne(meas_test,'Algorithm','exact','Distance','euclidean');
%subplot(2,2,4)
h = gscatter(Y(:,1),Y(:,2),labels , 'kmy') ;
title('Oak: test')
igroup = h(3);
igroup.Color = '#D95319';




%%  slide level density bar chart plotting

% Import cutoff confusion matrix, with man calls in columns, cutoff pred in
% rows


%% slide level histogram
desert = double(0);
excluded = double(0);
inflamed = double(0);
total = double(0);

% import CD8 densities into each by hand

edges = logspace(-4,1, 10) ;

dcount = histcounts(desert, edges) ;
ecount = histcounts(excluded, edges) ;
icount = histcounts(inflamed, edges) ;
tcount = histcounts(total, edges) ;
dcount = dcount(2:end);
ecount = ecount(2:end);
icount = icount(2:end);
tcount = tcount(2:end);
edges = edges(3:end);

edges_label = strsplit(num2str(edges, 2));

figure, hold on
b = bar([dcount; ecount; icount;tcount]', 1) 

b(1).FaceColor = '#808080';
b(2).FaceColor = '#0072BD';
b(3).FaceColor = '#D95319';
b(4).FaceColor = 'k';

xlabel('CD8 density')
ylabel('Patient count')
title('IMpassion130 Distribution')
set(b, {'DisplayName'}, {'Desert','Excluded','Inflamed','Total'}')
xticklabels(edges_label)
legend()

%%
binEdges = 12000:100:19000; 
h1 = histogram(desert, edges); 
hold on
h2 = histogram(ecount,edges); 




%% OS prediction high/low

% train SVM on OAK, OS, uncencsored and atezo treated only

os_data = logical(OakOSbins(:,1));
meas = OakOSbins(:,2:end);


rng default
Mdl = fitcsvm(meas,os_data,'OptimizeHyperparameters','all', ...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName', ...
    'expected-improvement-plus'))


% predict on train
labels_train = predict(Mdl,meas);


