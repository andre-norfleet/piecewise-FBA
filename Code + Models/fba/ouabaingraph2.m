[num, txt, raw] = xlsread("DiBAC_Ouabain_IntensityComparisona_10_11_2019Experiment.xlsx", 'Sheet1', 'A1:B3');

vals = num(:,1);
condition = raw(:,2);

testextract = condition{1,1};
testextract2 = condition{2,1};

x = char(testextract, testextract2, condition{3,1});
x2 = [1,2,3];

bar(x2,vals)
hold on
set(gca,'xticklabel',{'Steady-State','Ouabain','Gramicidin'})
hold off