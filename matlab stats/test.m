clear;

%confusionMatrix = processVideoFolder('H:\changeDetection2014\dataset2014\dataset\baseline\highway', 'H:\changeDetection2012\SOBS_20\results\baseline\highway');
%confusionMatrix = processVideoFolder('H:\changeDetection2014\dataset2014\dataset\baseline\highway', 'D:\sunf\Codes\OpencvTest\Result\Subsensex\baseline\input0'); 
 %confusionMatrix = confusionMatrix + processVideoFolder('H:\changeDetection2014\dataset2014\dataset\baseline\office', 'H:\changeDetection2012\SOBS_20\results\baseline\office');
 %confusionMatrix =  processVideoFolder('H:\changeDetection2014\dataset2014\dataset\baseline\pedestrians', 'D:\sunf\Codes\OpencvTest\Result\Subsensex\baseline\input3');
 %confusionMatrix = confusionMatrix + processVideoFolder('H:\changeDetection2014\dataset2014\dataset\baseline\PETS2006', 'H:\changeDetection2012\SOBS_20\results\baseline\PETS2006');
 %confusionMatrix = processVideoFolder('g:\changeDetection2014\dataset2014\dataset\ptz\continuousPan', 'D:\sunf\Codes\OpencvTest\Result\Subsensex\ptz\input0\gpu'); 
 %confusionMatrix = processVideoFolder('g:\changeDetection2014\dataset2014\dataset\ptz\zoominzoomout', 'D:\sunf\Codes\OpencvTest\Result\Subsensex\ptz\input3\gpu'); 
 %confusionMatrix = processVideoFolder('H:\changeDetection2014\dataset2014\dataset\ptz\zoominzoomout', 'D:\sunf\Codes\Download\fastMCD_v1\resultsZoom');
 %confusionMatrix = processVideoFolder('D:\sunf\Codes\OpencvTest\moseg\cars4', 'D:\sunf\Codes\OpencvTest\Result\Subsensex\moseg\cars4\gpu'); 
 confusionMatrix = processVideoFolder('D:\sunf\Codes\OpencvTest\moseg\cars1', 'D:\sunf\Codes\OpencvTest\Release\warpError\warp1\cars1'); 
 TP = confusionMatrix(1);
 FP = confusionMatrix(2);
 FN = confusionMatrix(3);
 TN = confusionMatrix(4);
 recall = TP / (TP + FN)
 specficity = TN / (TN + FP)
 FPR = FP / (FP + TN)
 FNR = FN / (TP + FN)
 PBC = 100.0 * (FN + FP) / (TP + FP + FN + TN)
 precision = TP / (TP + FP)
 FMeasure = 2.0 * (recall * precision) / (recall + precision)