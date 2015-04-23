clear;
clc;
file = fopen('result.txt','w+');
SumMatrix = [0 0 0 0 0]; % TP FP FN TN SE
 for i=1:8
     input = sprintf('D:\\sunf\\Codes\\OpencvTest\\moseg\\cars%1d',i);
     output = sprintf('D:\\sunf\\Codes\\OpencvTest\\Result\\Subsensex\\moseg\\cars%1d\\gpu',i);
     confusionMatrix = processVideoFolder(input, output);
     TP = confusionMatrix(1);
     FP = confusionMatrix(2);
     FN = confusionMatrix(3);    
     recall = TP / (TP + FN);
     precision = TP / (TP + FP);
     FMeasure = 2.0 * (recall * precision) / (recall + precision);
     fprintf(file,'moseg %d\n',i); 
     fprintf(file,'recall %.3f precision %.3f FMeasure %.3f\n',recall,precision,FMeasure);
     
     SumMatrix = SumMatrix + confusionMatrix;
 end

 
 
 confusionMatrix = processVideoFolder('H:\changeDetection2014\dataset2014\dataset\ptz\continuousPan', 'D:\sunf\Codes\OpencvTest\Result\Subsensex\ptz\input0\gpu');
 TP = confusionMatrix(1);
 FP = confusionMatrix(2);
 FN = confusionMatrix(3);
 recall = TP / (TP + FN);
 precision = TP / (TP + FP);
 FMeasure = 2.0 * (recall * precision) / (recall + precision);
 fprintf(file,'continuousPan\n');
 fprintf(file,'recall %.3f precision %.3f FMeasure %.3f\n',recall,precision,FMeasure);
 SumMatrix = SumMatrix + confusionMatrix;
 confusionMatrix = confusionMatrix + processVideoFolder('H:\changeDetection2014\dataset2014\dataset\ptz\zoominzoomout', 'D:\sunf\Codes\OpencvTest\Result\Subsensex\ptz\input3\gpu');
 TP = confusionMatrix(1);
 FP = confusionMatrix(2);
 FN = confusionMatrix(3);
 recall = TP / (TP + FN);
 precision = TP / (TP + FP);
 FMeasure = 2.0 * (recall * precision) / (recall + precision);
 fprintf(file,'zoominzoomout\n');
 fprintf(file,'recall %.3f precision %.3f FMeasure %.3f\n',recall,precision,FMeasure);
 SumMatrix = SumMatrix + confusionMatrix;
 
 
 TP = SumMatrix(1);
 FP = SumMatrix(2);
 FN = SumMatrix(3);
 TN = SumMatrix(4);
 recall = TP / (TP + FN);
 precision = TP / (TP + FP);
 FMeasure = 2.0 * (recall * precision) / (recall + precision);
 fprintf(file,'average \n');
 fprintf(file,'recall %.3f precision %.3f FMeasure %.3f\n',recall,precision,FMeasure);
 fclose(file);
 