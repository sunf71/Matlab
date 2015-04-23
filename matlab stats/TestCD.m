 clear;
 confusionMatrix = [5823	7037260	3481	41145076];
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