%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
%FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
%DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
%SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
%CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
%OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% Nil Goyette
% University of Sherbrooke
% Sherbrooke, Quebec, Canada. April 2012

classdef Stats < handle
    %STATS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = private, SetAccess = private)
        path = '';
        categories = 0;
    end
    
    methods
        function this = Stats(path)
            this.path = path;
            this.categories = containers.Map();
        end
        
        function this = addCategories(this, category)
            if this.categories.isKey(category) == false,
                this.categories(category) = containers.Map();
            end
        end
        
        function this = update(this, category, video, confusionMatrix)
            currentCategory = this.categories(category);
            currentCategory(video) = confusionMatrix;
            % it should update this.categories as well
            
            % Save the confusion matrix as well
            [TP FP FN TN SE] = confusionMatrixToVar(confusionMatrix);
            
            f = fopen([this.path '\' category '\' video '\cm.txt'], 'wt');
            fprintf(f, 'cm video %s %s %u %u %u %u %u', category, video, TP, FP, FN, TN, SE);
            fclose(f);
        end
        
        function this = writeCategoryResult(this, category)
            categoryStats = [];
            f = fopen([this.path '\' category '\cm.txt'], 'wt');
            
            currentCategory = this.categories(category);
            for video = keys(currentCategory),
                video = video{1};
                currentVideo = currentCategory(video);
                [TP FP FN TN SE stats] = confusionMatrixToVar(currentVideo);
                categoryStats = [categoryStats; stats];
                fprintf(f, 'cm video %s %s %u %u %u %u %u\n', category, video, TP, FP, FN, TN, SE);
            end
            
            confusionMatrix = sumCells(values(currentCategory));
            [TP FP FN TN SE] = confusionMatrixToVar(confusionMatrix);
            fprintf(f, 'cm category %s %u %u %u %u %u\n\n', category, TP, FP, FN, TN, SE);
            
            fprintf(f, '\nRecall\t\t\tSpecificity\t\tFPR\t\t\t\tFNR\t\t\t\tPBC\t\t\t\tPrecision\t\tFMeasure');
            fprintf(f, '\n%1.10f\t%1.10f\t%1.10f\t%1.10f\t%1.10f\t%1.10f\t%1.10f', mean(categoryStats));
            
            fclose(f);
        end
        
        function this = writeOverallResults(this)
            categoryStats = containers.Map();
            
            f = fopen([this.path '\cm.txt'], 'wt');
            
            for category = keys(this.categories),
                category = category{1};
                categoryStats(category) = [];
                currentCategory = this.categories(category);
                for video = keys(currentCategory),
                    video = video{1};
                    currentVideo = currentCategory(video);
                    [TP FP FN TN SE stats] = confusionMatrixToVar(currentVideo);
                    categoryStats(category) = [categoryStats(category); stats];
                    fprintf(f, 'cm video %s %s %u %u %u %u %u\n', category, video, TP, FP, FN, TN, SE);
                end
                
                confusionMatrix = sumCells(values(currentCategory));
                [TP FP FN TN SE] = confusionMatrixToVar(confusionMatrix);
                fprintf(f, 'cm category %s %u %u %u %u %u\n\n', category, TP, FP, FN, TN, SE);
            end
            
            sumInMatrix = mapToMatrix(this.categories);
            confusionMatrix = sum(sumInMatrix);
            [TP FP FN TN SE] = confusionMatrixToVar(confusionMatrix);
            fprintf(f, '\n\ncm overall %u %u %u %u %u', TP, FP, FN, TN, SE);
            
            overallStats = [];
            fprintf(f, '\n\n\n\n\t\t\tRecall\t\t\tSpecificity\t\tFPR\t\t\t\tFNR\t\t\t\tPBC\t\t\t\tPrecision\t\tFMeasure');
            for category = keys(this.categories),
                category = category{1};
                
                means = mean(categoryStats(category));
                overallStats = [overallStats; means];
                categoryName = category;
                if size(categoryName, 2) > 8,
                    categoryName = strcat(category(1:7), '..');
                end
                fprintf(f, '\n%s :\t%1.10f\t%1.10f\t%1.10f\t%1.10f\t%1.10f\t%1.10f\t%1.10f', categoryName, means);
            end

            fprintf(f, '\n\nOverall:\t%1.10f\t%1.10f\t%1.10f\t%1.10f\t%1.10f\t%1.10f\t%1.10f', mean(overallStats));
            if categoryStats.Count < 6,
                fprintf(f, '\nYour method will not be visible in the Overall section.\nYou need all 6 categories to appear in the overall section.');
            end
            
            fclose(f);
        end
    end
    
end

function total = sumCells(cells)
    total = [0 0 0 0 0];
    for cell = cells,
        total = total + cell{1};
    end
end


function mat = mapToMatrix(map)
    mat = [];
    maps = values(map);
    for idx = 1:map.length(),
        newPart = mapToRows(maps{idx});
        mat = [mat; newPart];
    end
end

function rows = mapToRows(map)
    rows = [];
    for value = values(map),
        rows = [rows; value{1}];
    end
end

function [TP FP FN TN SE stats] = confusionMatrixToVar(confusionMatrix)
    TP = confusionMatrix(1);
    FP = confusionMatrix(2);
    FN = confusionMatrix(3);
    TN = confusionMatrix(4);
    SE = confusionMatrix(5);
    
    recall = TP / (TP + FN);
    specficity = TN / (TN + FP);
    FPR = FP / (FP + TN);
    FNR = FN / (TP + FN);
    PBC = 100.0 * (FN + FP) / (TP + FP + FN + TN);
    precision = TP / (TP + FP);
    FMeasure = 2.0 * (recall * precision) / (recall + precision);
    
    stats = [recall specficity FPR FNR PBC precision FMeasure];
end
