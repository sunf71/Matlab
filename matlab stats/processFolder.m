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

function processFolder(datasetPath, binaryRootPath)
	% Call your method and the comparator on each video folder.

    if filesys('isRootFolder?', datasetPath) == false,
        disp(['The folder' datasetPath 'is not a valid root folder.']);
        return
    end

    stats = Stats(datasetPath);
    
    categoryList = filesys('getFolders', datasetPath);
    for strCategory = categoryList,
        category = strCategory{1};
        stats.addCategories(category);
        
        categoryPath = fullfile(datasetPath, category);
        filesys('mkdir', fullfile(binaryRootPath, category));

        videoList = filesys('getFolders', categoryPath);
        for strVideo = videoList,
            video = strVideo{1};
            
            videoPath = fullfile(categoryPath, video);
            binaryPath = fullfile(binaryRootPath, category, video);

            if filesys('isValidVideoFolder?', videoPath),
				filesys('mkdir', binaryPath);
                confusionMatrix = processVideoFolder(videoPath, binaryPath);
                stats.update(category, video, confusionMatrix);
            end
        end

        stats.writeCategoryResult(category);
    end
    
    %zip(fullfile(binaryRootPath, 'results.zip'), binaryRootPath);

    stats.writeOverallResults();
end
