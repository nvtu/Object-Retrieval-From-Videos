% COMPILE compiles the mex files

%% inverted file
fprintf('Compiling Inverted File\n');
mex -g -O -largeArrayDims mxBuildInvFile.cpp ccInvertedFile.cpp
mex -largeArrayDims -g -O mxComputeSelfSim.cpp ccInvertedFile.cpp
mex -largeArrayDims -g -O mxSearchInvFile.cpp ccInvertedFile.cpp
mex -largeArrayDims -g -O mxCleanInvFile.cpp ccInvertedFile.cpp
fprintf('All files compiled successfully\n');

