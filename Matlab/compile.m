if isequal(exist('dist', 'dir'),7)
    %do clean up later
else
    mkdir dist
end

cd src 
mex sampleindividuals.cpp
mex samplezHHwithHHnewv1_2HHvar.cpp
mex samplezmemberv1.cpp
mex randomsample.cpp
mex randomsample_new.cpp
mex checkingconstraints.cpp
mex groupcount.cpp

!mv *.mexmaci64 ../dist/.
cd ..
