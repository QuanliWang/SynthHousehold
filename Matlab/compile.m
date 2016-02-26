if isequal(exist('dist', 'dir'),7)
    %do clean up later
else
    mkdir dist
end

cd src 
mex samplezHHwithHHnewv1_2HHvar.cpp
mex samplehouseholds.cpp
mex mergeindividuals.cpp
mex households2individuals.cpp
mex samplezmemberv1.cpp
mex randomsample.cpp
mex randomsample_new.cpp
mex checkconstraints.cpp
mex groupcount.cpp

!mv *.mexmaci64 ../dist/.
cd ..
