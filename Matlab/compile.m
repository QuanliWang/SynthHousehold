if isequal(exist('dist', 'dir'),7)
    %do clean up later
else
    mkdir dist
end

cd src 
mex samplezHHwithHHnewv1_2HHvar.cpp
mex samplezmemberv1.cpp
mex randomsample.cpp
%mex checkingconstraints_size2_sorted.cpp
%mex checkingconstraints_size3_sorted.cpp
%mex checkingconstraints_size4_sorted_part1.cpp
%mex checkingconstraints_size4_sorted_part2.cpp
%mex checkingconstraints_size4_sorted_part3.cpp
mex checkingconstraints.cpp

!mv *.mexmaci64 ../dist/.
cd ..
