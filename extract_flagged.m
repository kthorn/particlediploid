%Instructions: replace data_to_analyze with name of variable to analyze
%save and hit F5 to run (or extract_flagged at command line)

%parameters of flagged dots are accumulated into results
%when done with one type of cell, rename results and fxn_2dots and it will start from
%scratch on next run.

temp_anal = data_to_analyze;

if ~exist('results')
    j=1;
end

for n=1:size(temp_anal,2)
    if any(strcmp(temp_anal(n).flags, 'disappearing'))
        result = [temp_anal(n).model1_I, temp_anal(n).model2_I];
        results(j,:)=result;
        j=j+1;
    end
    
end


