%acceptable values
%shoudl use subplot stuff instead
%alpha = [0 0.2 0.4 0.6 0.8 1]
%gamma = [0.1 0.5 1 5]
%T = [0 20 50 200]


alpha = [0.4];
gamma = [5];
T = [0 20 50 200];



for gamma1 = gamma
    for alpha1 = alpha
        for T1 = T
            figure
            imshow([pwd '\strainfigures\' num2str(alpha1) '-' num2str(gamma1) '-' num2str(T1) '.png'])
        end
    end
end