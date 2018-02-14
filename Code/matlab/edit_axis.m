function edit_axis(filename)
file2 = ['E:\Thesis\pictures\expfit\' filename '.fig']
open(file2)

title('')
%ylabel('Error')


SaveAsPngEpsAndFig(['C:\Users\Qazadex\Desktop\pics\' filename])
close all
file2 = ['/Users/reubenv/Thesis/pictures/expfit/' filename '.fig'];

open(file2)
%ylabel('Error')
%title('')
%xlabel('\delta')
SaveAsPngEpsAndFig(['/Users/reubenv/Desktop/expfit/' filename])
