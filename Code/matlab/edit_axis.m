function edit_axis(filename)
close all
file2 = ['/Users/reubenv/Thesis/pictures/expfit/' filename '.fig']

open(file2)
%ylabel('Error')
%title('')
%xlabel('\delta')
SaveAsPngEpsAndFig(['/Users/reubenv/Desktop/expfit/' filename])