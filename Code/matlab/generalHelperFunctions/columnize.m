function C = columnize(R,P)
%converts data from 2 sets of xy valuse into 1 column. Mainly used for
%feeding into inbuilt matlab solvers.
R = R(:);
P = P(:);
C = [R;P];