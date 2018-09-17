qqqq  = [];
for nn = 1:19
    kkk = find(q{k}>0);
    temp = q{k};
    qqq{nn} = temp(kkk);
    qqqq = [qqqq;qqq{nn}];
end
hist(qqqq)