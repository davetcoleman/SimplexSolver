var x1;
var x2;
var x3;
var w1;
var w2;
var w3;

maximize objVal: -0.2*x1 - 1.6*w3 - 0.2*w2;

c1: w1 <= 4;
c2: 0 <= x3;
c3: 0 <= x2;
c4: w1 = 1.2*x1 - 0.4*w3 + 0.2*w2;
c5: x3 = -0.2*x1 - 0.6*w3 - 0.2*w2;
c6: x2 = 0.6*x1 - 0.2*w3 - 0.4*w2;
c7: 0 <= x1;
c8: w3 <= -1;
c9: w2 <= -5;

solve;
display x1, x2, x3, w1, w2, w3,  objVal;
end;


