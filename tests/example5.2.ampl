var x1;
var x2;
var x3;
var e1;
var e2;
var w1;
var w2;
var w3;

maximize objVal: -1*e1 - e2;


c1: w1 <= 4;
c2: w2 <= -5;
c3: w3 <= -1;
c4: w1 = 2*x1 - x2 + x3 - e1;
c5: w2 = 2*x1 - 3*x2 + x3 - e1;
c6: w3 = -x1 + x2 - 2*x3 - e2;
c7: x1 >= 0;
c8: x2 >= 0;
c9: x3 >= 0;
c10: 0 <= e1 <= 5;
c11: 0 <= e2 <= 1;

solve;
display x1, x2, x3, e1, e2, w1, w2, w3,  objVal;
end;


