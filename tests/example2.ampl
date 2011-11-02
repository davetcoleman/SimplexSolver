var x1;
var x2;
var x3;
var x4;

maximize objVal: 2*x1 + -1*x2 + 8*x3 + -2*x4;

c1: 10 <= x1 + 2*x2 + x4 <= 30;
c2: 0 <= x1 - x2 - 3*x4;
c3: 6 <= -x2 + x3 <= 8;
c4: -10 <= x3 - x4 <= 10;
c5: 0 <= x1 <= 20;
c6: x2 <= 10;
c7: 0 <= x3;
c8: -10 <= x4 <= 0;

solve;
display x1, x2, x3, x4, objVal;
end;
