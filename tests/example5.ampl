var x1;
var x2;
var x3;

maximize objVal: x1 - x2 + x3;

c1: 2*x1 - x2 + x3 <= 4;
c2: 2*x1 - 3*x2 + x3 <= -5;
c3: -1*x1 + x2 - 2*x3 <= -1;
c4: x1 >= 0;
c5: x2 >= 0;
c6: x3 >= 0;

solve;
display x1, x2, x3, objVal;
end;