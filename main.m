clear;clc;

help revSimplex
%This is an example of LPP with unbounded solution.

A = [2,1,-3,1,1;1,4,0,1,2;3,0,0,4,2];
b = [10;20;15]';
c = [2;3;1;1;2]';

[x,z]=revSimplex(A,b,c);
