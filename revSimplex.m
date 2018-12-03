function [xres,zopt]=revSimplex(A,b,c)
%[xres,zopt]=revSimplex(A,b,c)
%   This function helps to solve the Linear Programing Problem of the standard form:
%   Maximize z = (c^T)X, 
%   subject to:
%   AX<=b,
%   X>=0,
%   (where b is a non-negative vector)
%   using Revised Simplex Method. 
%   The output, xres is the optimal solution and the zopt is the corresponding
%   value of objective function


% Initializing the variables
xres = null(1);
zopt = null(1);
sizeA = size(A);
[m,n] = size(A);
sizeB = size(b);
sizeC = size(c);

% Adjusting the dimension of b and c
if(sizeB(1)<sizeB(2))
    b = b';
    sizeB = size(b);
end
if(sizeC(1)<sizeC(2))
    c = c';
    sizeC = size(c);
end

% Testing if vector b is non-negative
if(~(min(b)>=0))
    fprintf("In is problem, b must be non-negative!\n")
    return
end

% Testing the dimensional conditions
if (sizeB(1)~=sizeA(1)||sizeB(2)~=1||sizeC(2)~=1||sizeC(1)~=sizeA(2))
    fprintf("Invalid Dimensional Input!\n");
    return
end


% Initializing the variables
c = [c;zeros(m,1)];
A = [A,eye(m)];
BInv = eye(m);
basicIdx = [zeros(1,n),ones(1,m)];
basicIdx = find(basicIdx); %The index of current basic variables
step = 0;

% Implementing the algorithm
while(1)
cB = c(basicIdx,:);
MInv = [[1,cB'*BInv];[zeros(m,1),BInv]];%the inverse of M
result = MInv*[0;b];%this is the combined vectof of z and xB
z = result(1);%current objective value
xB = result(2:m+1);%value of current basic varibales
z_c = (cB)'*BInv*A-c';
xres = zeros(m+n,1);
xres(basicIdx) = xB;

fprintf("step %d,the current feasible solution is:\n",step)
fprintf('%.4f, ',xres);
fprintf(',\nand the value of z is %.4f\n\n',z);

if(min(z_c)>=0)%This problem has optimal solution.
fprintf("\nOptimal Solution Found!\n")
    break;
end

entVar = find(z_c==min(z_c));%the index(es) of entering variable(s)
entVar = entVar(randperm(length(entVar),1));%if it is not unique, randomly pick one
tp = BInv*A(:,entVar);

if(max(tp)<=0)%The optimal solution of this problem is unbounded
    fprintf("The Optimal Solution is UNBOUNDED!\n")
    break;
end

tpRev = tp;
tpRev(tpRev==0) = -1;%the revised tp, in case the demominator could be 0,just make it negative
thetaRatio = xB./tpRev;
thetaRatio(thetaRatio<0) = Inf;% in case choosing the negative one
depVar = find(thetaRatio==min(thetaRatio));
depVar = depVar(randperm(length(depVar),1));%if it is not unique, randomly pick one
eta = tp;
eta = -eta/tp(depVar);
eta(depVar) =  1/tp(depVar);%compute the eta vector
E = eye(m);
E(:,depVar)=eta;%compute the E metrix
BInv = E*BInv;%update the inverse of B
%xB = E*xB;       This line is useless. It is another way to update xB.
basicIdx(depVar)=entVar;% update the index of basic variables
step = step+1;

end

zopt = z;
end
