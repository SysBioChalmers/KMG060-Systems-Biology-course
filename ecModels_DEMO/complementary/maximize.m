% Maximize:
function v = maximize(c,S,b,LB,UB)
[v,~,flag] = linprog(-c,[],[],S,b,LB,UB);   %-c to maximize instead of minimize
v          = v*(flag == 1);                 %Avoids plotting infeasible problems
end
