function [a, s] = getActivityRatio(responses)
n = length(responses);
top = ((1/n) * sum(responses))^2;
bottom = (1/n) * sum(responses.^2);

a = top / bottom;
s = (1-a) / (1 - 1/n);
end