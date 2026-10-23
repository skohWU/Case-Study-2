function [G] = make_star(n)

e_1 = zeros(n)';
e_1(1) = 1;
G = (e_1 * ones(n)') + (ones(n) * e_1') - 2 * (e_1 * e_1');

end