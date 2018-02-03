function y = solve_system_exact(f, t, y_start)

%determine size of the system
m = length(y_start);

%pre-allocate memory
y = zeros(m, length(t));

for n = 1:length(t)
    y(:,n) = f(t(n));
end
end

