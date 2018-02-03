function draw_disk(x,y,r)
n=100;
theta = linspace(0,2*pi,n);
X = x*ones(size(theta))+r*cos(theta);
Y = y*ones(size(theta))+r*sin(theta);
fill(X,Y,'w','LineWidth', 2, 'EdgeColor', 'none');
end
