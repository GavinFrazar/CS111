function draw_disk(x,y,r)
n=100;
theta = linspace(0,2*pi,n);
X = x*ones(size(theta))+r*cos(theta);
Y = y*ones(size(theta))+r*sin(theta);
fill(X,Y,'r','LineWidth', 1); % fills (in red 'r') the 2-D polygon defined by vectors X and Y
axis square;
% title('A Bouncing Ball','FontSize',14);
set(gca,'xtick',[]); set(gca,'xticklabel',[]);
set(gca,'ytick',[]); set(gca,'yticklabel',[]);
end
