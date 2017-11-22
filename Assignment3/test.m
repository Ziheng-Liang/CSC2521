X = [0 0; 16 0; 16 2; 2 2; 2 3; 8 3; 8 5; 0 5];
C = [1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 1];
dt = delaunayTriangulation(X, C);
subplot(2,1,1);
triplot(dt);
axis([-1 17 -1 6]);
xlabel('Constrained Delaunay triangulation', 'fontweight','b');
% Plot the constrained edges in red
hold on;
plot(X(C'),X(C'+size(X,1)),'-r', 'LineWidth', 2);
hold off;

% Now delete the constraints and plot the unconstrained Delaunay
dt.Constraints = [];
subplot(2,1,2);
triplot(dt);
axis([-1 17 -1 6]);
xlabel('Unconstrained Delaunay triangulation', 'fontweight','b');