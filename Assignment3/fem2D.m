function fem2D()

%create square mesh
% leave only one line from 5-7 uncomment to choose a model
[tri, v] = createSquareMesh(3, 3, 1,1);
% [tri, v] = createCustomMesh();
% [tri, v] = createCircleMesh(5, 16);
% trimesh(tri, v(:,1), v(:,2));
trimesh(tri, v(:,1), v(:,2));
numVerts = size(v,1);
pos = zeros(numVerts,2);

%floor is at zero
odeFun = @(time, state)(femOde(time,state,tri,v,pos));
outputFun = @(time, state, flag)(femOutputFcn(time,state,flag,pos, tri, v));
options = odeset('OutputFcn', outputFun);

y0 = [reshape(v', 2*numVerts,1); zeros(2*numVerts,1)];
[time, state] = ode45(odeFun, [0, 1000], y0, options);

pos(1:numVerts,:) = [state(end, 1:2:2*numVerts)' state(end, 2:2:2*numVerts)'];
trimesh(tri, pos(:,1), pos(:,2));
end

function y = femOde(time, state, tri, v, pos)
    numVerts = size(v,1);
    pos(1:numVerts,:) = [state(1:2:2*numVerts) state(2:2:2*numVerts)];
    
    %hack to fix object to the floor
    %just set y positions of vertices to zero
    floorPoints = v(:,2) < 1e-6; %some arbitrary toleranc
    underGround = pos(:,2) < 1e-6;
    %set positions to zero
    pos(underGround,2) = 0;
    %set velocities to zero
    vel = [state((2*numVerts+1):2:end) state((2*numVerts+2):2:end)];
    neg = vel(:,2) < 1e-6;
    op = underGround.*neg;
    op = op == 1;
    %only set velocity to zero if the point is below or on the ground
    %and the y velocity is negative
    vel(op,2) = 0;
    y = [reshape(vel', 2*numVerts,1) ; femAccelerations(tri, v, pos)];
    
end

%compute FEM accelerations
function a = femAccelerations(tri, v, pos)
%tri is ref, pos is world
numTris = size(tri, 1);
numVerts = size(v,1);
f = zeros(numVerts,2);
m = zeros(numVerts, 2);
density = 1.0;
g = [0 -9.8]';

for i=1:numTris
    
    %Edge vectors
    E1 = (v(tri(i,2), :)-v(tri(i,1), :))';
    E2 = (v(tri(i,3), :)-v(tri(i,1), :))';
    
    e1 = (pos(tri(i,2), :)-pos(tri(i,1), :))';
    e2 = (pos(tri(i,3), :)-pos(tri(i,1), :))';
    e3 = (pos(tri(i,3), :)-pos(tri(i,2), :))';
    area = 0.5*det([e1 e2]);
    me = density*area;
    
    %%%%%%%%%%ASSIGNMENT%%%%%%%%%%%%%%
    %compute F
    F = ([e1, e2]) / [E1, E2];
%      F = [e1, e2] * inv([E1, E2]);
    
    %%%%%%%%%%ASSIGNMENT%%%%%%%%%%%%%%
    %edit the cauchyStress method to add material models
    stress = cauchyStress(F);
    
    %compute forces for each edge
    
    %edge 1
    n1 = [e1(2) ; -e1(1)];
%     n1 = (n1 / norm(n1)) * norm(e1);
    %%%%%%%%%%ASSIGNMENT%%%%%%%%%%%%%%
    %Set fe1 to the force on the first edge of the triangle (from node 1 to
    %node 2)
    fe1 = stress * n1;
    
    f(tri(i,2), :) = f(tri(i,2), :) - 0.5*fe1';
    f(tri(i,1), :) = f(tri(i,1), :) - 0.5*fe1';
    
    %edge2
    n2 = [e3(2) ; -e3(1)];
%     n2 = (n2 / norm(n2)) * norm(e3);
    
    %%%%%%%%%%ASSIGNMENT%%%%%%%%%%%%%%
    %Set fe2 to the force on the first edge of the triangle (from node 2 to
    %node 3)
    fe2 = stress * n2;
    f(tri(i,3), :) = f(tri(i,3), :) - 0.5*fe2';
    f(tri(i,2), :) = f(tri(i,2), :) - 0.5*fe2';
    
    
    n3 = [-e2(2) ; e2(1)];
%     n3 = (n3 / norm(n3)) * norm(e2);
    
    %%%%%%%%%%ASSIGNMENT%%%%%%%%%%%%%%
    %Set fe3 to the force on the first edge of the triangle (from node 3 to
    %node 1)
    fe3 = stress * n3;
    f(tri(i,3), :) = f(tri(i,3), :) - 0.5*fe3';
    f(tri(i,1), :) = f(tri(i,1), :) - 0.5*fe3';
    
    
    %distribute mass to all vertices
    m(tri(i,1), : ) = m(tri(i,1), : ) + [me me]/3;
    m(tri(i,2), : ) = m(tri(i,2), : ) + [me me]/3;
    m(tri(i,3), : ) = m(tri(i,3), : ) + [me me]/3;
   
    
end
%reshape force vector and compute accelerations, add gravity here
% underGround = pos(:,2) < 1e-6;
% u2 = reshape([underGround ,zeros(numVerts, 1)]', 2*numVerts, 1);
% u2 = u2 == 1;
% f(underGround, 2) = max(0, -f(underGround, 2));
f = reshape(f', 2*numVerts,1);
m = reshape(m', 2*numVerts,1);
% fullG = repmat(g, numVerts,1);
% fullG(u2) = -1;
a = f./m + repmat(g, numVerts,1);
end

function status = femOutputFcn(time, state, flag, pos, tri, v)

if strcmp(flag, 'done') == 0
numVerts = size(v,1);
hold on
clf
pos(1:numVerts,:) = [state(1:2:2*numVerts, end) state(2:2:2*numVerts, end)];
trimesh(tri, pos(:,1), pos(:,2));
hold off
drawnow
end
status = 0;
end



function stress = cauchyStress(F)
    
    %%%%%%%%%%ASSIGNMENT%%%%%%%%%%%%%%
    %Set P = to the Piola-Kirchoff 1 stress for a neohookean model.
    %Use parameters mu = 50, lambda = 50 to test
    mu = 200;
    lambda = 200;
    J = det(F);
    invFT = inv(F');
    P = mu * (F - invFT) + lambda * log(J) .* invFT;
    stress = (P*F')./J;
    
end

function [tri, v] = createSquareMesh(width, height, dx, dy)

[X,Y] = meshgrid(1:dx:width, 1:dy:height);

v = [X(:) Y(:)];
tri = delaunay(v(:,1), v(:,2));
numVerts = size(v,1);
%center mesh
v(:,1) = v(:,1) - sum(v(:,1))./numVerts;
v(:,2) = v(:,2) - min(v(:,2));

end

function [tri, v] = createCircleMesh(r, l)
% v = [0, r];
v = [0, r];
for i = 1:l
    dg = 2 * pi * i / l;
    v = [cos(dg)*r/2, sin(dg)*r/2+r;v];
end

for i = 1:l
    dg = 2 * pi * i / l;
    v = [cos(dg)*r, sin(dg)*r+r;v];
end

tri = delaunay(v(:,1), v(:,2));
numVerts = size(v,1);
%center mesh
v(:,1) = v(:,1) - sum(v(:,1))./numVerts;
v(:,2) = v(:,2) - min(v(:,2));

end

function [tri, v] = createCustomMesh()
    r2 = sqrt(2) / 2;
    v = [1 0;2 0;2 4;2 0;3 0;3 5;3 9;6 9;6 10;3 10;1 10;-2 10;-2 9;1 9;1 5;
        2 11;2 10;3 11;2+r2 11+r2;2 12;2-r2 11+r2;1 11];
    tri = [1 2 3;
        3 4 5;
        5 6 3;
        6 15 3;
        15 1 3;
        6 7 15;
        7 14 15;
        7 8 9;
        9 10 7;
        17 11 14;
        14 10 17;
        17 16 11;
        14 7 10;
        11 12 14;
        12 13 14;
        16 17 10;
        18 16 10;
        19 16 18;
        20 16 19;
        21 16 20;
        22 16 21;
        11 16 22];
end
