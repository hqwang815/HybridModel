function [Isr]=calculate_mirror_reflection(R,S)

%% Plane reflection 

% take three points from the plane R 
A=R(1,:);
B=R(2,:);
C=R(3,:);

%Performing the reflection
N=cross(A-B,B-C);
N=N/norm(N);
Isr=S+2*dot(A-S,N)*N;
%Round the reflection because of computer accuracy.
Isr=round(Isr,4);
                
end