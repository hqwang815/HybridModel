function [mesh,Velements,Belements,Cell_center,Vcell] = Input_mesh(msh_file,lc)
%% Mesh importation from an ASCII file generated by Gmsh %%

% Gmsh mesh file syntax for tetrahedral elements is: 
%
% $MeshFormat$
% version_number file_type data_format
% $EndMeshFormat$
% $Nodes$
% Node_number
% Node_id x y z 
% $EndNodes$
% $Elements$
% Elements_number
% Element_id Element_type Tag_number <tags> Node_id_0 Node_id_1 Node_id_2 Node_id_3
% $EndElements$

fmid=fopen(msh_file,'r');
tline = fgets(fmid); %Opening of the first line of the .msh file 

%% Analysis of the mesh file format %% 
if contains(tline,'$MeshFormat')
    tline = fgets(fmid); %Opening of the second line of the .msh file 
    file_format = 'MSH';
    version_number = tline(1:3);
    file_type = tline(5);
    if file_type == '0'
        file_format = [file_format,' ASCII'];
        if version_number == '2.1'
            num3 = 7;
        else
            num3 = 6;
        end
        num1 = 2;
        num2 = 4;
        
        tline = fgets(fmid); % Skipping the line $EndMeshFormat$
    elseif file_type == '1'
        file_format = [file_format,' binary'];
        disp( ['reading file format ', file_format ,' not implemented'] );
    else
        disp( ['file type ', file_type, ' for file format ', file_format, ' not recognized'] );
    end
elseif contains(tline,'$NOD')
    file_format = 'MSH';
    version_number = '1.0';
    num1 = 3;
    num3 = 6;
elseif contains(tline,'$PostFormat')
    tline = fgets(fmid);
    file_format = 'POS';
    version_number = tline(1:3);
    file_type = tline(5);
    if file_type == '0'
        file_format = [file_format,' ASCII'];
    elseif file_type == '1'
        file_format = [file_format,' binary'];
    end
    disp( ['reading file format ', file_format ,'not implemented'] );
else
    disp('file format not recognized');
end
disp( [file_format, ' file format version ', version_number] );
tline = fgets(fmid); %Skipping the line $Node or $PhysicalNames
%% Sorting of the Physical entities %%
if contains(tline,'$PhysicalNames')
    tline = fgets(fmid);
    mesh.NP = sscanf(tline,'%d'); % Get the number of physical types 
    PhyEnt = cell(mesh.NP,2);
    for i = 1:mesh.NP
        tline = fgets(fmid);
        current_ent = sscanf(tline,'%d %d')'; %Scanning two integers 
        PhyEnt{current_ent(2),1} = current_ent(2); %Extracting the number of the PhyEnt
        PhyEnt{current_ent(2),2} = extractBetween(tline,'"','"'); %Extracting the name of the Boundary
    end
    mesh.ENT = PhyEnt;
    tline = fgets(fmid); %Skipping the line $EndPhysicalNames
else
    mesh.ENT = cell(1,1);
    mesh.ENT{1,1} = 0;
end
%% Sorting of the nodes %%
tline = fgets(fmid);
if contains(tline,'$Nodes') %If loop enabling to skip $Nodes line if complete format
    tline = fgets(fmid);
end
mesh.NN=sscanf(tline,'%d'); % Get the number of nodes


Nid=zeros(mesh.NN,1);
Ncoor=zeros(mesh.NN,3);                   % matrix of cooordinates
for i=1:mesh.NN
    tline = fgets(fmid);
    current_el=sscanf(tline,'%d %f %f %f')'; % Scanning one integer and 3 floats
    Nid(i)=current_el(1);
    Ncoor(i,:)=current_el(2:4);
end
mesh.POS = Ncoor; %Storing the nodes coordinate in an array
%% Sorting of the elements %% 
tline = fgets(fmid); % Skipping the line $EndNodes$
tline = fgets(fmid); % Skipping the line $EndElements$
tline = fgets(fmid); % Accessing the line of the number of elements 
mesh.NE = sscanf(tline,'%d');  % gets number of elements
mesh.NEE = 0; %Initialisation of the number of edge elements
mesh.NEB = 0; %... of boundary elements
mesh.NEV = 0; %...of volume elements
% creates a vector of structures. each structure represents one element and ...
% contains the connectivity (nodes), the faces information, and the 1st
% tag which corresponds to the physical entity englobing the element
Velements = repmat(struct('nodes',[],'neighbours',[],'elset',0),1,mesh.NEV);
Belements = repmat(struct('nodes',[],'elset',0,'Abn',0),1,mesh.NEB);
for i = 1:mesh.NE                        
    tline = fgets(fmid);                     
    current_el = sscanf(tline,'%d')';
    if current_el(num1) == 1
        mesh.NEE = mesh.NEE + 1;
    end
    if current_el(num1) == 2 
        mesh.NEB = mesh.NEB + 1;
        Belements(mesh.NEB).elset = current_el(num2); %Take the num2 first value of h
        Belements(mesh.NEB).nodes = current_el(num3:end); %Take the node-numbers of the boundary area
    end
    if current_el(num1) == 4 
        mesh.NEV = mesh.NEV + 1;
        Velements(mesh.NEV).elset = current_el(num2); %Take the num2 first value of h
        Velements(mesh.NEV).nodes = current_el(num3:end); %Take the node-numbers of the element
    end
end
switch current_el(2)          % all the last elements in the mesh are volume el                     
    case{4}
        mesh.Eltype=4;                   % T4 elements (Tetrahedrons)
        mesh.nfc = 3;                    % number of nodes in one surface element
        mesh.nedg = 2;                   % number of nodes in one line element
        mesh.ng = 1;                     % number of gauss points
end
fclose(fmid);

%% Control volume definition %% 

%To access to the coordinate of one element nodes: 
%disp(mesh.POS(elements(6).nodes(1),:)) % local node "O" of the element 6 coordinates 

Vcell = zeros(mesh.NEV,1); % Vector of the elements volumes
Cell_center = zeros(mesh.NEV,3); % Vector of the elements centers coordinates
for i = 1:mesh.NEV
    % Global node coordinates extraction
    Nc0 = mesh.POS(Velements(i).nodes(1),:);
    Nc1 = mesh.POS(Velements(i).nodes(2),:);
    Nc2 = mesh.POS(Velements(i).nodes(3),:);
    Nc3 = mesh.POS(Velements(i).nodes(4),:); 
    % Cell center coordinate calculation 
    for j = 1:3
    Cell_center(i,j) = (Nc0(j)+Nc1(j)+Nc2(j)+Nc3(j))/4;
    end
    % Cell volume calculation
    Vcell(i) = abs(dot(Nc0-Nc3,cross(Nc1-Nc3,Nc2-Nc3)))/6;    
end


%% Neighbours connectivity %% 
% if two elements has 3 nodes in common, they are neighbours 
wt = waitbar(0,'Computing connectivity nodes');
tic
for i = 1:mesh.NEV
    current_el = Velements(i).nodes;
    current_neig = zeros(1,4);
    for j = 1:mesh.NEV
        if j ~= i
            if norm(Cell_center(i,:)-Cell_center(j,:)) < 0.9*lc
                int = intersect(current_el,Velements(j).nodes);
                if length(int) == mesh.nfc
                    % Looking for the local indices of the common elements
                    idx_1 = find(current_el == int(1));
                    idx_2 = find(current_el == int(2));
                    idx_3 = find(current_el == int(3));
                    list_idx = [idx_1,idx_2,idx_3];
                    list_idx_sorted = sort(list_idx);
                    if list_idx_sorted == [1,2,3]
                        % The neighbor corresponds to the local "first" face:
                        current_neig(1) = j;
                    elseif list_idx_sorted == [2,3,4]
                        % ... "Second" face
                        current_neig(2) = j;
                    elseif list_idx_sorted == [1,2,4]
                        % And so on
                        current_neig(3) = j;
                    elseif list_idx_sorted == [1,3,4]
                        current_neig(4) = j;
                    end
                end
            end
            % Evry elements has a list of his neighbours, if there is a zero in the
            % list, it means that the face is a boundary
            % The list of neighbours is sorted in a way that the first
            % neighbour of the list corresponds to the first local face of the
            % element.
        end
    end
    Velements(i).neighbours = current_neig;
    waitbar(i / mesh.NEV);
end
toc
close(wt)
