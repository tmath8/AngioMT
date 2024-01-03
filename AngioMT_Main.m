tic
clear
clc


%% DEFINE MASS TRANSPORT PARAMETERS---------------------------------------%
Flux = ; % flux of oxygen leaving vessels at edges (mol/s/m2)
s1 = ; % bulk consumption of oxygen in vessels (mol/s/m3)
s2 = ; % bulk consumption of oxygen in tissues
D1 = ; % diffusivity of oxygen in vessels
D2 = ; % diffusivity of oxygen in tissues
phi = ; % partition coefficient of oxygen at vessel-tissue interface
C_in = ; % inlet concentration of oxygen

%% DEFINE PATH TO DEPENDENCIES--------------------------------------------%
addpath('/Users/tanmaymathur/Library/CloudStorage/GoogleDrive-tmath8@tamu.edu/My Drive/TAMU/Graduate Student (Aug 2017-Oct 2022)/Research/Simulations/Vascular Networks/MATLAB scripts/github_repo')    % for poly2mesh
addpath('/Users/tanmaymathur/Library/CloudStorage/GoogleDrive-tmath8@tamu.edu/My Drive/TAMU/Graduate Student (Aug 2017-Oct 2022)/Research/Simulations/Vascular Networks/MATLAB scripts/im2mesh')    % for im2mesh

%% DEFINE PATH TO FOLDER WITH IMAGES--------------------------------------%
path_directory='/Users/tanmaymathur/Desktop/AngioMT/Images';
original_files=dir([path_directory '/*.tif']);

%% CREATE TABLE TO STORE ANCILLARY VARIABLES FOR EACH IMAGE---------------%
File_names = [];

T=array2table(zeros(length(original_files),4));
T.Properties.VariableNames = {'Vascular Oxygen', ...
    'Oxygen Delivery','Flux','Computation time (s)'};

%% START READING IMAGES---------------------------------------------------%
for k=1:length(original_files)
    fprintf('------------------------------Analyzing file #%d',k)
    fprintf('\n')

    t0=toc;
    filename=[path_directory '/' original_files(k).name];

    im = imread(filename); % import image

    im_bin = img_bin(im); % binarizing image

    L = bwlabel(im_bin); % extracting labels in binarized image

    labels_inlet = setdiff(unique([L(:,1);L(:,size(L,2))]),0);

    im2=im_bin-0;
    labels = unique(L);

    for i=1:numel(labels)
        if labels(i) == 0
            im2(L==labels(i)) = 0;
        elseif sum(ismember(labels_inlet,labels(i))) == 1
            im2(L==labels(i)) = 1;
        elseif sum(ismember(labels_inlet,labels(i))) == 0
            im2(L==labels(i)) = 0.5;
        end
    end

    img = {im, im_bin, im2};
    titles = {'Input image','Binarized image','Segmented image'};

    fig = figure();
    tlo = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
    for i = 1:numel(img)
        ax = nexttile(tlo);
        imshow(img{i},'Parent',ax)
        title(titles{i},'FontSize',18)
    end

    fig.Position = [100 100 1300 500];

    %% MESHING THE NETWORKS USING IM2MESH---------------------------------%
    % Define meshing parameters
    tf_avoid_sharp_corner = false;
    tolerance = 1;
    hmax = 500;
    mesh_kind = 'delaunay';
    grad_limit = 1;
    select_phase = [];

    [nodes,element,tnum] = im2mesh( im2, select_phase, tf_avoid_sharp_corner, tolerance, hmax, mesh_kind, grad_limit );
    t1=toc;
    fprintf('Time elapsed in meshing: %d mins %d secs',floor((t1-t0)/60),round((t1-t0)-60*floor((t1-t0)/60)))

    %% EXTRACTING ALL DOMAIN FACES AND NODES------------------------------%
    n_phase= unique(tnum);

    if numel(n_phase)==3
        element_inlet = element(tnum==3,:);
        element_island = element(tnum==2,:);
        element_tissue = element(tnum==1,:);

        unique_nodes_inlet = unique(element_inlet);
        unique_nodes_island = unique(element_island);
        unique_nodes_tissue = unique(element_tissue);

        nodes_inlet = [unique_nodes_inlet nodes(unique_nodes_inlet,:)];
        nodes_island = [unique_nodes_island nodes(unique_nodes_island,:)];
        nodes_tissue = [unique_nodes_tissue nodes(unique_nodes_tissue,:)];
    else
        element_inlet = element(tnum==2,:);
        element_island = [];
        element_tissue = element(tnum==1,:);

        unique_nodes_inlet = unique(element_inlet);
        unique_nodes_tissue = unique(element_tissue);

        nodes_inlet = [unique_nodes_inlet nodes(unique_nodes_inlet,:)];
        nodes_island = [];
        nodes_tissue = [unique_nodes_tissue nodes(unique_nodes_tissue,:)];
    end

    %% PLOTTING MESH AND MESH QUALITY-------------------------------------%
    fig2 = figure();
    tlo = tiledlayout(fig2,1,2,'TileSpacing','compact','Padding','compact');

    fig2.Position = [100 100 1300 500];
    
    
    subplot(1,2,1)
    greyscale = unique(im2);
    for i=1:numel(greyscale)
        if greyscale(i)==0
            patch('Faces', element(tnum==i,1:3), 'Vertices', nodes,'Facecolor','cyan','EdgeAlpha',0.1)
            hold on
        else
            patch('Faces', element(tnum==i,1:3), 'Vertices', nodes,'Facecolor',[greyscale(i) greyscale(i) greyscale(i)],'EdgeAlpha',0.1)
            hold on
        end
    end
    hold off
    axis([0 1001 0 1001])
    legend('Tissue','Disconnected networks','Vessel networks','FontSize',12)
    title('Meshed image with different domains','FontSize',18)
    
    Qual_ele = element_quality(nodes,element);

    subplot(1,2,2)
    patch('Faces', element, 'Vertices', nodes, 'FaceVertexCData',Qual_ele, ...
        'Facecolor','flat','EdgeColor','w','EdgeAlpha',0.05);
    axis([0 1001 0 1001])
    colormap('jet')
    colorbar
    title('Element quality','FontSize',18)

    hold on
    Vessel_edge = intersect([nodes_inlet(:,1);nodes_island(:,1)],nodes_tissue(:,1));
    X_edge=nodes(Vessel_edge,1);
    Y_edge=nodes(Vessel_edge,2);

    plot(X_edge,Y_edge,'g.','MarkerSize',5)
    legend('Element quality','Edge nodes','FontSize',12)
    hold off


    %% DEFINING INLET NODES------------------------------%
    Inlet=[nodes_inlet(nodes_inlet(:,2)==min(nodes(:,1)),1);nodes_inlet(nodes_inlet(:,2)==max(nodes(:,1)),1)];

    %% EXTRACTING ALL EDGE NODES IN THE VESSEL DOMAIN---------------------%
    Edge_only = intersect(nodes_inlet(:,1),nodes_tissue(:,1));
    X=nodes(Edge_only,1);
    Y=nodes(Edge_only,2);

    %% CALCULATING FLUX VECTOR FOR VESSEL DOMAIN--------------------------%
    ee1=zeros(length(element_inlet),1); % ee stores the number of nodes on the edge for each element in domain 1
    for i=1:length(element_inlet)
        ee1(i)=length(intersect(element_inlet(i,:),Edge_only));
    end

    loc1 = [find(ee1(:,1)==3);find(ee1(:,1)==2)];
    edge_element_d1 = element_inlet(loc1,:);

    Qn = zeros(size(element_inlet));

    for i=1:size(edge_element_d1,1)
        [ic,id]=ismember(Edge_only,edge_element_d1(i,:));
        nn=6-sum(id);

        b = find(element_inlet(:,1)==edge_element_d1(i,1) &...
            element_inlet(:,2)==edge_element_d1(i,2) & ...
            element_inlet(:,3)==edge_element_d1(i,3));

        if nn==1
            Qn(b,2)=Flux;
        elseif nn==2
            Qn(b,3)=Flux;
        elseif nn==3
            Qn(b,1)=Flux;
        elseif nn==0
            b12=find(element_inlet(:,2)==edge_element_d1(i,1));
            b13=find(element_inlet(:,3)==edge_element_d1(i,1));
            b21=find(element_inlet(:,1)==edge_element_d1(i,2));
            b23=find(element_inlet(:,3)==edge_element_d1(i,2));
            b31=find(element_inlet(:,1)==edge_element_d1(i,3));
            b32=find(element_inlet(:,2)==edge_element_d1(i,3));

            bb=unique([b12;b13;b21;b23;b31;b32]);
            ele_common=element_inlet(bb,:);
            temp=[];
            for j=1:size(ele_common,1)
                bbb=intersect(ele_common(j,:),edge_element_d1(i,:));
                if length(bbb)==2
                    temp=[temp,bbb];
                end
            end
            node_common=mode(temp);
            if node_common==edge_element_d1(i,1)
                Qn(b,2)=Flux;
            elseif node_common==edge_element_d1(i,2)
                Qn(b,3)=Flux;
            elseif node_common==edge_element_d1(i,3)
                Qn(b,1)=Flux;
            end
        end
    end

    t2=toc;
    fprintf('\n')
    fprintf('Time elapsed in identifying edge nodes and applying flux BCs: %d  mins %d secs',floor((t2-t1)/60),round((t2-t1)-60*floor((t2-t1)/60)))

    %% CALCULATING ELEMENT STIFFNESS AND FORCE MATRICES-------------------%
    [K_g,F_g] = global_matrices(nodes_inlet,element_inlet,D1,Qn,s1);

    t3=toc;
    fprintf('\n')
    fprintf('Time elapsed in calculating stiffness and force matrices for vessel domain: %d  mins %d secs',floor((t3-t2)/60),round((t3-t2)-60*floor((t3-t2)/60)))

    %% APPLYING INLET BOUNDARY CONDITION AND CALCULATING VESSEL OXYGENATION
    G=Inlet;
    g=C_in*ones(length(G),1);

    for i=1:length(G)
        x=find(nodes_inlet(:,1)==G(i));
        F_g(:) = F_g(:) - K_g(:,x)*g(i);
        K_g(:,x) = 0;
        K_g(x,:) = 0;
        K_g(x,x) = 1;
        F_g(x) = g(i);
    end

    C =  K_g\F_g;

    % filtering the negative values
    for i=1:length(C)
        if C(i)<0
            C(i)=0;
        elseif C(i)>C_in
            C(i)=0;
        end
    end

    % making global concentration list across all nodes
    C_vessel = zeros(size(nodes,1),1);
    for i=1:length(nodes_inlet)
        C_vessel(nodes_inlet(i,1))=C(i);
    end

    element_vessel = [element_island;element_inlet];
    nodes_vessel = [nodes_island;nodes_inlet];

    t4=toc;
    fprintf('\n')
    fprintf('Time elapsed in applying inlet BCs for vessel domain and calculating solution: %d  mins %d secs',floor((t4-t3)/60),round((t4-t3)-60*floor((t4-t3)/60)))

    %% CALCULATING ELEMENT STIFFNESS AND FORCE MATRICES FOR DOMAIN 2------%
    Q2 = zeros(size(element_tissue));

    [K_g2,F_g2] = global_matrices(nodes_tissue,element_tissue,D2,Q2,s2);

    t5=toc;
    fprintf('\n')
    fprintf('Time elapsed in calculating stiffness and force matrices for tissue domain: %d  mins %d secs',floor((t5-t4)/60),round((t5-t4)-60*floor((t5-t4)/60)))
    fprintf('\n')
    fprintf('---------')

    G2=Edge_only;
    g2=phi*C_vessel(Edge_only);

    for i=1:length(G2)
        x2 = find(nodes_tissue(:,1)==G2(i));
        F_g2(:) = F_g2(:) - K_g2(:,x2)*g2(i);
        K_g2(:,x2) = 0;
        K_g2(x2,:) = 0;
        K_g2(x2,x2) = 1;
        F_g2(x2) = g2(i);
    end

    tol = 1e-11;
    maxit = 1000;
    C2 = pcg(K_g2,F_g2,tol,maxit);

    % filtering the negative values
    for i=1:length(C2)
        if C2(i)<0
            C2(i)=0;
        elseif C2(i)>8.6
            C2(i)=0;
        end
    end

    % making global concentration list across all nodes of domain 2
    C_tissue = zeros(length(nodes),1);
    for i=1:length(nodes_tissue)
        C_tissue(nodes_tissue(i,1))=C2(i);
    end

    t6=toc;
    fprintf('\n')
    fprintf('Time elapsed in applying BCs for tissue domain and calculating solution: %d  mins %d secs',floor((t6-t5)/60),round((t6-t5)-60*floor((t6-t5)/60)))

    %% CALCULATING AUX VARIABLES------------------------------------------%


    [Flux_centroid1, C_centroid1, Centroid_coord1, Flux_vector1, Area_average_C1, Area_average_Flux1, Area1] = aux_var(nodes_vessel,element_vessel,C_vessel,D1);
    [Flux_centroid2, C_centroid2, Centroid_coord2, Flux_vector2, Area_average_C2, Area_average_Flux2, Area2] = aux_var(nodes_tissue,element_tissue,C_tissue,D2);

    t7=toc;
    fprintf('\n')
    fprintf('Time elapsed in calculating auxillary variables: %d  mins %d secs',floor((t7-t6)/60),round((t7-t6)-60*floor((t7-t6)/60)))

    %% PLOTTING OXYGEN CONCENTRATION AND FLUX OVER THE MESH---------------%
    fig3 = figure();
    tlo = tiledlayout(fig3,1,2,'TileSpacing','compact','Padding','compact');

    fig3.Position = [100 100 1300 500];
    
    subplot(1,2,1)
    patch('Faces', element_vessel, 'Vertices', nodes, 'FaceVertexCData',C_vessel/C_in, ...
        'Facecolor','interp','EdgeColor','none','EdgeAlpha',0.1)
    hold on
    patch('Faces', element_tissue, 'Vertices', nodes, 'FaceVertexCData',C_tissue/C_in, ...
        'Facecolor','interp','Marker','.','EdgeColor','none')
    axis([0 1001 0 1001])
    axis off
    colormap('jet');
    colorbar
    title('Normalized oxygen concentration','FontSize',18)

    subplot(1,2,2)
    F1 = sqrt((Flux_vector1(:,1).^2)+(Flux_vector1(:,2).^2));
    F2 = sqrt((Flux_vector2(:,1).^2)+(Flux_vector2(:,2).^2));

    max_flux = max(max(F1),max(F2));

    patch('Faces', element_vessel, 'Vertices', nodes, 'FaceVertexCData',F1, ...
        'Facecolor','interp','EdgeColor','none','EdgeAlpha',0.05)
    hold on
    patch('Faces', element_tissue, 'Vertices', nodes, 'FaceVertexCData',F2, ...
        'Facecolor','interp','EdgeColor','none','EdgeAlpha',0.05)
    axis([0 size(im2,2)+1 0 size(im2,1)+1])
    axis off
    colormap('jet');
    title('Oxygen flux','FontSize',18)

    hold on
    quiver(Centroid_coord1(:,1),Centroid_coord1(:,2),Flux_centroid1(:,1),Flux_centroid1(:,2),'LineWidth',2,'AutoScaleFactor',1.5)
    hold on
    quiver(Centroid_coord2(:,1),Centroid_coord2(:,2),Flux_centroid2(:,1),Flux_centroid2(:,2),'LineWidth',2,'AutoScaleFactor',1.5)

    %% STORING CALCULATED VARIABLES---------------------------------------%

    T(k,1)={Area_average_C1};
    T(k,2)={Area_average_C2};
    T(k,3)={Area_average_Flux2};

    t8=toc;
    T(k,4)={t8-t0};

    fprintf('\n')
    fprintf('Total time elapsed in analyzing image "%s": %d  mins %d secs',original_files(k).name,floor((t8-t0)/60),round((t8-t0)-60*floor((t8-t0)/60)))
    fprintf('\n')
    fprintf('\n')
end


t9=toc;
fprintf('\n')
fprintf('Total time elapsed in analyzing %d files: %d  mins %d secs',k,floor(t9/60),round(t9-60*floor(t9/60)))

T2={original_files.name}';
T2=array2table(T2);
T2.Properties.VariableNames = {'File name'};
T3=horzcat(T2,T)