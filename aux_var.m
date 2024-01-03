function [Flux_centroid, C_centroid, Centroid_coord, Flux_vector, Area_average_C, Area_average_Flux, Area] = aux_var(nodes,element,C,D)

Flux_vector=zeros(size(C,1),2);
C_centroid = zeros(length(element),1);
Conc_times_area = 0;
Flux_times_area = 0;
Area = 0;
Area_average_Flux = 0;
Flux_centroid = zeros(length(element),2);
Centroid_coord = zeros(length(element),2);

for k=1:size(element,1)
    n1_e = element(k,1);
    n2_e = element(k,2);
    n3_e = element(k,3);

    x1_e = unique(nodes(find(nodes(:,1)==element(k,1)),2));
    x2_e = unique(nodes(find(nodes(:,1)==element(k,2)),2));
    x3_e = unique(nodes(find(nodes(:,1)==element(k,3)),2));
    y1_e = unique(nodes(find(nodes(:,1)==element(k,1)),3));
    y2_e = unique(nodes(find(nodes(:,1)==element(k,2)),3));
    y3_e = unique(nodes(find(nodes(:,1)==element(k,3)),3));

    xc_e = (x1_e+x2_e+x3_e)/3;
    yc_e = (y1_e+y2_e+y3_e)/3;

    Centroid_coord(k,:)=[xc_e yc_e];

    M_e = [1 x1_e y1_e; 1 x2_e y2_e; 1 x3_e y3_e];
    A_e = 0.5*det(M_e);

    B_e = (0.5/A_e)*[y2_e-y3_e y3_e-y1_e y1_e-y2_e; x3_e-x2_e x1_e-x3_e x2_e-x1_e];

    N_ec = (0.5/A_e)*[x2_e*y3_e - x3_e*y2_e + (y2_e-y3_e)*xc_e + (x3_e-x2_e)*yc_e;
        x3_e*y1_e - x1_e*y3_e + (y3_e-y1_e)*xc_e + (x1_e-x3_e)*yc_e;
        x1_e*y2_e - x2_e*y1_e + (y1_e-y2_e)*xc_e + (x2_e-x1_e)*yc_e];

    C_ce = N_ec'*[C(n1_e);C(n2_e);C(n3_e)];
    Flux_ce = -D*[C(n1_e) C(n2_e) C(n3_e)]*B_e';
    C_centroid(k) = C_ce;
    Flux_centroid(k,:) = Flux_ce;

    Area = Area + A_e;
    Conc_times_area = Conc_times_area+(A_e*C_ce);
    Flux_times_area = Flux_times_area+((A_e)*sqrt((Flux_ce(1)^2)+(Flux_ce(2)^2)));
    Flux_vector(n1_e,:)=Flux_vector(n1_e,:)+Flux_ce;
    Flux_vector(n2_e,:)=Flux_vector(n2_e,:)+Flux_ce;
    Flux_vector(n3_e,:)=Flux_vector(n3_e,:)+Flux_ce;

end
Area_average_C  = Conc_times_area/Area;
Area_average_Flux = Flux_times_area/Area;
end

