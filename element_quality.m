function [Qual_ele] = element_quality(nodes,element)

Qual_ele = zeros(size(element,1),1);

for k=1:size(element,1)
    x1_e = nodes(element(k,1),1);
    x2_e = nodes(element(k,2),1);
    x3_e = nodes(element(k,3),1);
    y1_e = nodes(element(k,1),2);
    y2_e = nodes(element(k,2),2);
    y3_e = nodes(element(k,3),2);

    l1 = sqrt(((x2_e-x1_e)^2)+((y2_e-y1_e)^2));
    l2 = sqrt(((x3_e-x2_e)^2)+((y3_e-y2_e)^2));
    l3 = sqrt(((x1_e-x3_e)^2)+((y1_e-y3_e)^2));

    M_e = [1 x1_e y1_e; 1 x2_e y2_e; 1 x3_e y3_e];
    A_e = 0.5*det(M_e);

    RMS_edge_e = sqrt(((l1^2)+(l2^2)+(l3^2))/3);

    Qual_ele(k,1) = (4/sqrt(3))*A_e/((RMS_edge_e)^2);
end