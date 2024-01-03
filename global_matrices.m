function [K_g,F_g] = global_matrices(nodes,element,D,Q,s)

K_g= zeros(length(nodes));
F_g = zeros(length(nodes),1);

for k=1:size(element,1)
    x1_e = unique(nodes(find(nodes(:,1)==element(k,1)),2));
    x2_e = unique(nodes(find(nodes(:,1)==element(k,2)),2));
    x3_e = unique(nodes(find(nodes(:,1)==element(k,3)),2));
    y1_e = unique(nodes(find(nodes(:,1)==element(k,1)),3));
    y2_e = unique(nodes(find(nodes(:,1)==element(k,2)),3));
    y3_e = unique(nodes(find(nodes(:,1)==element(k,3)),3));
    
    l1 = sqrt(((x2_e-x1_e)^2)+((y2_e-y1_e)^2));
    l2 = sqrt(((x3_e-x2_e)^2)+((y3_e-y2_e)^2));
    l3 = sqrt(((x1_e-x3_e)^2)+((y1_e-y3_e)^2));
    
    M_e = [1 x1_e y1_e; 1 x2_e y2_e; 1 x3_e y3_e];
    A_e = 0.5*det(M_e);
    
    B_e = (0.5/A_e)*[y2_e-y3_e y3_e-y1_e y1_e-y2_e; x3_e-x2_e x1_e-x3_e x2_e-x1_e];
    
    K_e = D*A_e*B_e'*B_e;
    F_e_b = s*A_e*(1/3)*[1:1:1];
    F_e_t = Q(k,1)*l1*(1/2)*[1;1;0] + Q(k,2)*l2*(1/2)*[0;1;1] + Q(k,3)*l3*(1/2)*[1;0;1];
    
    F_e = F_e_b+F_e_t;
    
    n_npe = size(element,2);
    
    for i=1:n_npe;
        for j=1:n_npe;
            K_g(find(nodes(:,1)==element(k,i)),find(nodes(:,1)==element(k,j))) = K_g(find(nodes(:,1)==element(k,i)),find(nodes(:,1)==element(k,j))) + K_e(i,j);
        end
        
        F_g(find(nodes(:,1)==element(k,i))) = F_g(find(nodes(:,1)==element(k,i))) + F_e(i);
    end
   
    
end