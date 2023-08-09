function [map2] = dFC_utils_caps_transform_mask(map1,ind1,ind2)
n1 =length(map1);
n2 = size(ind2,2);
map2 = zeros(n2,1);
for v1 = 1:n1
    for v2 = 1:n2
        if ind1(:,v1)==ind2(:,v2)
            map2(v2) = map1(v1);
        end
    end
end

end %functon