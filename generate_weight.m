function [w] = generate_weight(a,k,theta)
    N = size(a,1);
    labels = 1:1:N;
    w = sparse(N*(N-1)/2,1);
    for i = 1:size(a,1)
        D = [sqrt(sum((a(i,:)-a).^2,2)),labels'];
        D_sorted = sortrows(D,1);
        D_sorted(1,:) = [];
        topk = D_sorted(1:k,:);
        for j = 1:size(topk,1)
            if i > topk(j,2)
                continue;
            end
            ind = (i-1)*(2*size(a,1)-i+2)/2 + topk(j,2) - i - i + 1;
            w(ind) = exp(-theta*D(topk(j,2)));
        end
    end
end