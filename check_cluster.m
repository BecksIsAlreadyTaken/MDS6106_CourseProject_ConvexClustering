function [result,result_xstar,label,outlier_label] = check_cluster(x,a,epsilon,handle_outlier)
    N = size(x,1);
    labels = zeros(N,1);
    for i = 1:N
        h = vecnorm(x-x(i,:),2,2);
        labels(h<epsilon) = i;
    end
    label = unique(labels);
    result = sortrows([labels,a],1);
    result_xstar = [result(:,1),x];
    outlier_label = [];
    if handle_outlier == true
        [outlier_label] = detect_outlier(result);
        mask = logical(sum(result(:,1)==outlier_label,2));
        result(mask,:) = [];
        result_xstar(mask,:) = [];
        label(logical(sum(label==outlier_label,2))) = [];
    end
end