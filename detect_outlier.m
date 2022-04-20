function [outlier_label] = detect_outlier(x)
    outlier_label = [];
    label = unique(x(:,1));
    for i = 1:size(label,1)
        mask = x(:,1) == label(i);
        if size(x(mask,:),1) == 1
            x(mask,:) = [];
            outlier_label = [outlier_label, label(i)];
        end
    end
end