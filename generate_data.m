function [X] = generate_data(d,p,n,c,sigma)
    rng('default');
    X = [];
    for i = 1:p
        ai = c(i,:) + sigma(i).*randn(n,2);
        X = [X;ai];
    end
end