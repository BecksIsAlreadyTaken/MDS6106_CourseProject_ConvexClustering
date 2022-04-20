function [Q] = generate_Q(N)
    Q = sparse(N,N*(N-1)/2);
    csum = [0,cumsum(flip(1:N-1))];
    for i = 1:N-1
        start = csum(i)+1;
        stop = csum(i+1);
        Q(i+1:end,start:stop) = -speye(N-i);
        Q(i,start:stop) = 1;
    end
end