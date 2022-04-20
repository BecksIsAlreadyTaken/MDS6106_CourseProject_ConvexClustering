function [P,K] = generate_K(h,Q,delta,lambda)
    N = size(Q,1);
    H = sparse(1./max(vecnorm(h,2,2),delta));
    P = sparse(speye(N) + sparse(lambda.*Q*(H.*Q')));
    d = size(h,2);
    K = kron(P,speye(d,d));
%     K = sparse(N*d,N*d);
%     for i = 1:N-1
%         p = [sparse(1,i) diag(P,i)'];
%         e = ones(d,1);
%         h = reshape(p.*e,[],1).*ones(1,N*d);
%         temp = spdiags(h,i*d,N*d,N*d);
%         K = K + temp;
%     end
%     K = K + K' + spdiags(reshape(diag(P,0)'.*ones(d,1),[],1).*ones(1,N*d),0,N*d,N*d);
end