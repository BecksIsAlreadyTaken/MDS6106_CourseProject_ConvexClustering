function [iter,ng,x] = newton_cg(d,Q,a,x0,lambda,delta,tol,options,print_output)
    x = x0;    
    iter = 0;
    ng = [];
    while 1
        iter = iter + 1;
        h = Q'*x;
        normh = vecnorm(h,2,2);
        i = normh <= delta;
        f = (1/2)*sum(vecnorm(x - a,2,2).^2) + lambda*sum((1/(2*delta)).*normh(i).^2) + lambda*sum(normh(~i)-delta/2);
        tmp(i,:) = h(i,:)./delta;
        if size(h(~i,:),1) ~= 0
            tmp(~i,:) = h(~i,:)./normh(~i);
        end
        [P,K] = generate_K(h,Q,delta,lambda);
        g = P*x - a;
        normg = norm(g,'fro');
        cg_tol = min(1,normg^(0.1))*normg;
        v = sparse(size(x,1)*d,1);
        r = reshape(g.',size(x,1)*d,[]);
        p = -r;
        t = -g;
        dlt = 0;
        for j = 1:options.maxit
            tmp1 = K*p;
            tmp2 = p'*tmp1;
            if tmp2 <= 0
                if j == 1
                    t = -g;
                else
                    t = reshape(v',d,size(x,1))';
                end
                break;
            end
            normr = norm(r);
            dlt = normr^2 / tmp2;
            v = v + dlt.*p;
            r = r + dlt.*tmp1;
            normr_n = norm(r);
            if normr_n <= cg_tol
                t = reshape(v',d,size(x,1))';
                break;
            end
            beta = normr_n^2 / normr^2;
            p = -r + beta.*p;
            if j == options.maxit
                t = reshape(v',d,size(x,1))';
            end
        end
        alpha = options.s;
        xn = x + alpha.*t;
        h = Q'*xn;
        normh = vecnorm(h,2,2);
        i = normh <= delta;
        fn = (1/2)*sum(vecnorm(xn - a,2,2).^2) + lambda*sum((1/(2*delta)).*normh(i).^2) + lambda*sum(normh(~i)-delta/2);
        while fn > f - options.gamma * alpha * reshape(g.',size(x,1)*d,[])' * reshape(t.',size(x,1)*d,[])
            alpha = options.sigma * alpha;
            xn = x + alpha.*t;
            h = Q'*xn;
            normh = vecnorm(h,2,2);
            i = normh <= delta;
            fn = (1/2)*sum(vecnorm(xn - a,2,2).^2) + lambda*sum((1/(2*delta)).*normh(i).^2) + lambda*sum(normh(~i)-delta/2);
        end
        x = xn;
        g = P*xn - a;
        normg = norm(g,'fro');
        ng = [ng;normg];
        if print_output == true
            fprintf('Iter = %d\tf = %f\tg = %d\tj = %d\n',iter,fn,normg,j);
        end
        if normg <= tol
            break;
        end
    end
end