function [iter,ng,x] = standard_gradient_armijo(d,Q,a,L,x0,lambda,delta,tol,options,print_output)
    iter = 0;
    ng = [];
    x = x0;
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
        g = x - a + lambda.*(Q*tmp);
        normg = norm(g,'fro');
        alpha = options.s;
        xn = x - alpha.*g;
        h = Q'*xn;
        normh = vecnorm(h,2,2);
        i = normh <= delta;
        fn = (1/2)*sum(vecnorm(xn - a,2,2).^2) + lambda*sum((1/(2*delta)).*normh(i).^2) + lambda*sum(normh(~i)-delta/2);
        while fn > f - options.gamma * alpha * normg
            alpha = options.sigma * alpha;
            xn = x - alpha.*g;
            h = Q'*xn;
            normh = vecnorm(h,2,2);
            i = normh <= delta;
            fn = (1/2)*sum(vecnorm(xn - a,2,2).^2) + lambda*sum((1/(2*delta)).*normh(i).^2) + lambda*sum(normh(~i)-delta/2);
        end
        x = xn;
        h = Q'*xn;
        normh = vecnorm(h,2,2);
        i = normh <= delta;
        tmp(i,:) = h(i,:)./delta;
        if size(h(~i,:),1) ~= 0
            tmp(~i,:) = h(~i,:)./normh(~i);
        end
        g = x - a + lambda.*(Q*tmp);
        normg = norm(g,'fro');
        ng = [ng;normg];
        if print_output == true    
            fprintf('Iter = %d\tf = %f\tg = %d\n',iter,fn,normg);
        end
        if normg <= tol
            break;
        end
        if iter > 10000
            break;
        end
    end
end