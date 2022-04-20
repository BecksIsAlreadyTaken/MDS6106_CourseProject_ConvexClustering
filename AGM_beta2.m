function [iter,ng,x] = AGM_beta2(d,Q,a,L,x0,lambda,delta,tol,mu,print_output)
    x = x0;
    xp = x;
    iter = 0;
    ng = [];
    while 1
        iter = iter + 1;
        beta = (sqrt(L)-sqrt(mu))/(sqrt(L)+sqrt(mu));
        y = x + beta*(x-xp);
        xp = x;
        h = Q'*y;
        normh = vecnorm(h,2,2);
        i = normh <= delta;
        f = (1/2)*sum(vecnorm(y - a,2,2).^2) + lambda*sum((1/(2*delta)).*normh(i).^2) + lambda*sum(normh(~i)-delta/2);
        tmp(i,:) = h(i,:)./delta;
        if size(h(~i,:),1) ~= 0
            tmp(~i,:) = h(~i,:)./normh(~i);
        end
        g = y - a + lambda.*(Q*tmp);
        x = y - g./L;
        normg = norm(g,'fro');
        ng = [ng;normg];
        if print_output == true
            fprintf('Iter = %d\t f = %f\t g = %f\n',iter,f,normg);
        end
        if normg < tol
            break;
        end
    end
end