function y = errols(data,p,classes)
    c = max(classes);
    y = 0;
    for i=1:c,
        temp = data(:,classes==i);
        temp = bsxfun(@minus, temp , mean(temp,2));
        D=svds(temp,p-1);
        y = y + norm(temp,'fro')^2 - norm(D,2)^2;
    end
end