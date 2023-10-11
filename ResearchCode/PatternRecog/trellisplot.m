function trellisplot(x,n)
figure;
min1 = min(min(min(x)));max1 = max(max(max(x)));
k=0;
for i=1:n,
    for j=1:n,
        k=k+1;    
        if j >= i 
            subplot(n,n,k);plot(x(i,:),x(j,:),'.');
            xlabel(num2str(i));ylabel(num2str(j));
            axis([min1, max1, min1, max1]);
        end
    end
end



