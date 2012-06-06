function T = mytoeplitz(r,c)
    if(~(r(1)==c(1)))
        T = NaN(length(r),length(c));
    end
    
    rt = truncate(r,eps); 
    if(length(rt)<length(r)/2)
        T1 = toeplitz(rt,rt);
        mult = ceil(length(r)/length(rt));
        T=zeros(mult*length(rt));
        for j=1:mult
           T((j-1)*length(rt)+1:j*length(rt),(j-1)*length(rt)+1:j*length(rt)) = T1;  
        end
        %T = kron(eye(mult),T1);
        H = hankel(rt(2:end));
        [tm,~]=size(T1); [hm,hn]=size(H);
        for j=1:mult-1
            T(j*tm+1:hm+j*tm,tm*j+1-hm:tm*j)=T(j*tm+1:hm+j*tm,tm*j+1-hm:tm*j) + fliplr(H);
            T(2+(j-1)*tm:hm+1+(j-1)*tm,tm*j+1:tm*j+hm)=T(2+(j-1)*tm:hm+1+(j-1)*tm,tm*j+1:tm*j+hm)+ flipud(H);
        end
        %T = tril(T) + tril(T).' - diag(diag(T))
        T = T(1:length(r),1:length(r));
    else
       T = toeplitz(r,r); 
    end
end