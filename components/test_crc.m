while 1
    
%     P = randi([2,20]);
%     
%     crc_generator_pattern = round(rand([1,P+1]));
%     crc_generator_pattern(1) = 1;
%     crc_generator_pattern(end) = 1;

     P=24;
     crc_generator_pattern = [1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1];
    
    
    A = randi([0,100]);
    [P,A]
    a = round(rand([1,A]));
    
    
    
    
    toggle = calculate_crc(ones(1,P),crc_generator_pattern)
    
    if A >= P
        a2 = [xor(a(1:P),toggle),a(P+1:end)];
        crc3 = calculate_crc(a2,crc_generator_pattern);
    else
        a2 = xor(a,toggle(1:A));
        crc3 = calculate_crc(a2,crc_generator_pattern);
        crc3 = [xor(crc3(1:P-A),toggle(A+1:end)),crc3(P-A+1:end)];
   end
        
       
        
    
    G_P = get_crc_generator_matrix(A+P, crc_generator_pattern);
    
        
    crc1 = mod([ones(1,P),a]*G_P,2);
    
    crc2 = calculate_crc([ones(1,P),a],crc_generator_pattern);
    
    
    
    
    
    if ~isequal(crc1,crc2)
        [crc1;crc2]
        error('Rob');
    end
    
        if ~isequal(crc1,crc3)
            [size(crc1),size(crc3)]
        [crc1;crc3]
        error('Rob2');
    end

end
