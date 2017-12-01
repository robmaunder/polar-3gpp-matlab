while 1
    
    P = randi([2,20]);
    
    crc_generator_pattern = round(rand([1,P]));
    crc_generator_pattern(1) = 1;
    crc_generator_pattern(end) = 1;
    
    
    A = randi([0,100]);
    
    [P,A]
    
    
    a = round(rand([1,A]));
    
    crc = calculate_crc(a,crc_generator_pattern);
    
    
   
    
    if  ~isequal(calculate_crc([a,fliplr(crc)],crc_generator_pattern),zeros(1,P-1))
        error('Rob');
    end
end
