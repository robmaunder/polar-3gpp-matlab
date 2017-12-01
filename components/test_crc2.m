


while 1
    
    P = randi([2,20]);
    
    crc_generator_pattern = round(rand([1,P]));
    crc_generator_pattern(1) = 1;
    crc_generator_pattern(end) = 1;
    
    
    A = randi([20,100]);
    
    [P,A]
    
         
    crc_ones = calculate_crc_ones(zeros(1,A),crc_generator_pattern);
    
    diff = crc_ones;
    
    
    for rep = 1:100
        
     a = round(rand([1,A]));
    
    crc_zeros = calculate_crc(a,crc_generator_pattern);
    crc_ones = calculate_crc_ones(a,crc_generator_pattern);
    
    
    if ~isequal(diff, xor(crc_zeros,crc_ones))
        [crc1;crc2]
        error('Rob');
    end
    end
end
