while 1
    
    P = randi([2,20]);
    
    crc_generator_pattern = round(rand([1,P]));
    crc_generator_pattern(1) = 1;
    crc_generator_pattern(end) = 1;
    
    
    A = randi([0,2000]);
    
    [P,A]
    
    G_P = get_crc_generator_matrix(A, crc_generator_pattern);
    
    a = round(rand([1,A]));
    
    crc1 = mod(a*G_P,2);
    crc2 = calculate_crc(a,crc_generator_pattern);
    
    if ~isequal(crc1,crc2)
        [crc1;crc2]
        error('Rob');
    end
end
