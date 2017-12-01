A_max = 100;


while 1
    
    P = randi([2,20]);
    
    crc_generator_pattern = round(rand([1,P]));
    crc_generator_pattern(1) = 1;
    crc_generator_pattern(end) = 1;

    G_P2 = get_crc_generator_matrix_ones(A_max, crc_generator_pattern);
    
    
    A = randi([20,A_max]);
    
    [P,A]
    
    G_P = get_crc_generator_matrix(A, crc_generator_pattern);
    
    a = round(rand([1,A]));
        
%    crc1 = xor(mod(a*G_P,2),calculate_crc_ones(zeros(1,A),crc_generator_pattern));

    crc1 = xor(mod(a*G_P,2),G_P2(A,:));
    crc2 = calculate_crc_ones(a,crc_generator_pattern);
    
    if ~isequal(crc1,crc2)
        [crc1;crc2]
        error('Rob');
    end
end
