A_max = 20;



P = randi([2,20]);

crc_generator_pattern = round(rand([1,P]));
crc_generator_pattern(1) = 1;
crc_generator_pattern(end) = 1;

G_P2 = get_crc_generator_matrix_ones(A_max, crc_generator_pattern)

G_P = get_crc_generator_matrix(A_max, crc_generator_pattern)

