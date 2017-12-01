function crc = calculate_crc_ones(bits, crc_polynomial_pattern)

crc = ones(1,length(crc_polynomial_pattern)-1);

for bit_index = 1:length(bits)
    crc = xor(xor(crc(1),bits(bit_index))*crc_polynomial_pattern(2:end),[crc(2:end),0]);
    
end
