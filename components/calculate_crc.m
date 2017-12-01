function crc = calculate_crc(bits, crc_polynomial_pattern)

crc = zeros(1,length(crc_polynomial_pattern)-1);

for bit_index = 1:length(bits)
    crc = xor(xor(crc(1),bits(bit_index))*crc_polynomial_pattern(2:end),[crc(2:end),0]);
    
end
