# polar-3gpp-matlab
Matlab simulations of the encoder and SCL decoder for the New Radio polar code from 3GPP Release 15

At the time of writing, the most recent version of the relevant 3GPP standard is [TS38.212 V1.2.1](https://list.etsi.org/scripts/wa.exe?A3=ind1712B&L=3GPP_TSG_RAN_WG1&E=base64&P=176434682&B=--_004_543CF4C91C60E844AE997DC07D4CC64501A9F092DGGEML504MBSchi_&T=application%2Fx-zip-compressed;%20name=%22R1-1721342.zip%22&N=R1-1721342.zip&attachment=q&XSS=3).

Section of TS38.212 | Implemented in | Comment
--- | --- | ---
5.1 |  components/get_crc_generator_matrix.m | The CRC bits can be generated using b = [a, mod(a*G_P, 2)]
5.2.1 | PUCCH_encoder.m |
5.3.1 | components/get_3GPP_N.m |
5.3.1.1 | components/get_3GPP_crc_interleaver_pattern.m | Interleaving can be implemented using c_prime = c(Pi)
5.3.1.2 | components/PCCA_polar_encoder.m | Other components/\*_polar_encoder.m files are also useful for special cases without PC bits, without CRC bits or with distributed CRC bits
5.4.1.1 P | components/get_3GPP_rate_matching_pattern.m |
5.4.1.1 Q | components/get_3GPP_info_bit_pattern.m |
5.4.1.2 | components/get_3GPP_rate_matching_pattern.m |
5.4.1.3 | components/get_3GPP_channel_interleaver_pattern.m |
5.5 | PUCCH_encoder.m |
6.3.1.2.1 | PUCCH_encoder.m |
6.3.1.3.1 | PUCCH_encoder.m |
6.3.1.4.1 | PUCCH_encoder.m | Rate matching is implemented, but not the determination of E_UCI.
6.3.1.5 | PUCCH_encoder.m |
6.3.2.2.1 | PUCCH_encoder.m |
6.3.2.3.1 | PUCCH_encoder.m |
6.3.2.4.1 | PUCCH_encoder.m | Rate matching is implemented, but not the determination of E_UCI.
6.3.2.5 | PUCCH_encoder.m |
7.1.3 | PBCH_encoder.m |
7.1.4 | PBCH_encoder.m |
7.1.5 | PBCH_encoder.m |
7.3.1 | PDCCH_encoder.m | Only implements the zero padding of DCI formats, to increase their length to 12 bits.
7.3.2 | PDCCH_encoder.m |
7.3.3 | PDCCH_encoder.m |
7.3.4 | PDCCH_encoder.m |

Many thanks to my colleagues at [AccelerComm](http://www.accelercomm.com), who have spent lots of time double checking that this code matches the standard.

Have fun! Rob.

