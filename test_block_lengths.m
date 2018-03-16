for A = 12:19
    A
    bad = true;
    try
        PUCCH_encoder(round(rand(1,A)),A+8);
    catch
        bad = false;
    end
    if bad
        error('Rob');
    end
    
    PUCCH_encoder(round(rand(1,A)),A+9);
    
    PUCCH_encoder(round(rand(1,A)),8192);
        
    
     bad = true;
    try
        PUCCH_encoder(round(rand(1,A)),8193);
    catch
        bad = false;
    end
    if bad
        error('Rob');
    end

end

for A = 20:359
    A
    bad = true;
    try
        PUCCH_encoder(round(rand(1,A)),A+10);
    catch
        bad = false;
    end
    if bad
        error('Rob');
    end
    
    PUCCH_encoder(round(rand(1,A)),A+11);
    
    PUCCH_encoder(round(rand(1,A)),8192);
        
    
     bad = true;
    try
        PUCCH_encoder(round(rand(1,A)),8193);
    catch
        bad = false;
    end
    if bad
        error('Rob');
    end

end


for A = 360:1013
    A
    bad = true;
    try
        PUCCH_encoder(round(rand(1,A)),A+10);
    catch
        bad = false;
    end
    if bad
        [A,A+10]
        error('Rob');
    end
    
    PUCCH_encoder(round(rand(1,A)),A+11);
    
    PUCCH_encoder(round(rand(1,A)),16385);
        
    
     bad = true;
    try
        PUCCH_encoder(round(rand(1,A)),16386);
    catch
        bad = false;
    end
    if bad
        [A,16386]
        error('Rob');
    end

end


for A = 1014:1706
    A
    bad = true;
    try
        PUCCH_encoder(round(rand(1,A)),max(2*ceil(A/2)+22,1088)-1);
    catch
        bad = false;
    end
    if bad
        disp('Should not work but does');
        [A,max(2*ceil(A/2)+22,1088)-1]
        error('Rob');
    end
    
    bad = false;
    try
        PUCCH_encoder(round(rand(1,A)),max(2*ceil(A/2)+22,1088));
    catch
        bad = true;
    end
    if bad
        disp('Should work but does not');
        [A,max(2*ceil(A/2)+22,1088)]
        error('Rob');
    end
        
 
    
    bad = false;
    try
        PUCCH_encoder(round(rand(1,A)),16385);
    catch
        bad = true;
    end
    if bad
        disp('Should work but does not');
        [A,2*ceil(A/2)+22]
        error('Rob');
    end
    
    
        
    
     bad = true;
    try
        PUCCH_encoder(round(rand(1,A)),16386);
    catch
        bad = false;
    end
    if bad
        disp('Should not work but does');
        [A,16386]
        error('Rob');
    end

end