function C = jacobiConstantCR3BP(state, mu)

    for i = 1:length(state)
        x  = state(i,1);
        y  = state(i,2);
        z  = state(i,3);
        xd = state(i,4);
        yd = state(i,5);
        zd = state(i,6);
        
        r1 = sqrt((x + mu)^2 + y^2 + z^2);
        r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);
    
    
        C(i) = x^2 + y^2 + 2*((1-mu)/r1 + mu/r2) - (xd^2 + yd^2 + zd^2);
    end

end