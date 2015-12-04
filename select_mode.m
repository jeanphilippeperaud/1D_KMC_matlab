function index = select_mode(cumul,Nmodes)
%SELECT_MODE - This function randomly draws an index within a cumulative
%distribution. The distribution of the indexes follows the probability
%density function from which the cumulative distribution is generated
%   The cumulative distribution cumul must be a 1D array with increasing
%   values. Nmodes is the length of the array.

    R = rand();
    i1 = 0;
    i3 = Nmodes;
    i2 = floor((i1+i3)/2);
    
    while (i3-i1>1)
        
        if R<cumul(i2)/cumul(Nmodes)
            i3  = i2;
            i2 = floor((i1+i3)/2);
        else
            i1 = i2;
            i2 = floor((i1+i3)/2);
        end
    end
    index = i3;        
    
end

