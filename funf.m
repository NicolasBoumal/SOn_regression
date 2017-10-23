function [memory, f, gA, gB, hA, hB] = funf(memory, A, B, OA, OB)
% Part of cost function for DIGRESS --- see paper.
%
% OA, OB are 3D matrices with slices being skew symmetric matrices.
%
% Nicolas Boumal, Oct. 2017.

    import manopt.tools.*;
    
    assert(all(size(A) == size(B)));
    assert(ndims(A) <= 3);
    
    if ~isfield(memory, 'A_B')
        memory.A_B = A-B;
    end
    A_B = memory.A_B;
    
    if ~isfield(memory, 'f')
        memory.f = .5*squeeze(sum(sum(A_B.^2)));
    end
    f = memory.f;
    
    gA =  A_B;
    gB = -A_B;
    
    
    % Remember: there is no memory for OA, OB specific data.
    
    if nargin == 4
        assert(all(size(A) == size(OA)));
        HA = multiprod(A, OA);
        hA = HA;
    end
    
    if nargin == 5
        assert(all(size(A) == size(OA)));
        assert(all(size(B) == size(OB)));
        HA = multiprod(A, OA);
        HB = multiprod(B, OB);
        hA = HA-HB;
        hB = -hA;
    end

end
