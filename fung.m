function [memory, g, gA, gB, gC, hA, hB, hC] = fung(memory, A, B, C, OA, OB, OC)
% Part of cost function for DIGRESS --- see paper.
%
% OA, OB, OC are 3D matrices with slices being skew symmetric matrices.
%
% Nicolas Boumal, Oct. 2017.
    
    assert(all(size(A) == size(B)));
    assert(all(size(A) == size(C)));
    assert(ndims(A) <= 3);
    
    if ~isfield(memory, 'BpC')
        memory.BpC = B+C;
    end
    BpC = memory.BpC;
    
    if ~isfield(memory, 'At')
        memory.At = multitransp(A);
    end
    At = memory.At;
   
    if ~isfield(memory, 'R')
        memory.R = multiskew(multiprod(At, BpC));
    end
    R = memory.R;
    
    if ~isfield(memory, 'g')
        memory.g = .5*squeeze(sum(sum(R.^2)));
    end
    g = memory.g;
    
    if ~isfield(memory, 'gA') || ~isfield(memory, 'gB')
        memory.gA = -multiprod(BpC, R);
        memory.gB =  multiprod(A, R);
    end
    gA = memory.gA;
    gB = memory.gB;
    gC = memory.gB;
    
    if nargin == 7
        
        assert(all(size(A) == size(OA)));
        assert(all(size(B) == size(OB)));
        assert(all(size(C) == size(OC)));
        
        HA = multiprod(A, OA);
        HB = multiprod(B, OB);
        HC = multiprod(C, OC);
        
        HBpHC = HB + HC;
        HAt = multitransp(HA);
        
        HR = multiskew( multiprod(HAt,  BpC) + multiprod( At, HBpHC) );
        
        hA = -(multiprod(HBpHC, R) + multiprod(BpC, HR));
        hB =   multiprod(HA, R) + multiprod(A, HR);
        hC =   hB;
        
    end

end
