function [block_sequence] = make_seq(block_pattern, N)
% This function will take a pattern of sequence with 0 and 1, and will
% scale it to the length of the chain.
% Example block_pattern = [0 1] and chain length is 10,
% will return block_seq = [0 1 0 1 0 1 0 1 0 1], if chain is less than specified
% sequence, it will return N elements of block_pattern

chain_sz = N(end);
block_sz = size(block_pattern,2);

if chain_sz<=block_sz
    block_sequence = block_pattern(1:chain_sz);
else
    num_block = (chain_sz/block_sz);
    current_block = block_pattern;   
    for zz=1:1:num_block
        current_block = cat(2,current_block,block_pattern);        
    end  
    block_sequence = current_block(1:chain_sz);
end

end
