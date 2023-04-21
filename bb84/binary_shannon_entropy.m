function H_2 = binary_shannon_entropy(x)
    H_2 = -log2(x)*x - (1-x)*log2(1-x);
end