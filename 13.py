# Given a collection of strings Dna and an integer d, a k-mer is a (k,d)-motif 
# if it appears in every string from Dna with at most d mismatches.
# Find all (k, d)-motifs in a collection of strings.

def MOTIFENUMERATION(Dna, k, d):
    def neighbors(pattern, d):
        # Generate the d-neighborhood of the pattern
        if d == 0:
            return {pattern}
        if len(pattern) == 0:
            return {''}
        
        nucleotides = ['A', 'C', 'G', 'T']
        neighborhood = set()
        suffix_neighbors = neighbors(pattern[1:], d)
        
        for text in suffix_neighbors:
            if HammingDistance(pattern[1:], text) < d:
                for nucleotide in nucleotides:
                    neighborhood.add(nucleotide + text)
            else:
                neighborhood.add(pattern[0] + text)
        
        return neighborhood
    
    def HammingDistance(p, q):
        # Compute the Hamming distance between strings p and q
        distance = 0
        for p_char, q_char in zip(p, q):
            if p_char != q_char:
                distance += 1
        return distance
    
    patterns = set()
    
    # Iterate over each k-mer in each DNA string
    for dna_str in Dna:
        for i in range(len(dna_str) - k + 1):
            pattern = dna_str[i:i+k]
            neighborhood = neighbors(pattern, d)
            
            # Check if each neighbor appears in all DNA strings with at most d mismatches
            for neighbor in neighborhood:
                if all(any(HammingDistance(neighbor, dna_str[j:j+k]) <= d for j 
                           in range(len(dna_str) - k + 1)) for dna_str in Dna):
                    patterns.add(neighbor)
    
    return patterns

# Example usage
with open('rosalind_ba2a.txt', 'r')as file:
    # Read the first line and split it into two numbers
    first_line = file.readline().strip()
    k, d = map(int, first_line.split())

    # Read the remaining lines into a set of strings
    Dna = set()
    for line in file:
        Dna.add(line.strip())   

print(MOTIFENUMERATION(Dna, k, d))
