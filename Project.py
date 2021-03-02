#!/usr/bin/env python
# coding: utf-8

# In[31]:


def dna_count(dna):
    dna = dna.upper()
    count_A = dna.count('A')
    count_C = dna.count('C')
    count_G = dna.count('G')
    count_T = dna.count('T')
    return count_A,count_C,count_G,count_T


# In[28]:


def dna2rna(dna):
    rna = ''
    rna = rna.upper()
    for symbol in dna:
        if symbol == 'A':
            rna = rna + 'U'
        elif symbol == 'T':
            rna = rna + 'A'
        elif symbol == 'C':
            rna = rna + 'G'
        elif symbol == 'G':
            rna = rna + 'C' 
    return rna


# In[26]:


def rna2dna(rna):
    dna = ''
    for symbol in rna:
        if symbol == 'U':
            dna = dna + 'A'
        elif symbol == 'A':
            dna = dna + 'T'
        elif symbol == 'C':
            dna = dna + 'G'
        elif symbol == 'G':
            dna = dna + 'C'
    return dna


# In[32]:


def rna2codon(triplet):
    genetic_code = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',

        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',

        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',

        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    allowed_codons = set('ACGU')
    triplet.upper()
    amino = genetic_code.get(triplet)
    if amino != None:
        return amino
    else:
        return 'Invalid'


# In[29]:


def rna2codons(triplets):
    amino = ''
    for i in range( 0,int( len( triplets ) / 3 ) ):
        amino += rna2codon(triplets[3*i:3*i+3])
    return amino


# In[30]:


def dna2codons(dnaString):
    # Convert the string from DNA to RNA. (Which function does this?)
    rnaString = dna2rna(dnaString)
    
    # Convert the RNA string to its corresponding protein expression string. (Which function does this?)
    
    codons = rna2codons(rnaString)
    
    # Return the resulting string.
    return codons


# In[ ]:




