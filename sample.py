import pandas as pd

#find the longest common substring
def longest_common_substring(s1, s2):
    m, n = len(s1), len(s2)
    # Create a 2D array to store lengths of longest common suffixes
    lcsuff = [[0] * (n + 1) for _ in range(m + 1)]
    length = 0  # Length of longest common substring
    end_pos = 0  # Ending position of longest common substring in s1

    # Build the lcsuff array
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s1[i - 1] == s2[j - 1]:
                lcsuff[i][j] = lcsuff[i - 1][j - 1] + 1
                if lcsuff[i][j] > length:
                    length = lcsuff[i][j]
                    end_pos = i
            else:
                lcsuff[i][j] = 0

    # The longest common substring
    return s1[end_pos - length: end_pos]

#search IMGT-DB for protein sequence
def find_sequence(dataframe, name, species='Homo sapiens'):
    found = dataframe.loc[(dataframe['name'] == name) & (dataframe['species'] == species)]
    if len(found) > 0:
        return found['FASTA'].iloc[0]
    return ""

#get the full sequence from finding the longest common string in both v and j gene sequences
#return the replaced full sequence and the begin and end AA number(index)
def get_tcr_seq_cdr_region(gene_dataframe, v_gene, j_gene, cdr3, species="Homo sapiens"):
    v_seq = find_sequence(gene_dataframe, v_gene, species)
    j_seq = find_sequence(gene_dataframe, j_gene, species)
    #find the longest common string beween gene and cdr3
    v_common = longest_common_substring(cdr3, v_seq)
    j_common = longest_common_substring(cdr3, j_seq)
    if len(v_common) == 0 or len(j_common) == 0:
        raise AssertionError("Cannot find matching part")
    full_seq = ""
    begin = v_seq.index(v_common)
    end = begin + len(cdr3) - 1
    end_j = j_seq.index(j_common) + len(j_common)
    #assemble v gene (until the beginning of cdr3) + cdr3 + j gene (untile the end of cdr3)
    full_seq = v_seq[:begin] + cdr3 + j_seq[end_j:]
    return begin, end, full_seq

def main():
    v_b = "TRBV7-9*01"
    j_b = "TRBJ1-5*01"
    cdr3 = "CASSLAVGEGPQHF"
    df_gene = pd.read_csv('IMGT_GENE.csv', index_col=0)
    begin, end, full_seq = get_tcr_seq_cdr_region(df_gene, v_b, j_b, cdr3)
    print(f'CDR3 {cdr3} is at the region between {begin} and {end} AA in\n{full_seq}')
    
if __name__ == '__main__':
    main()
    