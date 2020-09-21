import numpy as np
import pandas as pd

def extract_seq(file_name):
    '''
    Given a text file, return a dictionary of headers and DNA sequences
    '''
    df = open(file_name,"r")
    lines = df.readlines()
    seq = []
    header = []
    a = ""
    for line in lines:
        if line[0]=='>':
            if len(a)!=0:
                seq.append(a);
            b = ""
            for i in range(1,len(line)):
                if line[i]!='\n':
                    b = b + line[i]
            header.append(b)
            a=""
        if line[0]!='>':
            for i in range(len(line)):
                if line[i]!='\n' and line[i]!=" ":
                    a = a + line[i]
    seq.append(a)
    return seq,header

def find_hamming(seq1 , seq2):
    '''
    Given two DNA sequence, return the hamming distance in O(n)
    '''
    dist = 0
    red_num = 0
    for i in range(len(seq1)):
        if seq1[i]!=seq2[i]:
            dist = dist + 1
        elif seq1[i]=='-' and seq2[i]=='-':
            red_num = red_num + 1
    return dist,red_num

if __name__ == "__main__":
    seq1,header1 = extract_seq('Nucleotide_alignment.txt')
    dist_matrix = np.zeros((len(seq1) , len(seq1)))
    for i in range(len(seq1)):
        for j in range(len(seq1)):
            if i!=j:
                dist,red_num = (find_hamming(seq1[i] , seq1[j]))
                dist_matrix[i][j] = float(dist)/(len(seq1[i])-red_num)
    data = {}
    for i in range(len(seq1)):
        data[header1[i]] = dist_matrix[i]
    df = pd.DataFrame(data, columns = header1, index=header1)
    df.to_csv('Ndistance.csv')
    print(df)
