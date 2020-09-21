import numpy as np
import pandas as pd


def extract_seq(file_name):
    '''
    Given a text file, return a dictionary of headers and DNA sequences
    '''
    f = open(file_name,"r")
    lines = f.readlines()
    seq = []
    header = []
    a = ""
    for line in lines:
        if line[0]=='>':
            if len(a)!=0:
                seq.append(a)
            b = ""
            for i in range(1,len(line)):
                if line[i]!='\n' and line[i]!='\r':
                    b = b + line[i]
            header.append(b)
            a=""
        if line[0]!='>':
            for i in range(len(line)):
                if line[i]!='\n' and line[i]!=" " and line[i]!='\r':
                    a = a + line[i]
    seq.append(a)
    return seq,header

def get_penalty(file_name):
    '''
    Retrieves cost matrix as a dictionary from BLOSUM62
    '''
    df = open('BLOSUM62.txt')
    lines = df.readlines()
    penalty = {}
    seq_list = []
    for i in range(1,len(lines[0])):
        if lines[0][i]!=' ' and lines[0][i]!="\n":
            seq_list.append(lines[0][i])
    for i in range(1,len(lines)):
        bp1 = lines[i][0]
        z = 0
        j = 1
        while j<len(lines[i]):
            num=0
            str1 = ""
            if (lines[i][j]=='-' or (lines[i][j]!='\n' and lines[i][j]!='\r' and lines[i][j]!=" ")):
                while lines[i][j]!='\n' and lines[i][j]!=' ' and lines[i][j]!='\r':
                    str1 = str1 + lines[i][j]
                    j=j+1
                    if j>=len(lines[i]):
                        break
                penalty[bp1+seq_list[z]] = str1
                z = z+1
            j=j+1
    return penalty,seq_list

def find_distance(seq1,seq2,penalty):
    '''
    Given the penalty matrix and sequences; the function returns their distance
    '''
    distance = 0
    for i in range(len(seq1)):
        if seq1[i]=='-' and seq2[i]=='-':
            continue
        elif seq1[i] == '-' and seq2[i] != '-':
            print(i,seq1[i],seq1[i-1])
            if seq1[i-1] != '-' or (seq1[i-1]=='-' and seq2[i-1]=='-'): 
                distance = distance - 11
            elif seq1[i-1] == '-' :
                distance = distance - 1
        elif seq2[i] == '-' and seq1[i] != '-':
            if seq2[i-1] != '-' or (seq1[i-1]=='-' and seq2[i-1]=='-'):
                distance = distance - 11
            elif seq2[i-1] == '-' :
                distance = distance - 1
       

        elif seq1[i]!='-' and seq2[i]!='-':
            distance = distance + int(penalty[seq2[i] + seq1[i]])
            
        
    return distance
            

if __name__ == "__main__":

    seq2,header2 = extract_seq('Protein_alignment.txt')
    print("Number of headers : " , len(header2))
    print("Number of sequences: " , len(seq2))
    penalty,seq_list = get_penalty('BLOSUM62.txt')
    dist_matrix = np.zeros((len(seq2) , len(seq2)))
    for i in range(len(seq2)):
        for j in range(len(seq2)):
            dist_matrix[i][j] = find_distance(seq2[i],seq2[j],penalty)
    data = {}
    for i in range(len(seq2)):
        data[header2[i]] = dist_matrix[i]
    df = pd.DataFrame(data, columns = header2, index=header2)
    df.to_csv('Pdistance.csv')
    print(df)