import numpy as np
import pandas as pd
import json

def create_dist_table(file_name):
    '''
    This function returns distance matrix along with headers
    '''
    Table = []
    df = pd.read_csv(file_name,index_col=0)
    headers1 = []
    for col in df.columns:
        headers1.append(col)
    for i in range(0,len(headers1)):
        dist = []
        for j in range(0,i):
            dist.append(df[headers1[i]][headers1[j]])
        Table.append(dist)
    #Strip commas and spaces in headers and spaces
    headers = []
    for i in headers1:
        i=i.replace(",","")
        i=i.replace(" ","")
        i=i.replace("(" , "")
        i=i.replace(")" , "")
        headers.append(i)
    return Table,headers


def find_max_dist(Table):
    '''
    Given a table, find the minimum distance in that table
    '''
    min_i = 0
    min_j = 0
    min_dist = -1*float("inf")
    for i in range(len(Table)):
        for j in range(len(Table[i])):
            if Table[i][j] > min_dist:
                min_dist = Table[i][j]
                min_i = i
                min_j = j
    return min_i,min_j,min_dist

def alter_table(Table,min_1,min_2):
    i = 0
    while i < min_1:
        avg = float(Table[min_1][i] + Table[min_2][i])/2
        Table[min_1][i] = avg
        i = i+1
    i = min_1+1
    while i < min_2:
        avg = float(Table[i][min_1] + Table[min_2][i])/2
        Table[i][min_1] = avg
        i = i + 1
    i = min_2+1
    while i<len(Table):
        avg = float(Table[i][min_1] + Table[i][min_2])/2
        Table[i][min_1] = avg
        Table[i].pop(min_2)
        i = i+1        
    Table.pop(min_2)

def merge_rows(Table , min_1 , min_2,min_dist,headers,branch_lengths):
    '''
    Given a table, merges two rows based on the indices of the least distances
    '''
    #Change branch_lenghts dictionary
    branch_lengths[headers[min_1] + headers[min_2]] = 1000/(min_dist)
    alter_table(Table,min_1,min_2)
    # Newick tree format for TA to input it into tree visualizer
    headers[min_1] = '(' + headers[min_1] + ',' + headers[min_2] + ':' + str(branch_lengths[headers[min_1] + headers[min_2]]) + ')'
    headers.pop(min_2)
    #If len(header) is 1, then append the string with root
    if len(headers)==1:
        headers[0] = headers[0]+'root' 

def UPMGA(Table,headers,branch_lengths):
    '''
    UPGMA for x iterations where during iteration x, length of header becomes 1
    '''
    while len(headers)>1 and len(Table)>1:
        max_i , max_j , max_dist = find_max_dist(Table)
        if max_i < max_j:
            merge_rows(Table , max_i , max_j,max_dist,headers,branch_lengths)
        elif max_j < max_i:
            merge_rows(Table , max_j , max_i,max_dist,headers,branch_lengths)
        print(len(Table))
        if len(headers) == 1:
            break
    return headers

def test_model(headers):
    '''
    Test if there is only one header present
    '''
    print(headers)
    print(len)

if __name__ == "__main__":

    Table,headers=create_dist_table('Pdistance.csv')
    branch_lengths = {}
    UPMGA(Table,headers,branch_lengths)
    test_model(headers)
    text_file = open('a22.txt',"w")
    text_file.write(headers[0])
    text_file.close()




