import numpy as np , pandas as pd ,os
from datetime import datetime


alpha = 0.9

def proteinSet(path):
    proteins = set()
    file = open(path, 'r')
    file.readline()
    for line in file:
        line = line.strip().split('\t')
        proteins.add(line[0])
        proteins.add(line[1])
    return proteins

def protein(path):
    df = pd.read_csv(path,sep='\t')
    pvector_v =np.sort( pd.concat([df['INTERACTOR_A'],df['INTERACTOR_B']]).unique() )

    pneighbour_df =df [['INTERACTOR_B','INTERACTOR_A']]
    pneighbour_df.columns = ['INTERACTOR_A','INTERACTOR_B']
    pneighbour_df = pd.concat([df ,pneighbour_df])
    pneighbour_df = pneighbour_df.groupby('INTERACTOR_A')['INTERACTOR_B'].apply(','.join ).reset_index()
    pneighbour_df.columns = ['protein','neighbours']
    return pvector_v,pneighbour_df

def NeighboursDict(path):
    Neighbours = {}
    file = open(path, 'r')
    file.readline()
    for line in file:
        line = line.strip().split('\t')

        if line[0] in Neighbours.keys():
            Neighbours[line[0]].add(line[1])
        else:
            Neighbours[line[0]]=set(line[1])

        if line[1] in Neighbours.keys():
            Neighbours[line[1]].add(line[0])
        else:
            Neighbours[line[1]]=set(line[0])
    return Neighbours

def blast(path , switch =False):
    blastDict = {}
    file = open(path, 'r')
    header = ()
    for line in file:
        line = line.strip().split('     ')
        if switch :
            blastDict[line[1] + '-' + line[0]] = float(line[2])
            header += (line[1] + '-' + line[0],)
        else:
            blastDict[line[0]+'-'+line[1]] = float(line[2])
            header += (line[0] + '-' + line[1],)

    blastVec = np.zeros(len(header))
    for i,pair in enumerate(header):
        blastVec[i] = blastDict[pair]

    normalizedBlastVec  = blastVec / np.linalg.norm(blastVec)
    return normalizedBlastVec #,blastVec,blastDict

def OptcalcStochasticMatrix1(network1path , network2path , printdetail = False):
    if os.path.exists("./preprocessed data/A[" + network1path[-6:-4] + '-' + network2path[-6:-4] + ".npy") :
        return np.load("./preprocessed data/A[" + network1path[-6:-4] + '-' + network2path[-6:-4] + ".npy")

    proteinsSet1 = proteinSet(network1path)
    NeighboursDict1 = NeighboursDict(network1path)

    proteinsSet2 = proteinSet(network2path)
    NeighboursDict2 = NeighboursDict(network2path)

    header = ()
    for p1 in proteinsSet1 :
        for p2 in proteinsSet2 :
            header = header + (p1+'_'+p2,)
            print(len(header) , '/' , len(proteinsSet1)*len(proteinsSet2))
    header = sorted(header)
    exit(-12)
    A = np.zeros((len(header),len(header)))

    for i,row_pair  in enumerate(header ):
        for j, col_pair in enumerate(header):
            row_pairs = row_pair.split('_')
            col_pairs = col_pair.split('_')
            if col_pairs[0] in NeighboursDict1[row_pairs[0]] and \
                col_pairs[1] in NeighboursDict2[row_pairs[1]] :
                A[header.index(row_pair),header.index(col_pair)]= 1 / (len(NeighboursDict1[col_pairs[0]]) * len (NeighboursDict2[col_pairs[1]]) )
                if printdetail :
                    print('pair (',row_pair,',',col_pair,') on N(',row_pairs[0],',',col_pairs[0],')=1 in G1 and',\
                                                              ' N(',row_pairs[1],',',col_pairs[1],')=1 in G2 )  => ' ,end='')
                    print('A[',header.index(row_pair),',',header.index(col_pair),'] = 1 / |N(',col_pairs[0],')| |N(',col_pairs[1],')| ' ,end='')
                    print( ' =1/',len(NeighboursDict1[col_pairs[0]]),'x',len (NeighboursDict2[col_pairs[1]]),' = ', (1 / (len(NeighboursDict1[col_pairs[0]]) * len (NeighboursDict2[col_pairs[1]]) ) ))

    if printdetail:
        print(A)
        print(A.sum(axis=0))
    np.save("preprocessed data/A["+network1path[-6:-4]+'-'+network2path[-6:-4]+".npy" , A)
    return A


def calcStochasticMatrix(network1path , network2path , printdetail = False):
    if os.path.exists("./preprocessed data/A[" + network1path[-6:-4] + '-' + network2path[-6:-4] + "].npy") :
        return np.load("./preprocessed data/A[" + network1path[-6:-4] + '-' + network2path[-6:-4] + "].npy")

    proteinsSet1 = proteinSet(network1path)
    NeighboursDict1 = NeighboursDict(network1path)

    proteinsSet2 = proteinSet(network2path)
    NeighboursDict2 = NeighboursDict(network2path)

    header = ()
    for p1 in proteinsSet1 :
        for p2 in proteinsSet2 :
            header = header + (p1+'_'+p2,)
    header = sorted(header)

    A = np.zeros((len(header),len(header)))

    for i,row_pair  in enumerate(header ):
        for j, col_pair in enumerate(header):
            row_pairs = row_pair.split('_')
            col_pairs = col_pair.split('_')
            if col_pairs[0] in NeighboursDict1[row_pairs[0]] and \
                col_pairs[1] in NeighboursDict2[row_pairs[1]] :
                A[header.index(row_pair),header.index(col_pair)]= 1 / (len(NeighboursDict1[col_pairs[0]]) * len (NeighboursDict2[col_pairs[1]]) )
                if printdetail :
                    print('pair (',row_pair,',',col_pair,') on N(',row_pairs[0],',',col_pairs[0],')=1 in G1 and',\
                                                              ' N(',row_pairs[1],',',col_pairs[1],')=1 in G2 )  => ' ,end='')
                    print('A[',header.index(row_pair),',',header.index(col_pair),'] = 1 / |N(',col_pairs[0],')| |N(',col_pairs[1],')| ' ,end='')
                    print( ' =1/',len(NeighboursDict1[col_pairs[0]]),'x',len (NeighboursDict2[col_pairs[1]]),' = ', (1 / (len(NeighboursDict1[col_pairs[0]]) * len (NeighboursDict2[col_pairs[1]]) ) ))
        print(i,'/',len(header))
    if printdetail:
        print(A)
        print(A.sum(axis=0))
    np.save("preprocessed data/A["+network1path[-6:-4]+'-'+network2path[-6:-4]+"].npy" , A)
    return A

def  calcSimilarityVector (network1path , network2path ,iteration =10 ,detail =False):
    if os.path.exists("./preprocessed data/R[" + network1path[-6:-4] + '-' + network2path[-6:-4] + ".npy") :
        return np.load("./preprocessed data/R[" + network1path[-6:-4] + '-' + network2path[-6:-4] + ".npy")
    proteinsSet1 = proteinSet(network1path)
    NeighboursDict1 = NeighboursDict(network1path)

    proteinsSet2 = proteinSet(network2path)
    NeighboursDict2 = NeighboursDict(network2path)

    header = ()
    for p1 in proteinsSet1 :
        for p2 in proteinsSet2 :
            header = header + (p1+'_'+p2,)
    header = sorted(header)

    similarityDict = {}
    total = len(proteinsSet1) + len(proteinsSet2)
    for p1 in proteinsSet1 :
        for p2 in proteinsSet2 :
            similarityDict[p1+'_'+p2] = 1/total

    for iter in range(iteration):
        for key in similarityDict.keys():
            keys = key.split('_')
            sum = 0

            for neighbour1 in NeighboursDict1[keys[0]]:
                for neighbour2 in NeighboursDict2[keys[1]]:
                    sum += (1 / (len(NeighboursDict1[neighbour1]) * len(NeighboursDict2[neighbour2]))) * similarityDict[neighbour1 + '_' + neighbour2]
                    # print('similarityDict[',key ,']=',sum , ' nei(',neighbour1,') = ',NeighboursDict1[neighbour1], ' nei(',neighbour2,') = ',NeighboursDict2[neighbour2])

            similarityDict[key] = sum
        if detail :print(iter,similarityDict)

    similarityVec = np.zeros(len(header))
    for i,pair in enumerate(header) :
        similarityVec[i] = similarityDict[pair]

    # key1 = ''
    # value1 = 0
    # for neighbour2 in NeighboursDict2[keys[1]]:
    #     if R[neighbour1+'_'+neighbour2] > value1:
    #         value1 = R[neighbour1+'_'+neighbour2]
    #         key1 = neighbour1+'_'+neighbour2
    #     R[neighbour1+'_'+neighbour2] = 1/total
    # R[key1] = value1
    # print('key:',key , 'value:',value)
    # print(R[neighbour1+'_A'],R[neighbour1+'_A'],R[neighbour1+'_B'],R[neighbour1+'_C'],R[neighbour1+'_D'],R[neighbour1+'_E'])
    np.save("preprocessed data/R[" + network1path[-6:-4] + '-' + network2path[-6:-4] + ".npy", similarityVec)
    return similarityVec  #,similarityDict

# A = calcStochasticMatrix('raw data/ppi_networks/ce.tab','raw data/ppi_networks/dm.tab',True)

# R = calcSimilarityVector('raw data/ppi_networks/ce.tab','raw data/ppi_networks/dm.tab',10)
# E = blast('raw data/BLAST_Bit_Scores/ce-dm.evals')

# res = alpha * A.dot(R) + (1-alpha) * E
# print(res)

#
#
#
# start=datetime.now()
#
# a = proteinSet('raw data/ppi_networks/ce.tab')
# a = np.array(a)
#
# print (datetime.now()-start)
#
# start=datetime.now()
#
# b = NeighboursDict('raw data/ppi_networks/ce.tab')
# print(b)
# b = pd.DataFrame(b.items())
#
# print (datetime.now()-start)
#
# start=datetime.now()
#
# c,d =protein('raw data/ppi_networks/ce.tab')
# print(b)
# print(d)
# print (datetime.now()-start)

# A = calcStochasticMatrix('raw data/ppi_networks/aa.tab','raw data/ppi_networks/bb.tab')
A = calcStochasticMatrix('raw data/ppi_networks/ce.tab','raw data/ppi_networks/dm.tab')