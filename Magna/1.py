import networkx as nx ,csv


def read_graph (path ,useKey = True):
    g = nx.Graph()
    dictkey = {}
    with open(path) as f:
        # content = f.readlines() ;    content = [x.strip() for x in content]
        for _ in range(4):
            next(f)
        for i,line in enumerate(f):
            if line[0] == '|' :
                dictkey[i] = line.strip()[2:-2]
                continue
            elif len(line.strip().split(' ')) == 4 :
                line = line.strip()[:-7].split(' ')
                if useKey :
                    g.add_edge(int(line[0]),int(line[1]))
                else:
                    g.add_edge(dictkey[ int(line[0]) ], dictkey[ int(line[1])])
            else:
                continue
    return g,dictkey

G1 , K1 = read_graph("data/synthetic_nets_known_node_mapping/0Krogan_2007_high.gw",False)
print(len(G1.nodes))
G1 , K1 = read_graph("data/synthetic_nets_known_node_mapping/low_confidence/0Krogan_2007_high+5e.gw",False)
print(len(G1.nodes))

def crossover():

    return 0
