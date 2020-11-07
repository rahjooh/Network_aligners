import matplotlib.pyplot as plt
import networkx as nx

file = open("raw data/ppi_networks/ce.tab", "rb")
file.readline()
G1 = nx.read_edgelist(file)
print(nx.info(G1))
file.close()

file = open("raw data/ppi_networks/dm.tab", "rb")
file.readline()
G2 = nx.read_edgelist(file)
print(nx.info(G2))
file.close()

G=nx.Graph()

def isorank(Mat1,Mat2):