"""NetworkX lane: bipartite RBP-receptor interaction graph analytics."""
import json
import networkx as nx
import pandas as pd

def main():
    pairs = pd.read_csv("data/processed/rbp_host_pairs.csv")
    panel = pd.read_csv("data/processed/host_receptor_panel.csv")
    G = nx.Graph()
    for _, r in pairs.iterrows():
        G.add_edge(f"phage:{r.accession}", f"host:{r.host}", kind="infects")
    B = nx.Graph()
    for _, r in panel.iterrows():
        B.add_node(f"receptor:{r.receptor}", species=r.species)
    comps = list(nx.connected_components(G))
    deg = dict(G.degree())
    host_deg = sorted(((k, v) for k, v in deg.items() if k.startswith("host:")),
                      key=lambda kv: -kv[1])
    phage_deg = [v for k, v in deg.items() if k.startswith("phage:")]
    json.dump({
        "n_phage_nodes": sum(1 for k in deg if k.startswith("phage:")),
        "n_host_nodes": sum(1 for k in deg if k.startswith("host:")),
        "n_edges": G.number_of_edges(),
        "n_connected_components": len(comps),
        "largest_component": max(len(c) for c in comps),
        "mean_hosts_per_phage": round(sum(phage_deg) / max(len(phage_deg), 1), 3),
        "host_degree_ranking": host_deg,
        "n_receptor_nodes": B.number_of_nodes()},
        open("results/networkx_corpus_graph.json", "w"), indent=2)
    print(len(comps), host_deg[:4])

if __name__ == "__main__":
    main()
