import networkx as nx
import heapq

# Example coverage data
coverage_data = {
    'edge_1': 5, 'edge_2': 37, 'edge_3': 4, 'edge_4': 10, 'edge_5': 6, 'edge_6': 7,
    'edge_7': 21, 'edge_8': 46, 'edge_9': 4, 'edge_10': 188, 'edge_11': 25, 'edge_12': 37,
    'edge_13': 5, 'edge_14': 20, 'edge_15': 8, 'edge_16': 6, 'edge_17': 19, 'edge_18': 48,
    'edge_19': 6, 'edge_20': 6, 'edge_21': 7, 'edge_22': 5, 'edge_23': 22, 'edge_24': 6,
    'edge_161': 7, 'edge_195': 30, 'edge_198': 37, 'edge_199': 29
}

# Nodes in the largest connected component as per your data
largest_component_contigs = [
    'edge_161', 'edge_2', 'edge_195', 'edge_4', 'edge_69', 'edge_198', 'edge_199',
    'edge_8', 'edge_10', 'edge_11', 'edge_12', 'edge_44', 'edge_14', 'edge_15',
    'edge_176', 'edge_18', 'edge_20', 'edge_23'
]

# Mapping contig names to integers for graph nodes
contig_to_node = {contig: i for i, contig in enumerate(largest_component_contigs)}
node_to_contig = {i: contig for contig, i in contig_to_node.items()}

# Creating the graph
G = nx.Graph()
for contig in largest_component_contigs:
    G.add_node(contig_to_node[contig])

# Add edges for demonstration, in a real scenario this should come from the GFA analysis
# For simplicity, let's assume each node is connected to the next, forming a simple cycle
for i in range(len(largest_component_contigs)):
    G.add_edge(contig_to_node[largest_component_contigs[i]],
               contig_to_node[largest_component_contigs[(i + 1) % len(largest_component_contigs)]])
def a_star_search(G, start_node, coverage_data, node_to_contig):
    def heuristic(node):
        contig = node_to_contig[node]
        return -coverage_data[contig]  # Negative for higher coverage priority

    open_set = set([start_node])
    came_from = {}
    g_score = {node: float('-inf') for node in G.nodes()}
    g_score[start_node] = coverage_data[node_to_contig[start_node]]

    open_heap = [(-g_score[start_node], start_node)]

    while open_set:
        _, current = heapq.heappop(open_heap)
        open_set.remove(current)
        print(f"open: {open_set}")

        if current == start_node and came_from:
            path = reconstruct_path(came_from, current)
            if len(path) > 1 and path[0] == start_node:
                return [node_to_contig[node] for node in path]

        for neighbor in G.neighbors(current):
            if neighbor == current or (neighbor == start_node and len(came_from) < len(G) - 1):
                continue  # Avoid trivial cycles

            tentative_g_score = g_score[current] + coverage_data[node_to_contig[neighbor]]

            if tentative_g_score > g_score[neighbor]:
                came_from[neighbor] = current
                g_score[neighbor] = tentative_g_score
                if neighbor not in open_set:
                    heapq.heappush(open_heap, (-g_score[neighbor], neighbor))
                    open_set.add(neighbor)

    return None

def reconstruct_path(came_from, start_node):
    path = [start_node]
    while start_node in came_from:
        start_node = came_from[start_node]
        path.append(start_node)
    return path[::-1]

# Find the start node (highest coverage)
start_contig = max(coverage_data, key=coverage_data.get)
start_node = contig_to_node[start_contig]

# Run A* search
circular_path = a_star_search(G, start_node, coverage_data, node_to_contig)
print(f'Circular path: {circular_path}')
