import networkx as nx
import heapq
def reconstruct_path(came_from, current):
    """
    Reconstructs the path from start to goal given the came_from map and the current (goal) node.
    """
    total_path = [current]
    while current in came_from:
        current = came_from[current]
        total_path.insert(0, current)  # Insert at the beginning to reverse the path
    return total_path

# Initialize came_from dictionary for tracking the path
came_from = {}

def extract_coverage_from_gfa(gfa_path):
    """
    Extracts coverage information from a GFA file. This assumes coverage
    is stored in a tag like 'RC:i:' (read count) in segment lines.
    Returns a dictionary with segment names as keys and coverage as values.
    """
    gfa = gfapy.Gfa.from_file(gfa_path)
    coverage_data = {}

    for segment in gfa.segments:
        segment_name = segment.name
        # Try to access the 'RC:i' tag for read count
        if "RC" in segment.tagnames:
            read_count = segment.get("RC")
            coverage_data[segment_name] = read_count
        else:
            # Handle segments without the RC tag
            coverage_data[segment_name] = 0  # Or some default/placeholder value

    return coverage_data
# Assuming you've created the graph `G` and have start and end nodes defined
coverage_data = extract_coverage_from_gfa("assembly_2graph_combined.gfa")


def coverage_func(current, neighbor, coverage_data):
    """
    Defines the cost of moving from current to neighbor based on coverage.
    Lower coverage results in a higher cost to prioritize high-coverage paths.
    """
    # Lookup coverage for the neighbor; add 1 to avoid division by zero
    return 1.0 / (coverage_data.get(str(neighbor), 1) + 1)


# Modified A* Function with 'came_from' handling
def a_star(graph, start, end, coverage_func):
    open_set = set([start])
    closed_set = set()

    g_score = {node: float('inf') for node in graph.nodes}
    g_score[start] = 0

    f_score = {node: float('inf') for node in graph.nodes}
    f_score[start] = coverage_func(start)

    open_heap = [(f_score[start], start)]

    while open_set:
        current = heapq.heappop(open_heap)[1]
        if current == end:
            return reconstruct_path(came_from, current)

        open_set.remove(current)
        closed_set.add(current)

        for neighbor in graph.neighbors(current):
            if neighbor in closed_set:
                continue

            tentative_g_score = g_score[current] + 1 / coverage_func(current,
                                                                     neighbor)  # Adjust coverage_func as needed

            if neighbor not in open_set:
                open_set.add(neighbor)
                heapq.heappush(open_heap, (f_score[neighbor], neighbor))
            elif tentative_g_score >= g_score[neighbor]:
                continue

            came_from[neighbor] = current
            g_score[neighbor] = tentative_g_score
            f_score[neighbor] = g_score[neighbor] + coverage_func(neighbor)  # Adjust coverage_func as needed

    return None


# Assuming coverage_data is a dictionary {segment_name: coverage}
# and gfa is your Gfa object loaded with gfapy
longest_segment = max(gfa.segments, key=lambda x: len(x.sequence))
start_node = longest_segment.name

# Modify the A* call to include coverage_data
path = a_star(G, start_node, end_node, lambda current, neighbor: coverage_func(current, neighbor, coverage_data))

# print("Optimal path based on coverage:", path)

