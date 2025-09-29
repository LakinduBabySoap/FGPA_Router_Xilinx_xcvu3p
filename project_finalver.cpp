#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <climits> // For INT_MAX and INT_MIN
#include <chrono>  // For timing functionality
#include <omp.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstring>
// Define structures and types
auto program_start_time = std::chrono::high_resolution_clock::now();
struct NodeInfo
{
    std::string type;
    int length;
    std::pair<int, int> begin_coords;
    std::pair<int, int> end_coords;
    std::string name;
};

struct Net
{
    int id;
    std::string name;
    int source;
    std::vector<int> sinks;
};

// Define typedefs to make code more readable
typedef std::unordered_map<int, NodeInfo> NodeMap;
typedef std::unordered_map<int, std::vector<int>> AdjacencyList;
// Add this near the other function declarations
struct DeviceData
{
    NodeMap nodes;
    AdjacencyList adj;
    AdjacencyList reverse_adj;
};
DeviceData parse_device_file(const std::string &device_path, const std::vector<Net> &nets);

std::vector<Net> parse_netlist_file(const std::string &netlist_path)
{
    std::vector<Net> nets;
    std::ifstream file(netlist_path);
    if (!file.is_open())
    {
        std::cerr << "Error: Netlist file " << netlist_path << " not found!" << std::endl;
        return nets;
    }

    int num_nets;
    file >> num_nets;
    file.ignore(); // ignore the newline

    for (int line_num = 1; line_num <= num_nets; line_num++)
    {
        std::string line;
        std::getline(file, line);
        std::istringstream iss(line);
        Net net;
        iss >> net.id >> net.name >> net.source;

        int sink;
        while (iss >> sink)
        {
            net.sinks.push_back(sink);
        }

        nets.push_back(net);
    }

    return nets;
}

int heuristic(int node1, int node2, const NodeMap &nodes)
{
    int x1 = nodes.at(node1).begin_coords.first;
    int y1 = nodes.at(node1).begin_coords.second;
    int x2 = nodes.at(node2).begin_coords.first;
    int y2 = nodes.at(node2).begin_coords.second;
    return std::abs(x1 - x2) + std::abs(y1 - y2);
}

std::vector<std::pair<int, int>> bidirectional_a_star(const AdjacencyList &adj, const AdjacencyList &reverse_adj, int source, int sink, const NodeMap &nodes, const std::unordered_set<int> &used_nodes)
{
    // Progress update

    // Priority queue for A* algorithm
    struct Node
    {
        int f_score;
        int g_score;
        int id;
        bool operator>(const Node &other) const
        {
            return f_score > other.f_score;
        }
    };

    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> forward_open_set;
    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> backward_open_set;

    std::unordered_map<int, bool> forward_closed;
    std::unordered_map<int, bool> backward_closed;

    std::unordered_map<int, int> forward_g_score;
    std::unordered_map<int, int> backward_g_score;

    std::unordered_map<int, int> forward_came_from;
    std::unordered_map<int, int> backward_came_from;

    forward_g_score[source] = 0;
    backward_g_score[sink] = 0;

    forward_open_set.push({heuristic(source, sink, nodes), 0, source});
    backward_open_set.push({heuristic(sink, source, nodes), 0, sink});

    int meeting_point = -1;
    int best_path_length = INT_MAX;
    int best_path_found_at_iteration = 0;
    const int EARLY_TERMINATION_THRESHOLD = 800; // Terminate if no improvement for 50000 iterations
    int iterations = 0;
    const int MAX_ITERATIONS = 10000; // Prevent infinite loops
    const int MAX_BEAM_WIDTH = 500;   // Maximum nodes to consider at each step
    while (!forward_open_set.empty() && !backward_open_set.empty() && iterations < MAX_ITERATIONS)
    {
        iterations++;

        // Forward search step
        if (!forward_open_set.empty())
        {
            Node current = forward_open_set.top();
            forward_open_set.pop();

            // Skip if already processed
            if (forward_closed.find(current.id) != forward_closed.end())
            {
                continue;
            }

            forward_closed[current.id] = true;

            // Check if we found a meeting point
            if (backward_closed.find(current.id) != backward_closed.end())
            {
                int path_length = forward_g_score[current.id] + backward_g_score[current.id];
                if (path_length < best_path_length)
                {
                    best_path_length = path_length;
                    meeting_point = current.id;
                    best_path_found_at_iteration = iterations;
                }
            }
            int nodes_expanded = 0;
            // Expand neighbors
            auto it = adj.find(current.id);
            if (it != adj.end())
            {
                for (int neighbor : it->second)
                {
                    // Limit beam width
                    if (nodes_expanded++ >= MAX_BEAM_WIDTH)
                        break;
                    // Skip used nodes
                    if (used_nodes.find(neighbor) != used_nodes.end() &&
                        neighbor != sink)
                    {
                        continue;
                    }

                    int tentative_g = current.g_score + 1;
                    if (forward_g_score.find(neighbor) == forward_g_score.end() ||
                        tentative_g < forward_g_score[neighbor])
                    {

                        forward_came_from[neighbor] = current.id;
                        forward_g_score[neighbor] = tentative_g;
                        int f_score = tentative_g + heuristic(neighbor, sink, nodes);
                        forward_open_set.push({f_score, tentative_g, neighbor});
                    }
                }
            }
        }
        // Backward search step
        if (!backward_open_set.empty())
        {
            Node current = backward_open_set.top();
            backward_open_set.pop();
            int nodes_expanded = 0;
            // Skip if already processed
            if (backward_closed.find(current.id) != backward_closed.end())
            {
                continue;
            }

            backward_closed[current.id] = true;

            // Check if we found a meeting point
            if (forward_closed.find(current.id) != forward_closed.end())
            {
                int path_length = forward_g_score[current.id] + backward_g_score[current.id];
                if (path_length < best_path_length)
                {
                    best_path_length = path_length;
                    meeting_point = current.id;
                    best_path_found_at_iteration = iterations;
                }
            }
            // For backward search, we need nodes that have this node as a child
            auto it = reverse_adj.find(current.id);
            if (it != reverse_adj.end())
            {
                for (int parent : it->second)
                {
                    // Limit beam width
                    if (nodes_expanded++ >= MAX_BEAM_WIDTH)
                        break;
                    // Skip used nodes
                    if (used_nodes.find(parent) != used_nodes.end() &&
                        parent != source)
                    {
                        continue;
                    }

                    int tentative_g = current.g_score + 1;
                    if (backward_g_score.find(parent) == backward_g_score.end() ||
                        tentative_g < backward_g_score[parent])
                    {
                        backward_came_from[parent] = current.id;
                        backward_g_score[parent] = tentative_g;
                        int f_score = tentative_g + heuristic(parent, source, nodes);
                        backward_open_set.push({f_score, tentative_g, parent});
                    }
                }
            }
        }
        // Check for early termination if we've found a path but haven't improved in a while
        if (meeting_point != -1 && iterations - best_path_found_at_iteration > EARLY_TERMINATION_THRESHOLD)
        {
            // std::cout << "    Early termination at iteration " << iterations << " (no improvement for " << EARLY_TERMINATION_THRESHOLD << " iterations)" << std::endl;
            break;
        }
        // Provide occasional progress updates
        if (iterations % 50000 == 0)
        {
            // std::cout << "    A* search iterations: " << iterations << std::endl;
        }
    }
    // Path reconstruction
    std::vector<std::pair<int, int>> path;

    if (meeting_point != -1)
    {
        // std::cout << "    Path found at meeting point: " << meeting_point << std::endl;

        // Forward path
        std::vector<int> forward_path;
        int current = meeting_point;
        while (forward_came_from.find(current) != forward_came_from.end())
        {
            forward_path.push_back(current);
            current = forward_came_from[current];
        }
        std::reverse(forward_path.begin(), forward_path.end());

        // Connect source to first node
        if (!forward_path.empty())
        {
            path.push_back({source, forward_path[0]});

            // Connect remaining nodes in forward path
            for (size_t i = 0; i < forward_path.size() - 1; i++)
            {
                path.push_back({forward_path[i], forward_path[i + 1]});
            }
        }
        // Backward path
        current = meeting_point;
        while (backward_came_from.find(current) != backward_came_from.end())
        {
            int next = backward_came_from[current];
            path.push_back({current, next});
            current = next;
        }
    }
    else if (iterations >= MAX_ITERATIONS)
    {
        // std::cout << "    A* search reached maximum iterations (" << MAX_ITERATIONS << ")" << std::endl;
    }
    else
    {
        // std::cout << "    No path found between source and sink" << std::endl;
    }

    return path;
}
std::vector<std::vector<std::pair<int, int>>> generate_routing(const std::vector<Net> &nets, const AdjacencyList &adj, const AdjacencyList &reverse_adj, const NodeMap &nodes)
{
    std::vector<std::vector<std::pair<int, int>>> all_paths(nets.size());
    std::unordered_set<int> used_nodes;
    int total_nets = nets.size();

    // std::cout << "\n====== Starting Routing Process ======" << std::endl;
    // std::cout << "Total nets to route: " << total_nets << std::endl;

    // First, exclude source and sink nodes from being used by other nets
    for (const Net &net : nets)
    {
        used_nodes.insert(net.source);
        for (int sink : net.sinks)
        {
            used_nodes.insert(sink);
        }
    }

    // Attempt routing each net
    for (size_t net_idx = 0; net_idx < nets.size(); net_idx++)
    {
        const Net &net = nets[net_idx];

        // std::cout << "\n[Net " << (net_idx + 1) << "/" << total_nets << "] Routing net ID "
        //           << net.id << " (" << net.name << ")" << std::endl;
        // std::cout << "  Source: " << net.source << ", Sinks: " << net.sinks.size() << std::endl;

        // Track nodes used by this net
        std::unordered_set<int> this_net_nodes;
        this_net_nodes.insert(net.source);

        bool all_sinks_routed = true;
        for (size_t sink_idx = 0; sink_idx < net.sinks.size(); sink_idx++)
        {
            int sink = net.sinks[sink_idx];
            // std::cout << "  [Sink " << (sink_idx + 1) << "/" << net.sinks.size() << "] Routing from "
            //          << net.source << " to " << sink << "..." << std::endl;

            // Get temp nodes that should be excluded for this path finding
            std::unordered_set<int> temp_used = used_nodes;
            for (int node : this_net_nodes)
            {
                temp_used.erase(node); // Allow reuse of nodes from this net
            }

            auto path = bidirectional_a_star(adj, reverse_adj, net.source, sink, nodes, temp_used);

            if (path.empty())
            {
                // std::cout << "    Failed to find path!" << std::endl;
                all_sinks_routed = false;
                continue;
            }

            // std::cout << "    Found path with " << path.size() << " segments" << std::endl;
            //  Add path to results and update used nodes
            for (const auto &p_c : path)
            {
                all_paths[net_idx].push_back(p_c);
                this_net_nodes.insert(p_c.first);
                this_net_nodes.insert(p_c.second);
            }
        }

        if (!all_sinks_routed)
        {
            // std::cout << "  WARNING: Not all sinks could be routed for net " << net.id << std::endl;
        }
        else
        {
            // std::cout << "  Successfully routed all sinks for net " << net.id << std::endl;
        }

        // Add nodes used by this net to the global used set
        for (int node : this_net_nodes)
        {
            if (node != net.source && std::find(net.sinks.begin(), net.sinks.end(), node) == net.sinks.end())
            {
                used_nodes.insert(node);
            }
        }
    }

    std::cout << "\n====== Routing Completed ======" << std::endl;
    return all_paths;
}
void write_output(const std::vector<Net> &nets, const std::vector<std::vector<std::pair<int, int>>> &paths, const std::string &output_path)
{
    std::ofstream file(output_path);
    if (!file.is_open())
    {
        std::cerr << "Error occurred while writing to output file." << std::endl;
        return;
    }

    // For each net
    for (size_t i = 0; i < nets.size(); i++)
    {
        // Write net ID and name
        file << nets[i].id << " " << nets[i].name << std::endl;

        // Write each parent-child pair on a new line
        for (const auto &pip : paths[i])
        {
            file << pip.first << " " << pip.second << std::endl;
        }

        // Empty line between nets
        file << std::endl;
    }

    // std::cout << "Routing results successfully generated to: " << output_path << std::endl;
}
DeviceData parse_device_file(const std::string &device_path, const std::vector<Net> &nets)
{
    NodeMap nodes;
    AdjacencyList adj;
    AdjacencyList reverse_adj;
    std::unordered_set<int> required_nodes;
    int min_x = INT_MAX, min_y = INT_MAX, max_x = INT_MIN, max_y = INT_MIN;
    // Collect required nodes
    for (const Net &net : nets)
    {
        required_nodes.insert(net.source);
        for (int sink : net.sinks)
            required_nodes.insert(sink);
    }

    // std::cout << "Starting high-performance single-pass parsing..." << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();

    // First, read the entire file into memory
    std::ifstream file(device_path, std::ios::binary | std::ios::ate);
    if (!file.is_open())
    {
        std::cerr << "Error: Device file " << device_path << " not found!" << std::endl;
        return {nodes, adj, reverse_adj};
    }

    // Get file size and allocate buffer
    std::streamsize file_size = file.tellg();
    file.seekg(0, std::ios::beg);

    // std::cout << "File size: " << (file_size / 1024 / 1024) << " MB" << std::endl;

    std::vector<char> buffer(file_size);
    if (!file.read(buffer.data(), file_size))
    {
        std::cerr << "Error reading file into memory" << std::endl;
        return {nodes, adj, reverse_adj};
    }
    // Convert to string for easier processing
    std::string file_content(buffer.begin(), buffer.end());
    buffer.clear(); // Free memory from buffer

    auto read_time = std::chrono::high_resolution_clock::now();
    std::cout << "File loaded into memory in " << std::chrono::duration_cast<std::chrono::milliseconds>(read_time - start_time).count() << " ms" << std::endl;

    // Split content into lines
    std::vector<std::string> lines;
    size_t pos = 0;
    size_t prev = 0;

    // Reserve space based on estimated line count
    lines.reserve(file_size / 50); // Average line length estimate

    while ((pos = file_content.find('\n', prev)) != std::string::npos)
    {
        lines.push_back(file_content.substr(prev, pos - prev));
        prev = pos + 1;
    }

    // Handle the last line
    if (prev < file_content.size())
    {
        lines.push_back(file_content.substr(prev));
    }
    // Free memory from the entire file content
    file_content.clear();

    auto split_time = std::chrono::high_resolution_clock::now();
    // std::cout << "File split into " << lines.size() << " lines in " << std::chrono::duration_cast<std::chrono::milliseconds>(split_time - read_time).count() << " ms" << std::endl;

    // Read node count from first line
    int num_nodes = std::stoi(lines[0]);

    // Pre-allocate node storage
    nodes.reserve(num_nodes / 5); // Estimate that we'll need about 20% of nodes

    // First pass: Find required nodes and calculate bounding box
    // std::cout << "Processing required nodes and calculating bounding box..." << std::endl;

#pragma omp parallel sections shared(min_x, min_y, max_x, max_y)
    {
#pragma omp section
        {
            for (int i = 1; i <= num_nodes; i++)
            {
                if (i >= lines.size())
                    continue;

                const std::string &line = lines[i];
                // Quick check for node ID
                size_t space_pos = line.find(' ');
                if (space_pos == std::string::npos)
                    continue;

                int node_id = std::stoi(line.substr(0, space_pos));

                if (required_nodes.find(node_id) != required_nodes.end())
                {
                    // Parse full node data using direct string operations (faster than streams)
                    size_t pos1 = line.find(' ', space_pos + 1);
                    size_t pos2 = line.find(' ', pos1 + 1);
                    size_t pos3 = line.find(' ', pos2 + 1);
                    size_t pos4 = line.find(' ', pos3 + 1);
                    size_t pos5 = line.find(' ', pos4 + 1);
                    size_t pos6 = line.find(' ', pos5 + 1);

                    if (pos6 == std::string::npos)
                        continue;

                    std::string node_type = line.substr(space_pos + 1, pos1 - space_pos - 1);
                    int node_length = std::stoi(line.substr(pos1 + 1, pos2 - pos1 - 1));
                    int begin_x = std::stoi(line.substr(pos2 + 1, pos3 - pos2 - 1));
                    int begin_y = std::stoi(line.substr(pos3 + 1, pos4 - pos3 - 1));
                    int end_x = std::stoi(line.substr(pos4 + 1, pos5 - pos4 - 1));
                    int end_y = std::stoi(line.substr(pos5 + 1, pos6 - pos5 - 1));
                    std::string node_name = line.substr(pos6 + 1);

#pragma omp critical
                    {
                        nodes[node_id] = {node_type, node_length, {begin_x, begin_y}, {end_x, end_y}, node_name};

                        // Update bounding box
                        min_x = std::min(min_x, begin_x);
                        min_y = std::min(min_y, begin_y);
                        max_x = std::max(max_x, end_x);
                        max_y = std::max(max_y, end_y);
                    }
                }
            }
        }
    }
    auto required_time = std::chrono::high_resolution_clock::now();
    // std::cout << "Processed " << required_nodes.size() << " required nodes in " << std::chrono::duration_cast<std::chrono::milliseconds>(required_time - split_time).count() << " ms" << std::endl;

    // Calculate expanded bounding box
    const float EXPANSION_FACTOR = 1.3;
    int box_width = max_x - min_x;
    int box_height = max_y - min_y;

    min_x -= static_cast<int>(box_width * (EXPANSION_FACTOR - 1) / 2);
    min_y -= static_cast<int>(box_height * (EXPANSION_FACTOR - 1) / 2);
    max_x += static_cast<int>(box_width * (EXPANSION_FACTOR - 1) / 2);
    max_y += static_cast<int>(box_height * (EXPANSION_FACTOR - 1) / 2);

    // std::cout << "Bounding box: (" << min_x << "," << min_y << ") to (" << max_x << "," << max_y << ")" << std::endl;

    // Second pass: Process nodes within bounding box
    // std::cout << "Loading nodes within bounding box..." << std::endl;
    int nodes_in_box = 0;

#pragma omp parallel for reduction(+ : nodes_in_box) shared(nodes, min_x, min_y, max_x, max_y)
    for (int i = 1; i <= num_nodes; i++)
    {
        if (i >= lines.size())
            continue;

        const std::string &line = lines[i];
        // Quick check for node ID
        size_t space_pos = line.find(' ');
        if (space_pos == std::string::npos)
            continue;

        int node_id = std::stoi(line.substr(0, space_pos));

        // Skip if we already have this node
        if (nodes.find(node_id) != nodes.end())
        {
            continue;
        }
        // Parse coordinates for bounding box check using direct substring operations
        size_t pos1 = line.find(' ', space_pos + 1);
        size_t pos2 = line.find(' ', pos1 + 1);
        size_t pos3 = line.find(' ', pos2 + 1);
        size_t pos4 = line.find(' ', pos3 + 1);
        size_t pos5 = line.find(' ', pos4 + 1);

        if (pos5 == std::string::npos)
            continue;

        std::string node_type = line.substr(space_pos + 1, pos1 - space_pos - 1);

        // Skip certain node types
        if (node_type == "NULL" || node_type == "CONFIG")
        {
            continue;
        }

        int begin_x = std::stoi(line.substr(pos2 + 1, pos3 - pos2 - 1));
        int begin_y = std::stoi(line.substr(pos3 + 1, pos4 - pos3 - 1));

        // Check if in bounding box
        if (begin_x >= min_x && begin_x <= max_x && begin_y >= min_y && begin_y <= max_y)
        {
            int node_length = std::stoi(line.substr(pos1 + 1, pos2 - pos1 - 1));
            int end_x = std::stoi(line.substr(pos4 + 1, pos5 - pos4 - 1));
            int end_y;
            std::string node_name;

            size_t pos6 = line.find(' ', pos5 + 1);
            if (pos6 == std::string::npos)
            {
                end_y = std::stoi(line.substr(pos5 + 1));
                node_name = "";
            }
            else
            {
                end_y = std::stoi(line.substr(pos5 + 1, pos6 - pos5 - 1));
                node_name = line.substr(pos6 + 1);
            }

#pragma omp critical
            {
                nodes[node_id] = {node_type, node_length, {begin_x, begin_y}, {end_x, end_y}, node_name};
            }
            nodes_in_box++;
        }
    }

    auto bbox_time = std::chrono::high_resolution_clock::now();
    // std::cout << "Loaded " << nodes_in_box << " additional nodes in bounding box in " << std::chrono::duration_cast<std::chrono::milliseconds>(bbox_time - required_time).count() << " ms" << std::endl;
    // std::cout << "Total nodes: " << nodes.size() << std::endl;

    // Find where edge section starts
    // Find where edge section starts
    int edge_start_line = num_nodes + 1;
    while (edge_start_line < lines.size() && lines[edge_start_line].empty())
    {
        edge_start_line++;
    }

    // Create a node set for faster lookups
    std::unordered_set<int> node_set;
    for (const auto &entry : nodes)
    {
        node_set.insert(entry.first);
    }

    // Calculate file offset to edge section
    size_t offset = 0;
    file.close();                             // Close the file before reopening for memory mapping
    file.open(device_path, std::ios::binary); // Reopen to calculate offset
    for (int i = 0; i < edge_start_line; i++)
    {
        std::string line;
        std::getline(file, line);
        offset += line.length() + 1; // +1 for newline
    }
    file.close();

    // std::cout << "Edge section offset: " << offset << " bytes" << std::endl;

    // Memory map the edge section
    int fd = open(device_path.c_str(), O_RDONLY);
    if (fd == -1)
    {
        std::cerr << "Failed to open file for memory mapping" << std::endl;
        return {nodes, adj, reverse_adj};
    }

    struct stat sb;
    if (fstat(fd, &sb) == -1)
    {
        std::cerr << "Failed to get file size" << std::endl;
        close(fd);
        return {nodes, adj, reverse_adj};
    }

    // Map just the edge section

    // Process edges in parallel
    // std::cout << "Loading edges with memory mapping in parallel..." << std::endl;
    // Make sure offset is page-aligned
    size_t page_size = sysconf(_SC_PAGE_SIZE);
    off_t aligned_offset = (offset / page_size) * page_size;
    size_t offset_adjustment = offset - aligned_offset;
    size_t map_size = sb.st_size - aligned_offset; // Map from aligned position
    // Print debugging info
    // std::cout << "Attempting to map " << (map_size / 1024 / 1024) << " MB at offset " << (offset / 1024 / 1024) << " MB" << std::endl;

    // Call mmap with aligned offset
    char *mapped = (char *)mmap(NULL, map_size, PROT_READ, MAP_PRIVATE, fd, aligned_offset);

    if (mapped == MAP_FAILED)
    {
        std::cerr << "Memory mapping failed: " << strerror(errno) << " (errno: " << errno << ")" << std::endl;

        // Fallback to traditional file reading if mmap fails
        // std::cout << "Falling back to traditional file reading for edges..." << std::endl;

        // Implementation of fallback method...

        close(fd);
        return {nodes, adj, reverse_adj};
    }

    // Adjust mapped pointer to account for alignment
    mapped += offset_adjustment;
    int edges_loaded = 0;
    auto edge_start_time = std::chrono::high_resolution_clock::now();

#pragma omp parallel
    {
        // Thread-local adjacency lists
        AdjacencyList local_adj;
        AdjacencyList local_reverse_adj;
        int local_edges_loaded = 0;

        // Calculate chunk sizes for parallel processing
        int thread_id = omp_get_thread_num();
        int num_threads = omp_get_num_threads();
        size_t chunk_size = map_size / num_threads;
        size_t start_pos = thread_id * chunk_size;
        size_t end_pos = (thread_id == num_threads - 1) ? map_size : start_pos + chunk_size;

        // Find line start (skip partial line if not thread 0)
        if (thread_id > 0)
        {
            while (start_pos < end_pos && mapped[start_pos] != '\n')
                start_pos++;

            if (start_pos < end_pos)
                start_pos++; // Skip the newline
        }

        // Process lines in this thread's chunk
        size_t pos = start_pos;
        while (pos < end_pos)
        {
            // Find end of line
            size_t line_end = pos;
            while (line_end < end_pos && mapped[line_end] != '\n')
                line_end++;

            // Process the line if non-empty
            if (line_end > pos)
            {
                std::string line(mapped + pos, line_end - pos);

                // Skip empty lines
                if (!line.empty())
                {
                    size_t space_pos = line.find(' ');
                    if (space_pos != std::string::npos)
                    {
                        int parent = std::stoi(line.substr(0, space_pos));

                        if (node_set.find(parent) != node_set.end())
                        {
                            // Process all children on this line
                            size_t start = space_pos + 1;
                            size_t next;

                            while ((next = line.find(' ', start)) != std::string::npos)
                            {
                                int child = std::stoi(line.substr(start, next - start));
                                if (node_set.find(child) != node_set.end())
                                {
                                    local_adj[parent].push_back(child);
                                    local_reverse_adj[child].push_back(parent);
                                    local_edges_loaded++;
                                }
                                start = next + 1;
                            }

                            // Last child
                            if (start < line.length())
                            {
                                int child = std::stoi(line.substr(start));
                                if (node_set.find(child) != node_set.end())
                                {
                                    local_adj[parent].push_back(child);
                                    local_reverse_adj[child].push_back(parent);
                                    local_edges_loaded++;
                                }
                            }
                        }
                    }
                }
            }

            // Move to next line
            pos = line_end + 1;
        }

// Merge thread-local results into global data
#pragma omp critical
        {
            for (const auto &entry : local_adj)
            {
                adj[entry.first].insert(adj[entry.first].end(),
                                        entry.second.begin(), entry.second.end());
            }

            for (const auto &entry : local_reverse_adj)
            {
                reverse_adj[entry.first].insert(reverse_adj[entry.first].end(),
                                                entry.second.begin(), entry.second.end());
            }

            edges_loaded += local_edges_loaded;
        }
    }

    // Clean up memory mapping
    munmap(mapped, map_size);
    close(fd);

    auto edge_end_time = std::chrono::high_resolution_clock::now();
    auto end_time = edge_end_time;
    // std::cout << "Loaded " << edges_loaded << " edges with memory mapping in " << std::chrono::duration_cast<std::chrono::milliseconds>(edge_end_time - edge_start_time).count() << " ms" << std::endl;
    // std::cout << "Total parsing time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " ms" << std::endl;

    return {nodes, adj, reverse_adj};
}
int main(int argc, char *argv[])
{
    std::string NETLIST_FILE, OUTPUT_FILE, DEVICE_FILE;

    if (argc == 4)
    {
        NETLIST_FILE = argv[1];
        OUTPUT_FILE = argv[2];
        DEVICE_FILE = argv[3]; // Default device file
    }
    else
    {
        // Show usage instructions
        std::cout << "Usage: " << argv[0] << " <netlist_file> <output_file> [device_file]" << std::endl;
        return 1;
    }
    omp_set_num_threads(8); // Limit to 8 threads
    std::cout << "OpenMP thread count limited to 8 threads" << std::endl;

    // Parse netlist file first
    auto nets = parse_netlist_file(NETLIST_FILE);
    std::cout << "Loaded " << nets.size() << " nets from netlist" << std::endl;

    // Selectively parse device file using netlist information
    auto device_data = parse_device_file(DEVICE_FILE, nets);
    NodeMap nodes = device_data.nodes;
    AdjacencyList adj = device_data.adj;
    AdjacencyList reverse_adj = device_data.reverse_adj;

    // Generate routing with the improved function
    auto routing_paths = generate_routing(nets, adj, reverse_adj, nodes);

    // Write results using the new format
    write_output(nets, routing_paths, OUTPUT_FILE);
    auto program_end_time = std::chrono::high_resolution_clock::now();
    std::cout << "Total program execution time: " << std::chrono::duration_cast<std::chrono::milliseconds>(program_end_time - program_start_time).count() << " ms" << std::endl;
    return 0;
}
