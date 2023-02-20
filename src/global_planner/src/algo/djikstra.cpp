#include "djikstra.h"

Djikstra::Node::Node()
        : g(INFINITY), visited(false), idx(-1, -1), parent(-1, -1) 
        {}
Djikstra::GOpen::GOpen()
        : g(0), idx(-1, -1) 
        {}
Djikstra::GOpen::GOpen(double g, bot_utils::Index idx)
        : g(g), idx(idx) 
        {}
Djikstra::Djikstra(MapData &map)
        : start(-1, -1), goal(-1, -1), map_(map), nodes(map.total_cells_), open_list() 
        {
            int k = 0;
            for (int i = 0; i < map.map_size_.i; ++i)
            {
                for (int j = 0; j < map.map_size_.j; ++j)
                {
                    nodes[k].idx.i = i;
                    nodes[k].idx.j = j;
                    ++k;
                }
            }
        }

void Djikstra::add_to_open(Node * node)
{
    //  Checking if its the first node
    if (open_list.empty())
    {
        open_list.emplace_front(node->g, node->idx);
    }
    insertionSort(open_list, node);
}

void Djikstra::insertionSort(std::deque<GOpen> &list, Node *node)
{
    //  List is already populated and find the correct position of the node using insertion sort
    for (int n = 0; n < list.size(); ++n)
    {
        if (node->g < list[n].g)
        {
            list.emplace(list.begin() + n, node->g, node->idx);
            break;
        }
    }
}

void Djikstra::insertionSort(std::vector<Node> &list, Node * node)
{
    //  List is already populated and find the correct position of the node using insertion sort
    for (int n = 0; n < list.size(); ++n)
    {
        if (node->g < list[n].g)
        {
            list.emplace(list.begin() + n, node);
            break;
        }
    }
}

Djikstra::Node * Djikstra::poll_from_open()
{
    // Retrive neccessary information from cell before removing it
    Node * node;
    GOpen &curr_cell = open_list.front();
    int k = flatten(curr_cell.idx);
    node = &(nodes[k]);

    open_list.pop_front();

    return node;
}

std::vector<bot_utils::Index> Djikstra::plan(bot_utils::Index idx_start, bot_utils::Index idx_goal , MapData& map_)
{
    std::vector<bot_utils::Index> path_idx;

    for (Djikstra::Node & node : nodes)
    {
        node.g = 1e5;
        node.visited = false;
    }

    // Start node, set g-cost to zero
    int k = Djikstra::flatten(idx_start);
    Djikstra::Node * node = &(nodes[k]);
    node->g = 0;

    // Add start node to GOpen list
    add_to_open(node);

    while (!open_list.empty())
    {
        node = poll_from_open();

        // Check if node has been visited
        if (node->visited)
        {
            continue;
        }
        node->visited = true;

        // Trace back and return the path if at the goal
        if (node->idx.i == idx_goal.i && node->idx.j == idx_goal.j)
        {
            path_idx.push_back(node->idx);

            while (node->idx.i != idx_start.i || node->idx.j != idx_start.j)
            {
                int k = flatten(node->parent);
                node = &(nodes[k]);
                path_idx.push_back(node->idx);
            }

            break;

        }
        // Check surrounding neighbours and update g-costs if cheaper
        bool is_cardinal = true;
        for (int dir = 0; dir < 8; ++dir)
        {
            // Get coordinates of the different neighbours
            bot_utils::Index & idx_nb_relative = NB_LUT[dir];
            bot_utils::Index idx_nb(
                node->idx.i + idx_nb_relative.i,
                node->idx.j + idx_nb_relative.j
            );

            // Ignore cell if it's not available
            if (!checkCell(idx_nb))
            {
                // Move to the next neighbouring cell
                continue;
            }
            // Update g-cost of node
            double g_cost = node->g;
            if (is_cardinal) 
                g_cost += 1;
            else
                g_cost += M_SQRT2;

            // compare the cost to any previous costs. If cheaper, mark the node as the parent
            int nb_k = flatten(idx_nb);
            Node & nb_node = nodes[nb_k]; // use reference so changing nb_node changes nodes[k]
            if (nb_node.g > g_cost + 1e-5)
            {   
                nb_node.g = g_cost;
                nb_node.parent = node->idx;

                // add to open
                add_to_open(&nb_node); // & a reference means getting the pointer (address) to the reference's object.
            }

            // toggle is_cardinal
            is_cardinal = !is_cardinal;
        }
        
    }
    open_list.clear();
    return path_idx;
}

bot_utils::Index Djikstra::e_plan(bot_utils::Index idx_curr, bot_utils::Index idx_bad_goal , MapData& map_)
{
    bot_utils::Index new_goal_idx;
    std::vector<Node> neighbour_idx;

    for (Djikstra::Node & node : nodes)
    {
        node.g = 1e5;
        node.visited = false;
    }

    // Start node is the bad goal, set g-cost to zero
    int k = Djikstra::flatten(idx_bad_goal);
    Djikstra::Node * node = &(nodes[k]);
    node->g = 0;

    // Add start node to GOpen list
    add_to_open(node);

    while (!open_list.empty())
    {
        node = poll_from_open();

        // Check if node has been visited
        if (node->visited)
        {
            continue;
        }
        node->visited = true;

        // Check surrounding neighbours for free cell
        for (int dir = 0; dir < 8; ++dir)
        {
            // Get coordinates of the different neighbours
            bot_utils::Index & idx_nb_relative = NB_LUT[dir];
            bot_utils::Index idx_nb(
                node->idx.i + idx_nb_relative.i,
                node->idx.j + idx_nb_relative.j
            );

            int k_nb = flatten(idx_nb);
            Node * nb_node = &(nodes[k_nb]);
            // if cell is not available, add the cell to the list and check its neighbours for free cells
            if (!checkCell(idx_nb))
            {
                add_to_open(nb_node);
                continue;
            }
            // Neighbouring cell is free, add it to the vector<Node>
            insertionSort(neighbour_idx, nb_node);
        }
        if (neighbour_idx.empty())
        {
            break;
        }
    }
    node = &(neighbour_idx[0]);
    new_goal_idx = node->idx;
    return new_goal_idx;
}   

bot_utils::Index Djikstra::pos2idx(bot_utils::Pos2D &pos)
{
    int i = round((pos.x - map_.origin_.x) / map_.cell_size_);
    int j = round((pos.y - map_.origin_.y) / map_.cell_size_);
    return bot_utils::Index(i,j);
}

bot_utils::Pos2D Djikstra::idx2pos(bot_utils::Index &idx)
{
    double x = idx.i * map_.cell_size_ + map_.origin_.x;
    double y = idx.j * map_.cell_size_ + map_.origin_.y;
    return bot_utils::Pos2D(x,y); 
}

bool Djikstra::oob(bot_utils::Index &idx)
{
    return (idx.i <=0 || idx.i >= map_.map_size_.i || idx.j < 0 && idx.j >= map_.map_size_.j);
}

bool Djikstra::oob(bot_utils::Pos2D &pos)
{
    bot_utils::Index idx = pos2idx(pos);
    return oob(idx);
}

int Djikstra::flatten(bot_utils::Index &idx)
{
    return (idx.i * map_.map_size_.j) + idx.j;
}

bool Djikstra::checkCell(bot_utils::Index &idx)
{
    if (oob(idx))
    {
        return false;
    }

    int k = flatten(idx);

    if (map_.grid_inflation_.at(k) > 0)
    {
        return false; //in map, inflated
    }

    else if (map_.grid_logodds_.at(k) > map_.lo_thresh_)
    {
        return false; //not inflated, lo occuped
    }

    else
    {
        return true;
    }
}

bool Djikstra::checkCell(bot_utils::Pos2D &pos)
{
    bot_utils::Index idx = pos2idx(pos);
    return checkCell(idx);
}