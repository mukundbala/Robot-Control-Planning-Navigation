#ifndef TBOT__DJIKSTRA_H
#define TBOT__DJIKSTRA_H
#include "map_data.h"
#include <deque>
class Djikstra
{
public:

    struct Node
    {
        double g;
        bool visited;
        bot_utils::Index idx;
        bot_utils::Index parent;
        Node();
    };
    struct GOpen
    {
        double g;
        bot_utils::Index idx;
        GOpen();
        GOpen(double g, bot_utils::Index idx);
    };

    Djikstra(MapData &map);

    bot_utils::Index start;

    bot_utils::Index goal;

    MapData &map_;

    std::vector<bot_utils::Index> plan(bot_utils::Index idx_start, bot_utils::Index idx_goal , MapData& map_);

    // std::vector<bot_utils::Pos2D> plan(bot_utils::Pos2D pos_start, bot_utils::Pos2D pos_goal , MapData& map_);

    bot_utils::Index e_plan(bot_utils::Index idx_curr, bot_utils::Index idx_bad_goal , MapData& map_);
    
    std::vector<bot_utils::Pos2D> path_smoothing(std::vector<bot_utils::Pos2D>& path);
    
    //Insertion sort algorithm to sort nodes in GOpen list based on ascending g-cost
    void insertionSort(std::deque<GOpen> &list, Node *node);
    //Insertion sort algorithm to sort neighbouring cells based on ascending g-cost
    void insertionSort(std::vector<Node> &list, Node *node);

private:

    std::vector<Node> nodes; // keeps a record of the cheapest cost of every cell in the grid, as well as their parents
    
    std::deque<GOpen> open_list; // List to sort the unvisited cells based on their g-cost

    bot_utils::Index NB_LUT[8] = {{1,0}, {1,1}, {0,1}, {-1,1}, {-1,0}, {-1,-1}, {0,-1}, {1,-1}};

    void add_to_open(Node * node);

    Node * poll_from_open();


    bot_utils::Index pos2idx(bot_utils::Pos2D &pos);

    bot_utils::Pos2D idx2pos(bot_utils::Index &idx);

    bool oob(bot_utils::Index &idx);

    bool oob(bot_utils::Pos2D &pos);

    int flatten(bot_utils::Index &idx);

    bool checkCell(bot_utils::Index &idx);

    bool checkCell(bot_utils::Pos2D &pos);
};
#endif //TBOT__DJIKSTRA_H