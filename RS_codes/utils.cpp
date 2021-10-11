#include "utils.h"

std::string print_table(const Table& data)
{
    std::string ret;
    // get dimensions
    uint32_t dim_x = floor(sqrt(data.size()));
    uint32_t dim_y = ceil((data.size()+0.0)/dim_x);

    uint32_t max_length[dim_y];

    for(uint32_t it = 0; it < dim_y; it++)
        max_length[it] = 0;
    
    for(uint32_t x = 0; x < dim_x; x++)
    {
        for(uint32_t y = 0; y < dim_y; y++)
        {
            uint32_t curr_len = std::to_string(data[x*dim_y+y]).size();
            if(max_length[y] < curr_len)
                max_length[y] = curr_len;
        }
    }
    
    for(uint32_t x = 0; x < dim_x; x++)
    {
        for(uint32_t y = 0; y < dim_y; y++)
        {
            std::string elem = std::to_string(data[x*dim_y+y]);
            while(elem.size() < max_length[y])
            {
                elem += " ";
            }
            ret += elem;
            ret += " ";
        }
        ret += "\n";   
    }
    return ret;
}

std::string print_poly(Polynomial* poly)
{
    std::string ret;
    ret = "[ ";
    for(uint32_t it = 0; it < poly->size(); it++)
    {
        ret += std::to_string((*poly)[it]) + " ";
    }
    ret += "]";
    return ret;
}

std::ostream& operator<< (std::ostream& stream, const Table& table)
{
    stream << "Table: ";
    stream << "size: " << table.size() << "\n";
    stream << print_table(table);
    return stream;
}