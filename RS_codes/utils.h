#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <vector>
#include <iostream>
#include <stdint.h>

typedef std::vector<int32_t> Polynomial;
typedef std::vector<int32_t> Table;
std::string print_table(const Table& data);
std::string print_poly(Polynomial* poly);
std::ostream& operator<< (std::ostream& stream, const Table& table);

#endif