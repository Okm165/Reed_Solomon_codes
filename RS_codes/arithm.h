#ifndef ARITHM_H
#define ARITHM_H

#include <vector>
#include <stdint.h>
#include <iostream>

#include "utils.h"


struct SymbF
{
    uint32_t size;                      // field size
    uint32_t primpoly;                  // primitive polynomial (generator)
    Table exp_table;                    // exponent table, [exponent] returns numeric value
    Table num_table;                    // numeric table, [numeric] returns exponent value

    SymbF(){};
    SymbF(uint32_t size, uint32_t primpoly)
    {
        this->size = size;
        this->primpoly = primpoly;
        // generate exp_table
        this->exp_table.resize(size-1);
        uint32_t tmp = 1;
        exp_table[0] = tmp;
        for(uint32_t i = 1; i < this->size-1; i++)
        {
            tmp = tmp << 1;
            if(tmp & 256)
                tmp ^= this->primpoly;
            this->exp_table[i] = tmp;
        }
        // generate num_table
        this->num_table.resize(size-1);
        this->num_table[0] = 0;
        for(uint32_t i = 0; i < this->size-1; i++)
            this->num_table[this->exp_table[i]] = i;
    }
    ~SymbF(){}

    bool operator==(const SymbF& val)
    {
        return (this->size == val.size)&&(this->primpoly == val.primpoly);
    }
    void operator=(const SymbF& val)
    {
        this->size = val.size;
        this->primpoly = val.primpoly;
        this->exp_table = val.exp_table;
        this->num_table = val.num_table;
    }

    friend std::ostream& operator<< (std::ostream& stream, const SymbF& symbf)
    {
        stream << "symbf: ";
        stream << "size: " << symbf.size << " " << "primpoly: " << symbf.primpoly << "\n";
        stream << "exp_table:\n" << print_table(symbf.exp_table) << "\n";
        stream << "num_table:\n" << print_table(symbf.num_table);
        return stream;
    }
};

enum SymbF_elem_mode{EXP,NUM};
struct SymbF_elem
{
    SymbF* symbf;           // symbol field pointer
    uint32_t exp;           // exponent value
    uint32_t num;           // numeric value

    SymbF_elem(){}
    SymbF_elem(SymbF* symbf, uint32_t val, SymbF_elem_mode mode)
    {
        this->symbf = symbf;
        if(mode == EXP)
        {
            this->exp = val;
            this->num = this->symbf->exp_table[val];
        }
        else if(mode == NUM)
        {
            this->num = val;
            this->exp = this->symbf->num_table[val];
        }
    }
    ~SymbF_elem(){}

    uint32_t norm(int32_t val)
    {
        while(val < 0)
            val += this->symbf->size;
        return val%(this->symbf->size);
    }

    bool operator==(SymbF_elem& val)
    {
        return (this->symbf == val.symbf)&&(this->exp == val.exp)&&(this->num == val.num);
    }
    void operator=(SymbF_elem& val)
    {
        this->symbf = val.symbf;
        this->exp = val.exp;
        this->num = val.num;
    }

    SymbF_elem operator^ (const SymbF_elem& val)
    {
        SymbF_elem elem;
        elem.symbf = this->symbf;
        elem.num = (this->num^val.num);
        elem.exp = this->symbf->num_table[elem.num];
        return elem;
    }
    SymbF_elem operator+ (const SymbF_elem& val)
    {
        SymbF_elem elem;
        elem.symbf = this->symbf;
        elem.num = elem.norm(this->num + val.num);
        elem.exp = this->symbf->num_table[elem.num];
        return elem;
    }
    SymbF_elem operator- (const SymbF_elem& val)
    {
        SymbF_elem elem;
        elem.symbf = this->symbf;
        elem.num = elem.norm(this->num - val.num);
        elem.exp = this->symbf->num_table[elem.num];
        return elem;
    }
    SymbF_elem operator* (const SymbF_elem& val)
    {
        SymbF_elem elem;
        elem.symbf = this->symbf;
        elem.exp = elem.norm(this->exp + val.exp);
        elem.num = elem.symbf->exp_table[elem.exp];
        return elem;
    }
    SymbF_elem operator/ (const SymbF_elem& val)
    {   
        return *this * ~(SymbF_elem&)val;
    }
    SymbF_elem operator~ ()
    {
        return SymbF_elem(this->symbf, (this->symbf->size-this->exp), EXP);
    }

    void operator^= (const SymbF_elem& val)
    {
        this->num = (this->num^val.num);
        this->exp = this->symbf->num_table[this->num];
    }
    void operator+= (const SymbF_elem& val)
    {
        this->num = this->norm(this->num + val.num);
        this->exp = this->symbf->num_table[this->num];
    }
    void operator-= (const SymbF_elem& val)
    {
        this->num = this->norm(this->num - val.num);
        this->exp = this->symbf->num_table[this->num];
    }
    void operator*= (const SymbF_elem& val)
    {
        this->exp = this->norm(this->exp + val.exp);
        this->num = this->symbf->exp_table[this->exp];
    }
    void operator/= (const SymbF_elem& val)
    {   
        *this *= ~(SymbF_elem&)val;
    }

    friend std::ostream& operator<< (std::ostream& stream, const SymbF_elem& symbf_elem)
    {
        stream << "symbf_elem: ";
        stream << "[" << symbf_elem.exp << " " << symbf_elem.num << "]";
        return stream;
    }
};

struct PolyF
{
    uint32_t size;          // field size
    PolyF(){}
    PolyF(uint32_t size)
    {
        this->size = size;
    }
    ~PolyF(){}

    bool operator==(const PolyF& val)
    {
        return (this->size == val.size);
    }
    void operator=(const PolyF& val)
    {
        this->size = val.size;
    }

    friend std::ostream& operator<< (std::ostream& stream, const PolyF& polyf)
    {
        stream << "polyf: ";
        stream << "size: " << polyf.size;
        return stream;
    }
};

struct PolyF_elem
{
    PolyF* polyf;           // poly field pointer
    Polynomial poly;        // polynomial obj

    PolyF_elem(){};
    PolyF_elem(PolyF* polyf, const Polynomial& poly)
    {
        this->polyf = polyf;
        this->poly = poly;
        this->norm();
    }
    ~PolyF_elem(){}

    uint32_t norm(int32_t val)
    {
        while(val < 0)
            val += this->polyf->size;
        return val%(this->polyf->size);
    }
    void norm()
    {
        for(uint32_t it = 0; it < this->poly.size(); it++)
            this->poly[it] = this->norm(this->poly[it]);
    }
    void red()
    {
        uint32_t index = this->poly.size()-1;
        while(this->poly[index] == 0 && index > 0)
            index--;
        if(this->poly.size()-1 != index)
            this->poly.resize(index+1);
    }

    bool operator==(const PolyF_elem& val)
    {
        return (this->polyf == val.polyf)&&(this->poly == val.poly);
    }
    bool operator==(uint32_t val)
    {
        return (this->poly.size() == 1)&&(this->poly[0] == val);
    }
    void operator=(const PolyF_elem& val)
    {
        this->polyf = val.polyf;
        this->poly = val.poly;
    }

    PolyF_elem operator<< (const uint32_t& val)
    {
        this->poly.insert(this->poly.end(), val, 0);
        return *this;
    }
    PolyF_elem operator>> (const uint32_t& val)
    {
        this->poly.insert(this->poly.begin(), val, 0);
        return *this;
    }
    PolyF_elem operator+ (const PolyF_elem& val)
    {
        PolyF_elem ret;
        ret.polyf = this->polyf;

        uint32_t len_a = this->poly.size();
        uint32_t len_b = val.poly.size();

        // select longer one
        if(len_a >= len_b)
        {
            ret.poly = this->poly;
            for(uint32_t it = 0; it < len_b; it++)
                ret.poly[it] = ret.norm(ret.poly[it] + val.poly[it]);
        }
        else
        {
            ret.poly = val.poly;
            for(uint32_t it = 0; it < len_a; it++)
                ret.poly[it] = ret.norm(ret.poly[it] + this->poly[it]);
        }
        return ret;
    }
    PolyF_elem operator+ (const int32_t& val)
    {
        PolyF_elem ret = *this;
        ret.poly[0] = ret.norm(ret.poly[0] + val);
        return ret;
    }
    PolyF_elem operator- (const PolyF_elem& val)
    {
        PolyF_elem ret;
        ret.polyf = this->polyf;

        uint32_t len_a = this->poly.size();
        uint32_t len_b = val.poly.size();
        
        // select longer one
        if(len_a >= len_b)
        {
            ret.poly = this->poly;
            for(uint32_t it = 0; it < len_b; it++)
                ret.poly[it] = ret.norm(ret.poly[it] - val.poly[it]);
        }
        else
        {
            ret.poly = val.poly;
            for(uint32_t it = 0; it < len_a; it++)
                ret.poly[it] = ret.norm(ret.poly[it] - this->poly[it]);
        }
        return ret;
    }
    PolyF_elem operator- (const int32_t& val)
    {
        PolyF_elem ret = *this;
        ret.poly[0] = ret.norm(ret.poly[0] - val);
        return ret;
    }
    PolyF_elem operator* (const PolyF_elem& val)
    {
        PolyF_elem ret;
        ret.polyf = this->polyf;

        uint32_t len_a = this->poly.size();
        uint32_t len_b = val.poly.size();
        ret.poly = Polynomial(len_a+len_b-1, 0);

        for(uint32_t i = 0; i < len_a; i++)
            for(uint32_t j = 0; j < len_b; j++)
                ret.poly[i+j] = ret.norm(ret.poly[i+j] + this->poly[i] * val.poly[j]);
        return ret;
    }
    PolyF_elem operator* (const int32_t& val)
    {
        PolyF_elem ret = *this;
        for(uint32_t i = 0; i < this->poly.size(); i++)
            ret.poly[i] = ret.norm(ret.poly[i] * val);
        return ret;
    }
    PolyF_elem operator% (const PolyF_elem& div)
    {
        // PolyF_elem quo;
        // quo.polyf = this->polyf;
        PolyF_elem rem = *this;
        while(rem.poly.size() >= div.poly.size())
        {
            int32_t r = rem.poly[rem.poly.size()-1] / div.poly[div.poly.size()-1];
            // quo.poly.insert(quo.poly.begin(), rval);
            rem -= (((PolyF_elem&)div * r) >> (rem.poly.size() - div.poly.size()));
            if(rem.poly[rem.poly.size()-1] != 0){break;}
            rem.poly.resize(rem.poly.size()-1);
        }
        return rem;
    }
    PolyF_elem operator/ (const PolyF_elem& div)
    {
        PolyF_elem quo;
        quo.polyf = this->polyf;
        PolyF_elem rem = *this;
        while(rem.poly.size() >= div.poly.size())
        {
            int32_t r = rem.poly[rem.poly.size()-1] / div.poly[div.poly.size()-1];
            quo.poly.insert(quo.poly.begin(), r);
            rem -= (((PolyF_elem&)div * r) >> (rem.poly.size() - div.poly.size()));
            if(rem.poly[rem.poly.size()-1] != 0){break;}
            rem.poly.resize(rem.poly.size()-1);
        }
        return quo;
    }

    void operator+= (const PolyF_elem& val)
    {
        uint32_t len_a = this->poly.size();
        uint32_t len_b = val.poly.size();
        Polynomial poly;
        // select longer one
        if(len_a >= len_b)
        {
            poly = this->poly;
            for(uint32_t it = 0; it < len_b; it++)
                poly[it] = this->norm(poly[it] + val.poly[it]);
        }
        else
        {
            poly = val.poly;
            for(uint32_t it = 0; it < len_a; it++)
                poly[it] = this->norm(poly[it] + this->poly[it]);
        }
        
        this->poly = poly;
    }
    void operator+= (const int32_t& val)
    {
        this->poly[0] = this->norm(this->poly[0] + val);
    }
    void operator-= (const PolyF_elem& val)
    {
        uint32_t len_a = this->poly.size();
        uint32_t len_b = val.poly.size();
        Polynomial poly;
        // select longer one
        if(len_a >= len_b)
        {
            poly = this->poly;
            for(uint32_t it = 0; it < len_b; it++)
                poly[it] = this->norm(poly[it] - val.poly[it]);
        }
        else
        {
            poly = val.poly;
            for(uint32_t it = 0; it < len_a; it++)
                poly[it] = this->norm(poly[it] - this->poly[it]);
        }
        this->poly = poly;
    }
    void operator-= (const int32_t& val)
    {
        this->poly[0] = this->norm(this->poly[0] - val);
    }
    void operator*= (const PolyF_elem& val)
    {
        uint32_t len_a = this->poly.size();
        uint32_t len_b = val.poly.size();
        Polynomial poly(len_a+len_b-1, 0);

        for(uint32_t i = 0; i < len_a; i++)
            for(uint32_t j = 0; j < len_b; j++)
                poly[i+j] = this->norm(poly[i+j] + this->poly[i] * val.poly[j]);
        this->poly = poly;
    }
    void operator*= (const int32_t& val)
    {
        uint32_t len_a = this->poly.size();
        for(uint32_t i = 0; i < len_a; i++)
            this->poly[i] = this->norm(this->poly[i] * val);
    }
    void operator%= (const PolyF_elem& div)
    {
        PolyF_elem rem = *this;
        while(rem.poly.size() >= div.poly.size())
        {
            int32_t rval = rem.poly[rem.poly.size()-1] / div.poly[div.poly.size()-1];
            PolyF_elem r (this->polyf, Polynomial{rval});
            // quo.poly.insert(quo.poly.begin(), rval);
            rem -= (((PolyF_elem&)div * r) >> (rem.poly.size() - div.poly.size()));
            if(rem.poly[rem.poly.size()-1] != 0){break;}
            rem.poly.resize(rem.poly.size()-1);
        }
        *this = rem;
    }
    void operator/= (const PolyF_elem& div)
    {
        PolyF_elem quo;
        quo.polyf = this->polyf;
        PolyF_elem rem = *this;
        while(rem.poly.size() >= div.poly.size())
        {
            int32_t r = rem.poly[rem.poly.size()-1] / div.poly[div.poly.size()-1];
            quo.poly.insert(quo.poly.begin(), r);
            rem -= (((PolyF_elem&)div * r) >> (rem.poly.size() - div.poly.size()));
            if(rem.poly[rem.poly.size()-1] != 0){break;}
            rem.poly.resize(rem.poly.size()-1);
        }
        *this = quo;
    }

    friend std::ostream& operator<< (std::ostream& stream, const PolyF_elem& val)
    {
        stream << print_poly((Polynomial*)&(val.poly));
        return stream;
    }
};

#endif