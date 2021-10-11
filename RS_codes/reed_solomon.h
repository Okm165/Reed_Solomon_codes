#ifndef REED_SOLOMON_H
#define REED_SOLOMON_H

#include "utils.h"
#include "arithm.h"

struct ReedSolomonCode
{
    uint32_t n;             // output codeword length in symbols
    uint32_t k;             // input codeword length in symbols
    uint32_t s;             // symbol length (number of bits per symbol)

    uint32_t t;             // error correcting copacity (number of correctable errors)

    SymbF symbfield;
    PolyF polyfield;

    ReedSolomonCode(){}
    ReedSolomonCode(uint32_t n, uint32_t k, uint32_t symbfield_primpoly, uint32_t symb_size = 8)
    {
        this->n = n;
        this->k = k;
        this->s = symb_size;
        this->t = (n-k)/2;

        this->symbfield = SymbF (pow(2,s), symbfield_primpoly);
        this->polyfield = PolyF (pow(2,s));
    }
    ~ReedSolomonCode(){}

    PolyF_elem composeGenPoly(uint32_t len, uint32_t eps = 0)
    {
        PolyF_elem ret(&polyfield, Polynomial{-symbfield.exp_table[eps], 1});
        for(uint32_t it = 1; it < len; it++)
        {
            ret *= PolyF_elem(&polyfield, Polynomial{-symbfield.exp_table[eps+it], 1});
        }
        return ret;
    }
    
    
    // encodePoly
    
    // decodeCodeWord
    // evalSyndromes

};

#endif