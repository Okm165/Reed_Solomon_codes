import numpy
import copy

input = [0,0,1,1,1,1,0,1,1,1,0,1,1,1,0,0,1,0]

def berlekamp_massey_algorithm(block_data):
    n = len(block_data)
    c = numpy.zeros(n)
    b = numpy.zeros(n)
    c[0], b[0] = 1, 1
    l, m, i = 0, -1, 0
    int_data = [int(el) for el in block_data]
    while i < n:
        v = int_data[(i - l):i]
        v = v[::-1]
        cc = c[1:l + 1]
        d = (int_data[i] + numpy.dot(v, cc)) % 2
        if d == 1:
            temp = copy.copy(c)
            p = numpy.zeros(n)
            for j in range(0, l):
                if b[j] == 1:
                    p[j + i - m] = 1
            c = (c + p) % 2
            if l <= 0.5 * i:
                l = i + 1 - l
                m = i
                b = temp
        i += 1
    return c[:l+1]

seq = (berlekamp_massey_algorithm(input))
print(seq)
seqlen = len(seq)
for i in range(len(input)-seqlen+1):
    v = input[i:i+seqlen]
    v = v[::-1]
    d = numpy.dot(v, seq) % 2
    print(d)