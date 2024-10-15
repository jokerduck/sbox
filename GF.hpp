

#include <iostream>
#include<vector>
#include<iomanip>

//std::vector<std::vector<uint16_t>> irr = {
//    {0x00},{0x00},{0x00},{0x00},
//    {0x13,0x19,0xf},//F24
//    {0x25,0x29,0x2f,0x37,0x3b,0x3d},//F25
//    {0x43,0x49,0x57,0x5b,0x61,0x67,0x6d,0x73,0x75},//F26
//    {0x83,0x89,0x8f,0x91,0x9d,0xa7,0xab,0xb9,0xbf,0xc1,0xcb,0xd3,0xd5,0x35,0xef,0xf1,0xf7,0xfd},//F27
//    {0x11b,0x11d,0x12b,0x12d,0x139,0x13f,0x14d,0x15f,0x163,0x165,0x169,0x171,0x177,0x17b,0x187,0x18b,0x18d,0x19f,0x1a3,0x1a9,0x1b1,0x1bb,0x1c3,0x1cf,0x1d7,0x1dd,0x1e7,0x1f3,0x1f5,0x1f9} };
static const std::vector<uint16_t> irr = { 0x00,0x00,0x00,0x0,0x13,0x25,0x43,0x83,0x11b };
template<int N>
class galois {
public:
    uint8_t x;
    uint8_t y;
    std::vector<uint8_t> Sbox;
    uint16_t Sboxirr = irr[N];
    galois()
    {
        x = 0;
        y = 0;
        std::vector<uint8_t> temp(1 << N);
        Sbox = temp;
        SboxGF();
    }
    galois(uint16_t inirr)
    {
        Sboxirr = inirr;
        SboxGF();
    }
    uint8_t gf_add(uint8_t a, uint8_t b) {
        return a ^ b;
    }

    uint8_t gf_subtract(uint8_t a, uint8_t b) {
        return a ^ b;
    }

    uint8_t gf_multiply(uint8_t a, uint8_t b) {
        uint8_t p = 0;
        uint8_t counter;
        uint8_t hi_bit_set;
        for (counter = 0; counter < N; counter++) {
            if (b & uint8_t(1)) {
                p ^= a;
            }
            hi_bit_set = (uint8_t)(a & (1<<(N-1)));
            a <<= 1;
            if (hi_bit_set) {
                a ^= uint8_t(Sboxirr);
            }
            b >>= 1;
        }
        return p;
    }

    int get_highest_bit(uint16_t a)
    {
        int index = -1;
        while (a != 0) {
            a >>= 1;
            index++;
        }
        return index;
    }

    uint16_t gf_divide(uint16_t a, uint16_t b)
    {
        if (b == 0)
            throw std::invalid_argument("Division by zero!");

        uint16_t quotient = 0;
        uint16_t remainder = a;

        while (true)
        {
            uint16_t temp = b;
            uint16_t count = 1;

            int indexA = get_highest_bit(remainder);
            int indexB = get_highest_bit(b);
            if (indexA < indexB)
                return quotient;
            temp <<= indexA - indexB;
            count <<= indexA - indexB;

            quotient ^= count;
            remainder = gf_subtract(remainder, temp);
        }

        return quotient;
    }

    uint8_t gf_inverse(uint8_t a) {
        if (a == 0)
            return 0;

        uint8_t u = 1, v = 0;
        uint8_t lastu = 0, lastv = 1;
        uint16_t x = Sboxirr, y = a;
        uint8_t temp, quotient, remainder;

        quotient = gf_divide(x, y);
        remainder = x ^ gf_multiply(quotient, y) ^ Sboxirr;
        u = gf_subtract(u, gf_multiply(quotient, lastu));
        v = gf_subtract(v, gf_multiply(quotient, lastv));

        temp = u; u = lastu; lastu = temp;
        temp = v; v = lastv; lastv = temp;
        x = y;
        y = remainder;

        while (y != 0)
        {
            quotient = gf_divide(x, y);
            remainder = gf_subtract(x, gf_multiply(quotient, y));

            u = gf_subtract(u, gf_multiply(quotient, lastu));
            v = gf_subtract(v, gf_multiply(quotient, lastv));

            temp = u; u = lastu; lastu = temp;
            temp = v; v = lastv; lastv = temp;
            x = y;
            y = remainder;
        }

        return v;
    }
    void SboxGF()
    {
        int exe_count = 0;
        for (uint8_t x = 0; exe_count < (1<<N); x++, exe_count++)
        {
            Sbox[x] = gf_inverse(x);
        }
    }
    void print_sbox()
    {
        if (N == 8)
        {
            for (int i = 0; i < 16; i++)
            {
                for (int j = 0; j < 16; j++)
                {
                    std::cout <<std::hex<<(int) Sbox[i * 16 + j] << " ";
                }
                std::cout << std::endl;
            }
        }
        else
        {
            for (int i = 0; i < (1 << N); i++)
            {
                std::cout <<std::setfill('0')<<std::setw(2) << std::hex << (int)Sbox[i];
            }
            std::cout << std::endl;
        }
    }

};

