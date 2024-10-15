#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include "GF.hpp"

using namespace std;

uint8_t mask[9] = { 0x00,0x00,0x00,0x00,0x0f,0x1f,0x3f,0x7f,0xff };

template<int N>
class gensbox {
public:
    vector<vector<uint8_t>> Affinematrix1={{0}};
    vector<vector<uint8_t>> Affinematrix2={{0}};
    vector<uint8_t> Affinevector1 = { 0 };
    vector<uint8_t> Affinevector2 = { 0 };
    vector<uint8_t> newSbox;
    galois<N> baseSbox;

    gensbox() 
    {
        genAffinepar();
        
        vector<uint8_t> space(1 << N);
        newSbox = space;
        int count = 0;
        for (uint8_t i = 0; count < (1 << N); i++,count++)
        {
            newSbox[i] = affineunit(baseSbox.Sbox[affineunit(i, Affinematrix1, Affinevector1)], Affinematrix2, Affinevector2);
        }
    }
    gensbox(vector<vector<uint8_t>> &a1, vector<vector<uint8_t>> &a2, vector<uint8_t> &b1, vector<uint8_t> &b2)
    {

        Affinematrix1 = a1;
        Affinematrix2 = a2;
        Affinevector1 = b1;
        Affinevector2 = b2;
        vector<uint8_t> space(1 << N);
        newSbox = space;
        int count = 0;
        for (uint8_t i = 0; count < (1 << N); i++, count++)
        {
            newSbox[i] = affineunit(baseSbox.Sbox[affineunit(i, Affinematrix1, Affinevector1)], Affinematrix2, Affinevector2);
        }
    }
    uint8_t random_bit() 
    {
        return rand() % 2;
    }

    vector<vector<uint8_t>> generate_matrix(int n) 
    {
        vector<vector<uint8_t>> matrix(n, vector<uint8_t>(n));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                matrix[i][j] = random_bit();
            }
        }
        return matrix;
    }
    vector<uint8_t> generate_vector(int n)
    {
        vector<uint8_t> matrix(n);
        for (int j = 0; j < n; j++) {
            matrix[j] = random_bit();
        }
        return matrix;
    }


    uint8_t determinant(const vector<vector<uint8_t>>& matrix,uint8_t n) 
    {
        if (n == 1) return matrix[0][0];

        uint8_t det = 0;
        vector<vector<uint8_t>> submatrix(n, vector<uint8_t>(n));
        for (uint8_t x = 0; x < n; x++) {
            uint8_t subi = 0;
            for (uint8_t i = 1; i < n; i++) {
                uint8_t subj = 0;
                for (uint8_t j = 0; j < n; j++) {
                    if (j == x)
                        continue;
                    submatrix[subi][subj] = matrix[i][j];
                    subj++;
                }
                subi++;
            }
            det += (x % 2 == 0 ? 1 : -1) * matrix[0][x] * determinant(submatrix, n - 1);
        }

        return det%2;
    }


    bool is_invertible(const vector<vector<uint8_t>>& matrix) 
    {
        return determinant(matrix,N);
    }


    void genAffinepar()
    {
        do {
            Affinematrix1 = generate_matrix(N);
        } while (!is_invertible(Affinematrix1));


        do {
            Affinematrix2 = generate_matrix(N);
        } while (!is_invertible(Affinematrix2));
        Affinevector1 = generate_vector(N);
        Affinevector2 = generate_vector(N);
    }


    uint8_t affineunit(uint8_t a,vector<vector<uint8_t>> matrix,vector<uint8_t> vector)
    {
        uint8_t result=0;
        a = a & mask[N];
        for (uint8_t i = 0; i < N; i++)
        {
            uint8_t temp=0;
            for (uint8_t j = 0; j < N; j++)
            {
                temp ^= ((a >> j) & 1) & matrix[i][j];
            }
            result ^= temp << i;
        }
        for (uint8_t k = 0; k < N; k++)
        {
            result ^= vector[k] << k;
        }
        return result;
    }

    void pruint8_t_matrix_vector()
    {
        for (uint8_t i = 0; i < N;i++) {
            for (uint8_t j = 0; j < N;j++) {
                cout <<(int) Affinematrix1[i][j] << ",";
            }
            cout << endl;
        }

        for (uint8_t i = 0; i < N; i++) {
            for (uint8_t j = 0; j < N; j++) {
                cout << (int)Affinematrix2[i][j] << ",";
            }
            cout << endl;
        }

        for (uint8_t j = 0; j < N; j++) {
            cout << (int)Affinevector1[j] << " ";
        }
        cout << endl;

        for (uint8_t j = 0; j < N; j++) {
            cout << (int)Affinevector2[j] << " ";
        }
        cout << endl;
    }
    string evaluate_LUT()
    {
        std::stringstream ss;

        ss <<std::hex<< std::setfill('0');
        for (int i = 0; i < (1 << N); i++)
        {
            ss<<setw(2)<<(int)newSbox[i];
        }
        return ss.str();
    }
    void print_sbox()
    {

        if (N == 8)
        {
            for (int i = 0; i < 16; i++)
            {
                for (int j = 0; j < 16; j++)
                {
                    std::cout << std::hex << (int)newSbox[i * 16 + j] << " ";
                }
                std::cout << std::endl;
            }
        }
        else
        {
            for (int i = 0; i < (1 << N); i++)
            {
                std::cout << std::setfill('0') << std::setw(2) << std::hex << (int)newSbox[i];
            }
            std::cout << std::endl;
        }
    }
};



uint8_t main() {
    srand(time(0));
    std::vector<vector<uint8_t>> a1 =
    {
        {1,0,0,0,0,0,0,0},
        {0,1,0,0,0,0,0,0},
        {0,0,1,0,0,0,0,0},
        {0,0,0,1,0,0,0,0},
        {0,0,0,0,1,0,0,0},
        {0,0,0,0,0,1,0,0},
        {0,0,0,0,0,0,1,0},
        {0,0,0,0,0,0,0,1}
    };
    std::vector<uint8_t> b1= {0,0,0,0,0,0,0,0};
    std::vector<vector<uint8_t>> a2 =
    {
        {1,0,0,0,1,1,1,1},
        {1,1,0,0,0,1,1,1},
        {1,1,1,0,0,0,1,1},
        {1,1,1,1,0,0,0,1},
        {1,1,1,1,1,0,0,0},
        {0,1,1,1,1,1,0,0},
        {0,0,1,1,1,1,1,0},
        {0,0,0,1,1,1,1,1}
    };
    std::vector<uint8_t> b2 = { 1,1,0,0,0,1,1,0 };
    gensbox<8> a(a1, a2, b1, b2);
    //a.baseSbox.print_sbox();
    //a.print_sbox();
    string ss = a.evaluate_LUT();
    cout << ss << endl;
    return 0;
}
