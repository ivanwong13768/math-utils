/*

--- most schoolbook way to do matrix calculation ---
* this only supports calculations between matrices with the same type, so use double as type where it is possible
* WARNING: this is extremely slow and computationally expensive for large matrices!

*/

// TODO: make multiple types compatible with one another (more operator overloading) (idk how to do cuz sizeof() is probably not possible)

#include <cmath>
#include <cassert>
#include "vector_operations.h"

#ifndef MATRIX
#define MATRIX
template <typename T>
class matrix2d {
    private:
        std::vector<std::vector<T>> m{};
        std::vector<size_t> s{0, 0};   // row, column
        /* helper functions */
        double det_helper(std::vector<std::vector<T>> m) {
            if ((m.size() == 1) && (m[0].size() == 1)) {
                return m[0][0];
            } else if ((m.size() == 2) && (m[0].size() == 2)) {
                return ((m[0][0] * m[1][1]) - (m[0][1] * m[1][0]));
            }
            T determinant = 0;
            for (size_t col = 0; col < m[0].size(); col++) {
                std::vector<std::vector<T>> cut_m{};
                for (size_t row = 1; row < m.size(); row++) {
                    std::vector<T> v = m[row];
                    v.erase(v.begin() + col);
                    cut_m.push_back(v);
                }
                determinant += m[0][col] * det_helper(cut_m) * pow(-1, col);
            }
            return determinant;
        }
        void pad() {
            size_t max_length = 0;
            for (size_t i = 0; i < m.size(); i++) {
                if (m[i].size() > max_length) {
                    max_length = m[i].size();
                }
            }
            for (size_t i = 0; i < m.size(); i++) {
                for (size_t j = m[i].size(); j < max_length; j++) {
                    m[i].push_back(static_cast<T>(0));
                }
            }
            s[0] = m.size();
            s[1] = max_length;
        }
    public:
        /* constructors */
        matrix2d(std::vector<std::vector<T>> m) {
            this->m = m;
            pad();
        }
        matrix2d(size_t row, size_t col, T default_value = static_cast<T>(0)) {
            s[0] = row;
            s[1] = col;
            std::vector<T> v{};
            for (size_t i = 0; i < col; i++) {
                v.push_back(default_value);
            }
            for (size_t i = 0; i < row; i++) {
                m.push_back(v);
            }
        }
        matrix2d(std::vector<size_t> size, T default_value) {
            assert(size.size() >= 2);
            if (size.size() > 2) {
                std::cerr << "Warning: length of size vector in initilization is not 2, using the first 2 elements only" << std::endl;
            }
            s[0] = size[0];
            s[1] = size[1];
            std::vector<T> v{};
            for (size_t i = 0; i < size[1]; i++) {
                v.push_back(default_value);
            }
            for (size_t i = 0; i < size[0]; i++) {
                m.push_back(v);
            }
        }
        matrix2d(size_t size[2], T default_value) {
            s[0] = size[0];
            s[1] = size[1];
            std::vector<T> v{};
            for (size_t i = 0; i < size[1]; i++) {
                v.push_back(default_value);
            }
            for (size_t i = 0; i < size[0]; i++) {
                m.push_back(v);
            }
        }
        matrix2d() {
            // default values
        }
        /* operators overloading */ 
        void operator= (matrix2d<T> m) {
            for (size_t i = 0; i < s[0]; i++) {
                this->m.pop_back();
            }
            for (size_t i = 0; i < m.shape()[0]; i++) {
                this->m.push_back(m[i]);
            }
            pad();
            // s[0] = m.shape()[0];
            // s[1] = m.shape()[1];
        }
        std::vector<T>& operator[] (size_t i) {
            return m[i];
        }
        /* matrix operations */
        std::vector<size_t> shape() {
            pad();
            return s;
        }
        matrix2d<T> transpose() {
            pad();
            matrix2d<T> new_m(s[1], s[0]);
            for (size_t col = 0; col < s[1]; col++) {
                for (size_t row = 0; row < s[0]; row++) {
                    new_m[row][col] = m[col][row];
                }
            }
            return new_m;
        }
        double det() {
            pad();
            assert(s[0] == s[1]);
            return det_helper(m);
        }
        matrix2d<T> inverse() {
            pad();
            assert(s[0] == s[1]);
            double d = det_helper(m);
            assert(d != 0);
            std::vector<std::vector<T>> coj{};
            for (size_t row = 0; row < s[0]; row++) {
                std::vector<T> det_v{};
                for (size_t col = 0; col < s[1]; col++) {
                    std::vector<std::vector<T>> temp_m = m;
                    temp_m.erase(temp_m.begin() + row);
                    for (size_t r = 0; r < temp_m.size(); r++) {
                        temp_m[r].erase(temp_m[r].begin() + col);
                    }
                    det_v.push_back(det_helper(temp_m) / d * pow(-1, ((row + col) % 2)));
                }
                coj.push_back(det_v);
            }
            matrix2d<T> inv(coj);
            return inv.transpose();
        }
        /* functions for convenience */
        template <typename A>
        matrix2d<A> to_type() {
            pad();
            matrix2d<A> new_m(s, 0);
            for (size_t row = 0; row < s[0]; row++) {
                for (size_t col = 0; col < s[1]; col++) {
                    new_m[row][col] = static_cast<T>(m[row][col]);
                }
            }
            return new_m;
        }
        void push_back_row(std::vector<T> v) {
            m.push_back(v);
            pad();
        }
        void pop_back_row() {
            m.pop_back();
            pad();
        }
        void push_back_col(std::vector<T> v) {
            pad();
            std::vector<T> v_padded = v;
            for (size_t row = 0; row < s[0]; row++) {
                if (row < v.size()) {
                    m[row].push_back(v_padded[row]);
                } else {
                    m[row].push_back(static_cast<T>(0));
                }
            }
            s[1] += 1;
        }
        void pop_back_col() {
            pad();
            for (size_t row = 0; row < s[0]; row++) {
                m[row].pop_back();
            }
            s[1] -= 1;
        }
};

template <typename T>
matrix2d<T> operator+ (matrix2d<T> m1, matrix2d<T> m2) {
    assert(m1.shape() == m2.shape());
    matrix2d<T> new_m = m1;
    for (size_t row = 0; row < m1.shape()[0]; row++) {
        for (size_t col = 0; col < m1[row].size(); col++) {
            new_m[row][col] += m2[row][col];
        }
    }
    return new_m;
}

template <typename T>
void operator+= (matrix2d<T>& m1, matrix2d<T> m2) {
    assert(m1.shape() == m2.shape());
    m1 = m1 + m2;
    return;
}

template <typename T>
matrix2d<T> operator- (matrix2d<T> m1, matrix2d<T> m2) {
    assert(m1.shape() == m2.shape());
    matrix2d<T> new_m = m1;
    for (size_t row = 0; row < m1.shape()[0]; row++) {
        for (size_t col = 0; col < m1[row].size(); col++) {
            new_m[row][col] -= m2[row][col];
        }
    }
    return new_m;
}

template <typename T>
void operator-= (matrix2d<T>& m1, matrix2d<T> m2) {
    assert(m1.shape() == m2.shape());
    m1 = m1 - m2;
    return;
}

template <typename T>
matrix2d<T> operator* (matrix2d<T> m1, matrix2d<T> m2) {
    assert(m1.shape()[1] != m2.shape()[0]);
    matrix2d<T> new_m(m1.shape()[0], m2.shape()[1]);
    matrix2d<T> m2_t = m2.transpose();
    for (size_t i = 0; i < m1.shape[0]; i++) {
        for (size_t j = 0; i < m2.shape[1]; j++) {
            new_m[i][j] = static_cast<T>(vector_sum(m1[i] * m2_t[j]));
        }
    }
    return new_m;
}

template <typename T>
matrix2d<T> operator* (matrix2d<T> m, T c) {
    matrix2d<T> new_m = m;
    for (size_t row = 0; row < m.shape()[0]; row++) {
        for (size_t col = 0; col < m[row].size(); col++) {
            new_m[row][col] *= c;
        }
    }
    return new_m;
}

template <typename T>
void operator*= (matrix2d<T>& m, T c) {
    m = m * c;
    return;
}

template <typename T>
matrix2d<T> operator/ (matrix2d<T> m, T c) {
    assert(c != 0);
    matrix2d<T> new_m = m;
    for (size_t row = 0; row < m.shape()[0]; row++) {
        for (size_t col = 0; col < m[row].size(); col++) {
            new_m[row][col] /= c;
        }
    }
    return new_m;
}

template <typename T>
void operator/= (matrix2d<T>& m, T c) {
    assert(c != 0);
    m = m / c;
    return;
}

template <typename T>
bool operator== (matrix2d<T> m1, matrix2d<T> m2) {
    if (m1.shape() != m2.shape()) {
        return false;
    }
    for (size_t i = 0; i < m1.shape()[0]; i++) {
        if (m1[i] != m2[i]) {
            return false;
        }
    }
    return true;
}

template <typename T>
bool operator!= (matrix2d<T> m1, matrix2d<T> m2) {
    return !(m1 == m2);
}

template <typename T>
std::ostream& operator<< (std::ostream& out, matrix2d<T> m) {
    for (size_t i = 0; i < m.shape()[0]; i++) {
        out << m[i] << std::endl;
    }
    return out;
}

#endif