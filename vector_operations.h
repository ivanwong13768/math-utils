#include <iostream>
#include <vector>

#ifndef VECTOR
#define VECTOR

template <typename T>
std::vector<T> operator+ (const std::vector<T>& v1, const std::vector<T>& v2) {
    std::vector<T> out{};
    if (v1.size() >= v2.size()) {
        for (int i = 0; i < v2.size(); i++) {
            out.push_back(v1[i] + v2[i]);
        }
        for (int i = v2.size(); i < v1.size(); i++) {
            out.push_back(v1[i]);
        }
    } else {
        for (int i = 0; i < v1.size(); i++) {
            out.push_back(v1[i] + v2[i]);
        }
        for (int i = v1.size(); i < v2.size(); i++) {
            out.push_back(v2[i]);
        }
    }
    return out;
}

template <typename T>
std::vector<T> operator- (const std::vector<T>& v1, const std::vector<T>& v2) {
    std::vector<T> out{};
    if (v1.size() >= v2.size()) {
        for (int i = 0; i < v2.size(); i++) {
            out.push_back(v1[i] - v2[i]);
        }
        for (int i = v2.size(); i < v1.size(); i++) {
            out.push_back(v1[i]);
        }
    } else {
        for (int i = 0; i < v1.size(); i++) {
            out.push_back(v1[i] - v2[i]);
        }
        for (int i = v1.size(); i < v2.size(); i++) {
            out.push_back(-v2[i]);
        }
    }
    return out;
}

template <typename T>
std::vector<T> operator* (const std::vector<T>& v1, const std::vector<T>& v2) {
    std::vector<T> out{};
    if (v1.size() >= v2.size()) {
        for (int i = 0; i < v2.size(); i++) {
            out.push_back(v1[i] * v2[i]);
        }
        for (int i = v2.size(); i < v1.size(); i++) {
            out.push_back(0);
        }
    } else {
        for (int i = 0; i < v1.size(); i++) {
            out.push_back(v1[i] * v2[i]);
        }
        for (int i = v1.size(); i < v2.size(); i++) {
            out.push_back(0);
        }
    }
    return out;
}

template <typename T>
std::vector<T> operator* (const std::vector<T>& v1, const double& c) {
    std::vector<T> out{};
    for (int i = 0; i < v1.size(); i++) {
        out.push_back(v1[i] * c);
    }
    return out;
}

template <typename T>
std::vector<T> operator/ (const std::vector<T>& v1, const std::vector<T>& v2) {
    std::vector<T> out{};
    if (v1.size() >= v2.size()) {
        for (int i = 0; i < v2.size(); i++) {
            out.push_back(v1[i] / v2[i]);
        }
        for (int i = v2.size(); i < v1.size(); i++) {
            out.push_back(INFINITY);    // placeholder to keep length constant
        }
    } else {
        for (int i = 0; i < v1.size(); i++) {
            out.push_back(v1[i] / v2[i]);
        }
        for (int i = v1.size(); i < v2.size(); i++) {
            out.push_back(0);
        }
    }
    return out;
}

template <typename T>
std::vector<T> operator/ (const std::vector<T>& v1, const double& c) {
    std::vector<T> out{};
    for (int i = 0; i < v1.size(); i++) {
        out.push_back(v1[i] / c);
    }
    return out;
}

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& vec) {
    for (T i : vec) {
        out << i << '\t';
    }
    return out;
}

template <typename T>
T vector_sum(const std::vector<T>& v) {
    T sum = 0;
    for (T i : v) {
        sum += i;
    }
    return sum;
}

template <typename T>
T vector_mul(const std::vector<T>& v) {
    T sum = 1;
    for (T i : v) {
        sum *= i;
    }
    return sum;
}

#endif