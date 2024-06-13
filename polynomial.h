#include <cmath>
#include <utility>
#include <cassert>
#include "vector_operations.h"

#ifndef POLYNOMIAL
#define POLYNOMIAL
template <typename T>
class polynomial {
    private:
        std::vector<std::pair<double, double>> terms{};    // structure: <coeff, power>; power from large to small
        /* helper functions */
        void add_term(T coeff, T power) {
            size_t i = 0;
            for (; i < terms.size(); i++) {
                if (terms[i].second < power) {
                    break;
                } else if (terms[i].second == power) {
                    terms[i].first += coeff;
                    return;
                }
            }
            std::pair<double, double> p(static_cast<double>(coeff), static_cast<double>(power));
            terms.insert(terms.begin() + i, p);
            remove_zero_terms();
        }
        void add_term(std::pair<double, double> p) {
            add_term(p.first, p.second);
        }
        void remove_zero_terms() {
            std::vector<std::pair<double, double>>::iterator it = terms.begin();
            while (it < terms.end() && (terms.size() > 1)) {
                // std::cout << it->first << ' ' << it->second << std::endl;
                if (it->first == 0) {
                    it = terms.erase(it);
                } else {
                    it++;
                }
            }
        }
    public:
        /* constructors */
        polynomial(std::vector<std::pair<double, double>> terms) {
            for (std::pair<double, double> i : terms) {
                add_term(i);
            }
            if (this->terms.size() == 0) {
                add_term(0, 0);
            }
            remove_zero_terms();
        }
        polynomial(std::vector<double> terms) {
            for (size_t i = 0; i < terms.size(); i++) {
                add_term(terms[i], terms.size() - 1 - i);
            }
            if (this->terms.size() == 0) {
                add_term(0, 0);
            }
            remove_zero_terms();
        }
        polynomial() {
            add_term(0, 0);
        }
        double calculate(T x) {
            double sum = 0;
            for (std::pair<double, double> term : terms) {
                sum += term.first * pow(x, term.second);
            }
            return sum;
        }
        /* operators overloading */
        void operator= (polynomial<T> p) {
            terms = p.terms;
        }
        friend std::ostream& operator<< (std::ostream& out, polynomial<T> p) {
            bool leading_zeros = true;
            for (size_t i = 0; i < p.terms.size(); i++) {
                if ((p.terms[i].first == 0) && (i < p.terms.size() - 1)) {
                    continue;
                } else if ((p.terms[i].first > 0) && (i > 0) && !(leading_zeros)) {
                    out << "+";
                }
                leading_zeros = false;
                if ((p.terms[i].first == -1) && (p.terms[i].second != 0)) {
                    out << "-";
                } else if (p.terms[i].first == -0) {
                    out << "0";
                } else if ((p.terms[i].first != 1) || (((p.terms[i].first == -1) || (p.terms[i].first == 1)) && (p.terms[i].second == 0))) {
                    out << p.terms[i].first;
                }
                if ((p.terms[i].first == 0) && (i >= p.terms.size() - 1)) {
                    continue;
                }
                if (p.terms[i].second == 1) {
                    out << "x";
                } else if (p.terms[i].second != 0) {
                    out << "x^" << ((p.terms[i].second < 0) ? "(" : "") << p.terms[i].second << ((p.terms[i].second < 0) ? ")" : "");
                }
            }
            return out;
        }
        std::pair<double, double>& operator[] (size_t index) {
            assert(index < terms.size());
            assert(index >= 0);
            return terms[index];
        }
        polynomial<T> operator+ (polynomial<T> p) {
            polynomial new_p;
            for (std::pair t : terms) {
                new_p.add_term(t);
            }
            for (std::pair t : p.terms) {
                new_p.add_term(t);
            }
            return new_p;
        }
        polynomial<T> operator+ (T c) {
            polynomial new_p(terms);
            new_p.add_term(c, 0);
            return new_p;
        }
        polynomial<T> operator- (polynomial<T> p) {
            polynomial new_p;
            for (std::pair t : terms) {
                t.first *= -1;
                new_p.add_term(t);
            }
            for (std::pair t : p.terms) {
                t.first *= -1;
                new_p.add_term(t);
            }
            return new_p;
        }
        polynomial<T> operator- (T c) {
            polynomial new_p(terms);
            new_p.add_term(-c, 0);
            return new_p;
        }
        polynomial<T> operator* (T c) {
            polynomial new_p(terms);
            for (size_t i = 0; i < terms.size(); i++) {
                new_p[i].first *= c;
            }
            new_p.remove_zero_terms();
            return new_p;
        }
        polynomial<T> operator* (polynomial<T> p) {
            polynomial new_p;
            for (std::pair t : terms) {
                for (size_t i = 0; i < p.num_terms(); i++) {
                    new_p.add_term(t.first * p[i].first, t.second + p[i].second);
                }
            }
            return new_p;
        }
        polynomial<T> operator/ (T c) {
            assert(c != 0);
            polynomial new_p(terms);
            for (size_t i = 0; i < terms.size(); i++) {
                new_p[i].first /= c;
            }
            new_p.remove_zero_terms();
            return new_p;
        }
        /* functions */
        void differentiate(size_t order = 1) {
            for (size_t o = 0; o < order; o++) {
                for (size_t i = 0; i < terms.size(); i++) {
                    terms[i].first *= terms[i].second;
                    if (terms[i].second != 0) {
                        terms[i].second--;
                    }
                }
                remove_zero_terms();
            }
        }
        void integrate(size_t order = 1) {
            for (size_t o = 0; o < order; o++) {
                for (size_t i = 0; i < terms.size(); i++) {
                    assert(terms[i].second != -1);
                    terms[i].second++;
                    terms[i].first /= terms[i].second;
                }
                remove_zero_terms();
            }
        }
        size_t num_terms() {
            return terms.size();
        }
};

template <typename T1, typename T2>
std::ostream& operator<< (std::ostream& out, std::pair<T1, T2> p) {
    out << p.first << ' ' << p.second;
    return out;
}

#endif