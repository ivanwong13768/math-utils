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
        double find_term(double power) {
            double coeff = 0;
            for (std::pair t : terms) {
                if (t.second == power) {
                    coeff = t.first;
                    break;
                }
            }
            return coeff;
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
            polynomial<T> new_p(terms);
            for (std::pair t : p.terms) {
                new_p.add_term(t);
            }
            new_p.remove_zero_terms();
            return new_p;
        }
        polynomial<T> operator+ (T c) {
            polynomial<T> new_p(terms);
            new_p.add_term(c, 0);
            return new_p;
        }
        polynomial<T> operator- (polynomial<T> p) {
            polynomial<T> new_p(terms);
            for (std::pair t : p.terms) {
                new_p.add_term(-t.first, t.second);
            }
            new_p.remove_zero_terms();
            return new_p;
        }
        polynomial<T> operator- (T c) {
            polynomial<T> new_p(terms);
            new_p.add_term(-c, 0);
            return new_p;
        }
        polynomial<T> operator* (T c) {
            polynomial<T> new_p(terms);
            for (size_t i = 0; i < terms.size(); i++) {
                new_p[i].first *= c;
            }
            new_p.remove_zero_terms();
            return new_p;
        }
        polynomial<T> operator* (polynomial<T> p) {
            polynomial<T> new_p;
            for (std::pair t : terms) {
                for (size_t i = 0; i < p.num_terms(); i++) {
                    new_p.add_term(t.first * p[i].first, t.second + p[i].second);
                }
            }
            return new_p;
        }
        polynomial<T> operator/ (T c) {
            assert(c != 0);
            polynomial<T> new_p(terms);
            for (size_t i = 0; i < terms.size(); i++) {
                new_p[i].first /= c;
            }
            new_p.remove_zero_terms();
            return new_p;
        }
        polynomial<T> operator/ (std::pair<T, T> t) {
            assert(t.first != 0);
            polynomial<T> new_p;
            for (std::pair term : terms) {
                new_p.add_term(term.first / t.first, term.second - t.second);
            }
            return new_p;
        }
        polynomial<T> operator% (polynomial<T> p) {
            polynomial<T> new_p(terms);
            for (size_t i = p.num_terms() - 1; i < num_terms(); i++) {
                polynomial<T> minus(p);
                polynomial<T> power(std::vector<std::pair<double, double>>{std::make_pair(1, static_cast<double>(terms[i].second))});
                minus *= power;
                // std::cout << "power: " << power << std::endl;
                // std::cout << "minus: " << minus << std::endl;
                double coeff = new_p.find_term(minus[0].second);
                // std::cout << "coeff: " << coeff << std::endl;
                if (coeff == 0) {
                    continue;
                }
                // std::cout << "product: " << minus * coeff << std::endl;
                new_p -= minus * coeff;
                // std::cout << "new_p: " << new_p << std::endl << std::endl;
            }
            new_p.remove_zero_terms();
            return new_p;
        }
        polynomial<T> operator% (T c) {
            polynomial const_p(std::vector<T>{c});
            return (*this % const_p);
        }
        polynomial<T> operator/ (polynomial<T> p) {
            for (std::pair t : p.terms) {
                assert(t.first != 0);
            }
            polynomial<T> new_p;
            polynomial<T> temp(terms);
            polynomial<T> remainder = *this % p;
            // std::cout << "remainder: " << remainder << std::endl;
            if (!((remainder.num_terms() == 1) && (remainder[0].first == 0))) {
                std::cerr << "Warning: remainder truncated -> ";
            }
            // temp -= remainder;
            temp.remove_zero_terms();
            // std::cout << "debug: " << temp << ' ' << new_p << std::endl;
            size_t loop_times = temp.num_terms();
            for (size_t i = p.num_terms() - 1; i <= loop_times; i++) {
                polynomial<T> minus(p);
                polynomial<T> power(std::vector<std::pair<double, double>>{std::make_pair(1, temp.terms[0].second - p[0].second)});
                minus *= power;
                // std::cout << "power: " << power << std::endl;
                // std::cout << "minus: " << minus << std::endl;
                double coeff = temp.find_term(minus[0].second);
                // std::cout << "coeff: " << coeff << std::endl;
                if (coeff == 0) {
                    continue;
                }
                temp -= minus * coeff;
                new_p.add_term(coeff, power[0].second);
                // std::cout << "debug: " << temp << ' ' << new_p << std::endl << std::endl;
            }
            new_p.remove_zero_terms();
            return new_p;
        }
        bool operator== (polynomial<T> p) {
            if (num_terms() != p.num_terms()) {
                return false;
            }
            for (size_t i = 0; i < num_terms(); i++) {
                if ((terms[i].first != p[i].first) || (terms[i].second != p[i].second)) {
                    return false;
                }
            }
            return true;
        }
        bool operator== (T c) {
            if ((num_terms() != 1) || terms[0].first != c) {
                return false;
            }
            return true;
        }
        bool operator!= (polynomial<T> p) {
            return !(*this == p);
        }
        bool operator!= (T c) {
            return !(*this == c);
        }
        bool operator> (polynomial<T> p) {
            /* long double min_power = std::min(terms[num_terms() - 1].second, p[p.num_terms() - 1].second);
            long double max_power = std::max(terms[0].second, p[0].second); */
            // TODO: make a set that contains all power -> loop through the whole set and compare
            for (long double i = max_power; i >= min_power; i--) {
                if (this->find_term(i) < p.find_term(i)) {
                    return false;
                }
                std::cout << i << ' ';
            }
            return true;
        }
        bool operator<= (polynomial<T> p) {
            return !(*this > p);
        }
        bool operator>= (polynomial<T> p) {
            return ((*this > p) || (*this == p));
        }
        bool operator< (polynomial<T> p) {
            return ((*this <= p) && (*this != p));
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

template <typename T>
polynomial<T> operator+ (T c, polynomial<T> p) {
    return (p + c);
}

template <typename T>
polynomial<T> operator- (T c, polynomial<T> p) {
    return (p - c);
}

template <typename T>
polynomial<T> operator* (T c, polynomial<T> p) {
    return (p * c);
}

template <typename T>
polynomial<T> operator/ (T c, polynomial<T> p) {
    return (p / c);
}

template <typename T>
void operator+= (polynomial<T>& p1, polynomial<T> p2) {
    p1 = p1 + p2;
}

template <typename T>
void operator+= (polynomial<T>& p, T c) {
    p = p + c;
}

template <typename T>
void operator-= (polynomial<T>& p1, polynomial<T> p2) {
    p1 = p1 - p2;
}

template <typename T>
void operator-= (polynomial<T>& p, T c) {
    p = p - c;
}

template <typename T>
void operator*= (polynomial<T>& p1, polynomial<T> p2) {
    p1 = p1 * p2;
}

template <typename T>
void operator*= (polynomial<T>& p, T c) {
    p = p * c;
}

template <typename T>
void operator/= (polynomial<T>& p1, polynomial<T> p2) {
    p1 = p1 / p2;
}

template <typename T>
void operator/= (polynomial<T>& p, T c) {
    p = p / c;
}

template <typename T>
void operator%= (polynomial<T>& p1, polynomial<T> p2) {
    p1 = p1 % p2;
}

template <typename T>
bool operator== (T c, polynomial<T> p) {
    return (p == c);
}

template <typename T>
bool operator!= (T c, polynomial<T> p) {
    return (p != c);
}

#endif