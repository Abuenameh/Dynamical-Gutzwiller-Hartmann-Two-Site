/* 
 * File:   mathematica.hpp
 * Author: Abuenameh
 *
 * Created on November 6, 2011, 12:58 PM
 */

#ifndef MATHEMATICA_HPP
#define	MATHEMATICA_HPP

template<typename T>
class mathematica {
public:

    mathematica(T& v_) : v(v_) {
    }

    T& v;
};

template<>
class mathematica<int> {
public:

    mathematica(int i_) : i(i_) {
    }

    int i;
};

template<>
class mathematica<double> {
public:

    mathematica(double d_) : d(d_) {
    }

    double d;
};

template<>
class mathematica<bool> {
public:

    mathematica(bool b_) : b(b_) {
    }

    bool b;
};

template<>
class mathematica<complex<double> > {
public:

    mathematica(complex<double> c_) : c(c_) {
    }

    complex<double> c;
};

ostream& operator<<(ostream& out, const mathematica<int> m) {
    out << m.i;
    return out;
}

ostream& operator<<(ostream& out, const mathematica<double> m) {
    double d = m.d;
    ostringstream oss;
    oss << setprecision(numeric_limits<double>::digits10) << d;
    out << replace_all_copy(oss.str(), "e", "*^");
    return out;
}

ostream& operator<<(ostream& out, const mathematica<bool> m) {
    out << (m.b ? "True" : "False");
    return out;
}

ostream& operator<<(ostream& out, const mathematica<complex<double> > m) {
    complex<double> c = m.c;
    out << "(" << mathematica<double>(c.real()) << ")+I(" << mathematica<double>(c.imag()) << ")";
    return out;
}

template<typename Scalar, int Rows>
ostream& operator<<(ostream& out, const mathematica<Matrix<Scalar, Rows, 1 > >& m) {
    Matrix<Scalar, Rows, 1 > & v = m.v;
    out << "{" << mathematica<Scalar > (v[0]);
    for (int i = 1; i < Rows; i++) {
        out << "," << mathematica<Scalar > (v[i]);
    }
    out << "}";
    return out;
}

template<typename Scalar, int Rows, int Cols>
ostream& operator<<(ostream& out, const mathematica<Matrix<Scalar, Rows, Cols> >& m) {
    Matrix<Scalar, Rows, Cols>& mat = m.v;
    out << "{{" << mathematica<Scalar > (mat(0, 0));
    for (int j = 1; j < Cols; j++) {
        out << "," << mathematica<Scalar > (mat(0, j));
    }
    out << "}";
    for (int i = 1; i < Rows; i++) {
        out << ",{" << mathematica<Scalar > (mat(i, 0));
        for (int j = 1; j < Cols; j++) {
            out << "," << mathematica<Scalar > (mat(i, j));
        }
        out << "}";
    }
    out << "}";
    return out;
}

template<typename Scalar, int Rows>
ostream& operator<<(ostream& out, const mathematica<Array<Scalar, Rows, 1 > >& m) {
    Array<Scalar, Rows, 1 > & v = m.v;
    out << "{" << mathematica<Scalar > (v[0]);
    for (int i = 1; i < Rows; i++) {
        out << "," << mathematica<Scalar > (v[i]);
    }
    out << "}";
    return out;
}

template<typename Scalar, int Rows, int Cols>
ostream& operator<<(ostream& out, const mathematica<Array<Scalar, Rows, Cols> >& m) {
    Array<Scalar, Rows, Cols>& mat = m.v;
    out << "{{" << mathematica<Scalar > (mat(0, 0));
    for (int j = 1; j < Cols; j++) {
        out << "," << mathematica<Scalar > (mat(0, j));
    }
    out << "}";
    for (int i = 1; i < Rows; i++) {
        out << ",{" << mathematica<Scalar > (mat(i, 0));
        for (int j = 1; j < Cols; j++) {
            out << "," << mathematica<Scalar > (mat(i, j));
        }
        out << "}";
    }
    out << "}";
    return out;
}

template<typename Scalar>
ostream& operator<<(ostream& out, const mathematica<Array<Scalar, Dynamic, 1 > >& m) {
    Array<Scalar, Dynamic, 1 > & a = m.v;
    out << "{" << mathematica<Scalar > (a[0]);
    for (int i = 1; i < a.size(); i++) {
        out << "," << mathematica<Scalar > (a[i]);
    }
    out << "}";
    return out;
}

template<typename T, typename Alloc>
ostream& operator<<(ostream& out, const mathematica<deque<T, Alloc> >& m) {
    deque<T, Alloc>& d = m.v;
    out << "{" << mathematica<T > (d[0]);
    for (int i = 1; i < d.size(); i++) {
        out << "," << mathematica<T > (d[i]);
    }
    out << "}";
    return out;
}

template<typename T, typename Alloc>
ostream& operator<<(ostream& out, const mathematica<multi_array<T, 2, Alloc > >& m) {
    multi_array<T, 2, Alloc> & ma = m.v;
    int r = ma.shape()[0];
    int c = ma.shape()[1];
    out << "{{" << mathematica<T > (ma[0][0]);
    for (int j = 1; j < c; j++) {
        out << "," << mathematica<T > (ma[0][j]);
    }
    out << "}";
    for (int i = 1; i < r; i++) {
        out << ",{" << mathematica<T > (ma[i][0]);
        for (int j = 1; j < c; j++) {
            out << "," << mathematica<T > (ma[i][j]);
        }
        out << "}";
    }
    out << "}";
    return out;
}

template<typename T, typename Alloc>
ostream& operator<<(ostream& out, const mathematica<multi_array<T, 3, Alloc > >& m) {
    multi_array<T, 3, Alloc> & ma = m.v;
    int r = ma.shape()[0];
    int c = ma.shape()[1];
    int d = ma.shape()[2];
    out << "{{{" << mathematica<T > (ma[0][0][0]);
    for (int k = 1; k < d; k++) {
        out << "," << mathematica<T > (ma[0][0][k]);
    }
    out << "}";
    for (int j = 1; j < c; j++) {
        out << ",{" << mathematica<T > (ma[0][j][0]);
        for (int k = 1; k < d; k++) {
            out << "," << mathematica<T > (ma[0][j][k]);
        }
        out << "}";
    }
    out << "}";
    for (int i = 1; i < r; i++) {
        out << ",{{" << mathematica<T > (ma[i][0][0]);
        for (int k = 1; k < d; k++) {
            out << "," << mathematica<T > (ma[i][0][k]);
        }
        out << "}";
        for (int j = 1; j < c; j++) {
            out << ",{" << mathematica<T > (ma[i][j][0]);
            for (int k = 1; k < d; k++) {
                out << "," << mathematica<T > (ma[i][j][k]);
            }
            out << "}";
        }
        out << "}";
    }
    out << "}";
    return out;

}

template<typename T, typename U>
ostream& operator<<(ostream& out, const mathematica<pair<T, U> >& p) {
    out << "{" << mathematica<T>(p.v.first) << "," << mathematica<U>(p.v.second) << "}";
    return out;
}

template<typename T, size_t N>
ostream& operator<<(ostream& out, const mathematica<boost::array<T, N> >& arr) {
    boost::array<T, N>& a = arr.v;
    out << "{" << mathematica<T>(a[0]);
    for (int i = 1; i < N; i++) {
        out << "," << mathematica<T>(a[i]);
    }
    out << "}";
    return out;
}

template<typename T>
mathematica<T> math(T& t) {
    return mathematica<T > (t);
}

mathematica<double> math(double d) {
    return mathematica<double>(d);
}

mathematica<complex<double> > math(complex<double> c) {
    return mathematica<complex<double> >(c);
}

#endif	/* MATHEMATICA_HPP */

