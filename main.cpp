/* 
 * File:   main.cpp
 * Author: Abuenameh
 *
 * Created on November 4, 2011, 12:22 AM
 */

//#include <windows.h>

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <complex>
#include <deque>
#include <string>
#include <sstream>
#include <iomanip>
#include <limits>
#include <fstream>
#include <utility>
#include <map>
#include <cmath>

//using namespace std;

using std::complex;
using std::cout;
using std::cerr;
using std::endl;
using std::deque;
using std::ostream;
using std::numeric_limits;
using std::setprecision;
using std::string;
using std::ostringstream;
using std::allocator;
using std::flush;
using std::pair;
using std::make_pair;
using std::max;
using std::map;
using std::ptr_fun;

#define BOOST_THREAD_USE_LIB

#define BOOST_DISABLE_ASSERTS

#include <boost/lexical_cast.hpp>
#include <boost/multi_array.hpp>
#include <boost/algorithm/string.hpp>
#define timer   timer_class
#include <boost/progress.hpp>
#undef timer
#include <boost/random.hpp>
#include <boost/thread.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/ref.hpp>
#include <boost/bind.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>
#include <boost/timer/timer.hpp>

using boost::lexical_cast;
using boost::multi_array;
using boost::extents;
using boost::progress_display;
using boost::random::mt19937;
using boost::random::uniform_real_distribution;
using boost::filesystem::path;
using boost::filesystem::exists;
using boost::filesystem::ofstream;
using boost::thread;
using boost::ref;
using boost::reference_wrapper;
using boost::bind;
using boost::mutex;
using boost::lock_guard;
using boost::thread_group;

using namespace boost::algorithm;

#include <Eigen/Dense>

using Eigen::Matrix;
using Eigen::Array;
using Eigen::Dynamic;
using Eigen::aligned_allocator;

#define L 2
#define nmax 7
#define dim (nmax+1)

complex<double> I = complex<double>(0, 1);

typedef Matrix<double, dim, 1 > vector;
typedef Matrix<double, dim, dim> matrix;
typedef Array<double, L, 1 > array;
typedef Matrix<double, L, L> matrixL;
typedef Matrix<double, L, 1 > vectorL;

typedef boost::array<complex<double>, dim> fiarray;
typedef boost::array<fiarray, L> farray;
typedef boost::array<complex<double>, L> carrayL;
typedef boost::array<double, L> arrayL;
typedef deque<farray> farraydeque;
typedef boost::array<double, L*L> Jmat;
typedef boost::array<double, dim> arraydim;
typedef boost::array<double, L*L> matrixLL;

#define ij(i, j) L * i + j

/*typedef Matrix<complex<double>, dim, 1 > cvector;
typedef Matrix<complex<double>, dim, dim> cmatrix;
typedef Array<complex<double>, L, 1 > carrayL;*/
typedef Matrix<complex<double>, L, L> cmatrixL;

/*typedef Array<double, Dynamic, 1 > darray;
typedef Array<double, Dynamic, Dynamic > darray2d;*/

typedef map<pair<int, int>, double> pairmap;

/*typedef multi_array<pairmap, 2 > multimap;
typedef multi_array<array, 2, aligned_allocator<array> > multiarray;

typedef deque<pairmap> mapdeque;
typedef deque<array, aligned_allocator<array> > arraydeque;

typedef multi_array<vector, 2, aligned_allocator<vector> > multivector;*/

typedef Array<double, dim, 1 > rarray;
typedef Matrix<complex<double>, dim, 1 > cvector;
typedef Array<complex<double>, dim, 1 > carray;
typedef Matrix<complex<double>, dim, dim> cmatrix;
typedef Array<complex<double>, L, 1 > aarray;

/*typedef Array<complex<double>, L*dim, 1 > farray;
typedef deque<farray, aligned_allocator<farray> > farraydeque;
typedef Array<double, L*dim, 1 > tarray;
typedef deque<tarray, aligned_allocator<tarray> > tarraydeque;

typedef deque<carray, aligned_allocator<carray> > carraydeque;*/

typedef boost::array<deque<complex<double> >, L*dim> fpoints;
typedef boost::array<deque<complex<double> >, L> funcpoints;

/*typedef deque<matrix, aligned_allocator<matrix> > matrixdeque;*/

using Eigen::SelfAdjointEigenSolver;
using Eigen::Success;

#include "mathematica.hpp"

template<typename T> void printMath(ostream& out, string name, T& t) {
    out << name << "=" << ::math(t) << ";" << endl;
}

template<typename T> void printMath(ostream& out, string name, int i, T& t) {
    out << name << "[" << i << "]" << "=" << ::math(t) << ";" << endl;
}

inline int mod(int i) {
    return (i + L) % L;
}

struct operators {
    matrix& ni;
    matrix& bai;
    matrix& bci;
    matrix& n2i;
    matrix& n2ni;
};

struct inputs {
    //    double W;
    double mu;
    array& U;
    pairmap& J;

    double eps;
    double theta;
};

struct results {
    double& fc;
    double& fs;

    farray& f0;
};

double chopeps = 0;

double chop(double x) {
    return (fabs(x) < chopeps) ? 0 : x;
}

complex<double> ai(fiarray& fi) {
    complex<double> a = 0;
    for (int n = 1; n <= nmax; n++) {
        a += sqrt(n) * conj(fi[n - 1]) * fi[n];
    }
    return a;
}

double N = 1000;
double g13 = 2.5e9;
double g24 = 2.5e9;
double delta = 1.0e12;
double Delta = -2.0e10;
double alpha = 1.1e7;

double g = sqrt(N) * g13;
double g2 = N * g13 * g13;
double g242 = g24 * g24;

double JW(double W) {
    return alpha * (W * W) / (g * g + W * W);
}

double JWij(double Wi, double Wj) {
    return alpha * (Wi * Wj) / (sqrt(g * g + Wi * Wi) * sqrt(g * g + Wj * Wj));
}

void JW(double W2, matrixLL& xij, matrixLL& J) {
    double aW2 = alpha * W2;
    arrayL gW;
    for (int i = 0; i < L; i++) {
        gW[i] = sqrt(g2 + W2 * xij[ij(i, i)]);
    }
    for (int j = 0; j < L; j++) {
        int l1 = mod(j - 1);
        int l2 = mod(j + 1);

        J[ij(j, l1)] = aW2 * xij[ij(j, l1)] / (gW[j] * gW[l1]);
        J[ij(j, l2)] = aW2 * xij[ij(j, l2)] / (gW[j] * gW[l2]);
    }
}

array JW(array W) {
    array v = W / sqrt(g * g + W * W);
    array J;
    for (int i = 0; i < L - 1; i++) {
        J[i] = v[i] * v[i + 1];
    }
    J[L - 1] = v[L - 1] * v[0];
    J *= alpha;
    return J;
}

void JW(arrayL& W, matrixLL& J) {
    for (int j = 0; j < L; j++) {
        int l1 = mod(j - 1);
        int l2 = mod(j + 1);

        J[ij(j, l1)] = JWij(W[j], W[l1]);
        J[ij(j, l2)] = JWij(W[j], W[l2]);
    }

}

double UW(double W) {
    return -(g24 * g24) / Delta * (g * g * W * W) / ((g * g + W * W) * (g * g + W * W));
}

array UW(array W) {
    return -(g24 * g24) / Delta * (g * g * W * W) / ((g * g + W * W) * (g * g + W * W));
}

void UW(arrayL& W, arrayL& U) {
    for (int i = 0; i < L; i++) {
        U[i] = -(g24 * g24) / Delta * (g * g * W[i] * W[i]) / ((g * g + W[i] * W[i]) * (g * g + W[i] * W[i]));
    }
}

void UW(double W2, matrixLL& xij, arrayL& U) {
    for (int i = 0; i < L; i++) {
        double xii = xij[ij(i, i)];
        double W2i = W2 * xii;
        double gW = g2 + W2i;
        U[i] = -g242 / Delta * (g2 * W2i) / (gW * gW);
    }
}

class Wfunc {
public:

    virtual double operator()(double t) {
        return 0;
    }
};

class dfi {
public:

    dfi(int i_, carrayL& a_, arrayL& xi_, Wfunc* Wt_, double mu_) : i(i_), a(a_), xi(xi_), Wt(Wt_), mu(mu_) {
    }

    void operator()(double t, fiarray& f, fiarray& dfdt) {
        double W = (*Wt)(t);
        double U = UW(W * xi[i]);
        complex<double> aj = JWij(W*xi[i], W*xi[mod(i + 1)]) * a[mod(i + 1)];
        complex<double> ajc = conj(aj);
        dfdt[0] = -I * (-1.0 * ajc * f[1]);
        dfdt[1] = -I * (-mu * f[1] - aj * f[0] - sqrt(2) * ajc * f[2]);
        dfdt[2] = -I * ((2 * U - 2 * mu) * f[2] - sqrt(2) * aj * f[1] - sqrt(3) * ajc * f[3]);
        dfdt[3] = -I * ((6 * U - 3 * mu) * f[3] - sqrt(3) * aj * f[2] - 2.0 * ajc * f[4]);
        dfdt[4] = -I * ((12 * U - 4 * mu) * f[4] - 2.0 * aj * f[3] - sqrt(5) * ajc * f[5]);
        dfdt[5] = -I * ((20 * U - 5 * mu) * f[5] - sqrt(5) * aj * f[4] - sqrt(6) * ajc * f[6]);
        dfdt[6] = -I * ((30 * U - 6 * mu) * f[6] - sqrt(6) * aj * f[5] - sqrt(7) * ajc * f[7]);
        dfdt[7] = -I * ((42 * U - 7 * mu) * f[7] - sqrt(7) * aj * f[6]);
    }

private:

    int i;
    carrayL& a;
    arrayL& xi;
    Wfunc* Wt;
    double mu;
};

void rk4(fiarray& f, fiarray& dfdt, double t, double h, fiarray& fout, dfi& derivs) {
    double th, hh, h6;
    fiarray dfn, dft, ft;

    hh = 0.5 * h;
    h6 = h / 6.0;
    th = t + hh;
    for (int i = 0; i < dim; i++) {
        ft[i] = f[i] + hh * dfdt[i];
    }
    derivs(th, ft, dft);
    for (int i = 0; i < dim; i++) {
        ft[i] = f[i] + hh * dft[i];
    }
    derivs(th, ft, dfn);
    for (int i = 0; i < dim; i++) {
        ft[i] = f[i] + h * dfn[i];
    }
    for (int i = 0; i < dim; i++) {
        dfn[i] += dft[i];
    }
    derivs(t + h, ft, dft);
    for (int i = 0; i < dim; i++) {
        fout[i] = f[i] + h6 * (dfdt[i] + dft[i] + 2.0 * dfn[i]);
    }
}

void find_state(double eps, double J, array& U, double mu, matrix& n2ni, matrix& ni, matrix& bai, farray& f0) {

    array a = array::Zero();
    array n = array::Zero();
    array n2n = array::Zero();
    array anew = array::Zero();
    double Da = numeric_limits<double>::infinity();
    array da = array::Zero();


    matrixL rho = matrixL::Zero();

    vector ffin;

    a = array::Constant(0.5);
    Da = numeric_limits<double>::infinity();
    while (Da > eps) {
        for (int i = 0; i < L; i++) {
            int j1 = mod(i - 1);
            int j2 = mod(i + 1);

            double aj = J * a[mod(i + 1)];
            matrix H = U[i] * n2ni - mu * ni - aj * (bai.transpose() + bai);

            SelfAdjointEigenSolver<matrix> es(H);
            if (es.info() != Success) {
                cerr << "Diagonalization failed!" << endl << "J = " << J << "\t" << "mu = " << mu << endl;
                exit(1);
            }

            vector fi = es.eigenvectors().col(0);
            fiarray& f0i = f0[i];
            for (int m = 0; m <= nmax; m++) {
                f0i[m] = fi[m];
            }

            anew[i] = fi.adjoint() * bai * fi;
            n[i] = (fi.adjoint() * ni * fi).value();
            n2n[i] = (fi.adjoint() * n2ni * fi).value();
        }

        Da = fabs((a.abs() - anew.abs()).maxCoeff());
        a = anew;
    }

}

//double Wt(double t) {
//}

class Wconst : public Wfunc {
public:

    Wconst(double W0) {
        W = W0;
    }

    virtual double operator()(double t) {
        return W;
    }

private:
    double W;
};

class Wlin : public Wfunc {
public:

    Wlin(double W0, double W1, double t0) {
        Wi = W0;
        Wf = W1;
        tf = t0;
    }

    virtual double operator()(double t) {
        return (t > tf) ? Wf : Wi + t * (Wf - Wi) / tf;
    }

private:
    double Wi;
    double Wf;
    double tf;
};

class WS : public Wfunc {
public:

    WS(double W0, double W1, double t0, double t1, double t2) {
        Wi = W0;
        Wf = W1;
        ti2 = t0;
        tf2 = t1;
        tf = t2;
    }

    virtual double operator()(double t) {
        if (t < ti2) {
            return Wi;
        }
        double tc = 0.5 * (ti2 + tf2);
        if (t <= tc) {
            return 2 * ((Wf - Wi) / ((tf2 - ti2) * (tf2 - ti2))) * (t - ti2) * (t - ti2) + Wi;
        }
        if (t < tf2) {
            return -2 * ((Wf - Wi) / ((tf2 - ti2) * (tf2 - ti2))) * (t - tf2) * (t - tf2) + Wf;
        }
        //        if (t < tf) {
        return Wf;
        //        }
    }

private:
    double Wi;
    double Wf;
    double ti2;
    double tf2;
    double tf;
};

class Wpw : public Wfunc {
public:

    Wpw(double W0, double W1, double t0, double t1, double t2) {
        Wi = W0;
        Wf = W1;
        ti2 = t0;
        tf2 = t1;
        tf = t2;
    }

    virtual double operator()(double t) {
        if (t < ti2) {
            return Wi;
        }
        if (t < tf2) {
            return Wi + (t - ti2) * (Wf - Wi) / (tf2 - ti2);
        }
        return Wf;
    }

private:
    double Wi;
    double Wf;
    double ti2;
    double tf2;
    double tf;
};

class Wpwcyc : public Wfunc {
public:

    Wpwcyc(double W0, double W1, double t0, double t1, int ncyc, bool f) {
        Wi = W0;
        Wf = W1;
        ti = t0;
        tf = t1;
        n = ncyc;
        endf = f;
        if (endf) {
            dt = (tf - ti) / (2 * n + 1);
        } else {
            dt = (tf - ti) / (2 * n);
        }
        m = (Wf - Wi) / dt;
    }

    virtual double operator()(double t) {
        if (t < ti) {
            return Wi;
        }
        if (t > tf) {
            return endf ? Wf : Wi;
        }
        int i1 = 2 * n;
        if (endf) {
            i1++;
        }
        for (int i = 0; i < i1; i++) {
            if (t < (ti + (i + 1) * dt)) {
                if (i % 2 == 0) {
                    return m * (t - (ti + i * dt)) + Wi;
                } else {
                    return -m * (t - (ti + i * dt)) + Wf;
                }
            }
        }
        return Wf;
    }

private:
    double Wi;
    double Wf;
    double ti;
    double tf;
    int n;
    bool endf;
    double m;
    double dt;
};

class Wpwsqrt : public Wfunc {
public:

    Wpwsqrt(double W0, double W1, double t0, double t1, double t2) {
        Wi = W0;
        Wf = W1;
        ti2 = t0;
        tf2 = t1;
        tf = t2;
        p = (Wf * Wf - Wi * Wi) / (tf2 - ti2);
        q = (tf2 * Wi * Wi - ti2 * Wf * Wf) / (tf2 - ti2);
    }

    virtual double operator()(double t) {
        if (t < ti2) {
            return Wi;
        }
        if (t < tf2) {
            return sqrt(p * t + q);
            //            return Wi + (t - ti2) * (Wf - Wi) / (tf2 - ti2);
        }
        return Wf;
    }

private:
    double Wi;
    double Wf;
    double ti2;
    double tf2;
    double tf;
    double p;
    double q;
};

class Wramp : public Wfunc {
public:

    Wramp(double Wi_, double Wf_, double tau_) : Wi(Wi_), Wf(Wf_), tau(tau_) {
    }

    double operator()(double t) {
        if (t <= tau) {
            return Wi + (Wf - Wi) * t / tau;
        } else {
            return Wf + (Wi - Wf) * (t - tau) / tau;
        }
    }

private:
    double Wi;
    double Wf;
    double tau;
};

void fsolve(Wfunc* Wt, double mu, arrayL& xi, farray& fstart, /*double t1,*/ double tf, int nstep, deque<double>& tt, deque<farray>& ft/*, matrixLL& xij*/) {
    double t, h;
    fiarray fout, df, fi;
    farray f;
    carrayL a;

    f = fstart;
    t = 0;
    ft.push_back(f);
    tt.push_back(t);
    //    ft[0] = f;
    //    tt[0] = t1;
    h = tf / nstep;
    for (int i = 0; i < L; i++) {
        fiarray fi = f[i];
        a[i] = ai(fi);
    }
    //    progress_display progress(nstep);

    for (int k = 1; k <= nstep; k++) {
        for (int i = 0; i < L; i++) {
            int site = (k % 2 == 0) ? i : L - 1 - i;

            dfi derivs(site, a, xi, Wt, mu);
            fiarray& fi = f[site];
            derivs(t, fi, df);
            rk4(fi, df, t, h, fout, derivs);
            a[site] = ai(fout);
            f[site] = fout;
        }
        t += h;
        //        ++progress;
    }
    tt.push_back(t);
    ft.push_back(f);
}

/*
 * 
 */
int main(int argc, char** argv) {

    //        cout << Eigen::SimdInstructionSetsInUse() << endl;
    //
    //    return 0;

    mt19937 rng;
    uniform_real_distribution<> uni(-1, 1);

    int seed = lexical_cast<int>(argv[1]);

    double Wi = lexical_cast<double>(argv[2]);
    double Wf = lexical_cast<double>(argv[3]);

    double mu = lexical_cast<double>(argv[4]);

    double Ui = UW(Wi);

    double D = lexical_cast<double>(argv[5]);

    double tau = lexical_cast<double>(argv[6]);
    double tf = 2 * tau;
    double dt = lexical_cast<double>(argv[7]);
    int dnsav = lexical_cast<int>(argv[8]);

    int nsteps = (int) ceil(2 * tau / dt);
    double h = 2 * tau / nsteps;


    int resi = 0;
    if (argc > 9) {
        resi = lexical_cast<int>(argv[9]);
    }

    double theta = 0;
    double eps = 1e-8;

#ifdef AMAZON
    path resdir("/home/ubuntu/Results/Dynamical Gutzwiller Hartmann Two Site");
#else
    path resdir("/Users/Abuenameh/Documents/Simulation Results/Dynamical Gutzwiller Hartmann Two Site Uniform");
    //        path resdir("/Users/Abuenameh/Documents/Simulation Results/Dynamical Gutzwiller Hartmann Comparison");
#endif
    if (!exists(resdir)) {
        cerr << "Results directory " << resdir << " does not exist!" << endl;
        exit(1);
    }
    ostringstream oss;
    oss << "res." << resi << ".txt";
    path resfile = resdir / oss.str();
#ifndef AMAZON
//    while (exists(resfile)) {
//        resi++;
//        oss.str("");
//        oss << "res." << resi << ".txt";
//        resfile = resdir / oss.str();
//    }
#endif
    //        if (seed == -1) {
    //            resi = -1;
    if (seed < 0) {
        resi = seed;
        oss.str("");
        oss << "res." << resi << ".txt";
        resfile = resdir / oss.str();
    }
    array xi = array::Constant(1);
    rng.seed(seed);
    if (seed > -1) {
        for (int j = 0; j < L; j++) {
            xi[j] = (1 + D * uni(rng));
        }
    }
    /*matrixLL xij;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            xij[ij(i, j)] = xi[i] * xi[j];
        }
    }*/

    //        double Ui = UW(Wi);
    double mui = mu * 2 * Ui;

    ofstream os(resfile);
    printMath(os, "seed", resi, seed);
    //        printMath(os, "theta", resi, theta);
    //        printMath(os, "eps", resi, eps);
    printMath(os, "Delta", resi, D);
    printMath(os, "Wres", resi, Wi);
    printMath(os, "mures", resi, mui);
    printMath(os, "Ures", resi, Ui);
    printMath(os, "xires", resi, xi);
    os << flush;

    printMath(os, "taures", resi, tau);
    printMath(os, "dtres", resi, dt);
    printMath(os, "dnsavres", resi, dnsav);
    os << flush;

    printMath(os, "Wires", resi, Wi);
    printMath(os, "Wfres", resi, Wf);
    os << flush;

    //        cout << "Res: " << resi << endl;

    matrix ni = matrix::Zero();
    matrix bai = matrix::Zero();
    matrix bci = matrix::Zero();

    for (int n = 0; n < dim; n++) {
        ni(n, n) = n;
        if (n > 0) {
            bai(n - 1, n) = sqrt(n);
        }
        if (n < nmax) {
            bci(n + 1, n) = sqrt(n + 1);
        }
    }

    matrix n2i = ni * ni;
    matrix n2ni = /*0.5 * */ni * (ni - matrix::Identity());


    array Um = array::Zero();

    array Wrand = Wi * xi;
    double Uav = 2 * Ui;
    Um = UW(Wrand) / Uav;
    double J = JWij(Wrand[0], Wrand[1]) / Uav;

    farray f0;

    find_state(eps, J, Um, mu, n2ni, ni, bai, f0);
    
    //        cout << endl << math(f0) << "//Transpose//MatrixForm" << endl << endl;
    //        exit(0);
    //        printMath(os, "f0res", resi, f0);

    /*double perturb = 0;
    if (argc > 18) {
        perturb = lexical_cast<double>(argv[18]);
    }

    if (perturb > 0) {
        for (int i = 0; i < L; i++) {
            //f0[i*dim + 0] += perturb;
        }
    }*/

    //        printMath(os, "f1res", resi, f0);

    //        printMath(os, "fc0res", resi, fc0);
    //        printMath(os, "fs0res", resi, fs0);

    deque<double> tt; //(nsteps + 1/*, 0*/);
    farraydeque ft; //(nsteps + 1/*, farray::Zero()*/);

    Wconst Wcons(Wi);
    Wramp Wt(Wi, Wf, tau);

    arrayL xiL;
    for (int i = 0; i < L; i++) {
        xiL[i] = xi[i];
    }
    fsolve(&Wt, mui, xiL, f0, tf, nsteps, tt, ft);
//    fsolve(&Wcons, mui, xiL, f0, tf, nsteps, tt, ft);

    cmatrix nic = ni.cast<complex<double> >();
    cmatrix baic = bai.cast<complex<double> >();
    cmatrix n2nic = n2ni.cast<complex<double> >();

    arrayL W0rand;
    for (int i = 0; i < L; i++) {
        W0rand[i] = Wi * xi[i];
    }
    arrayL U0;
    UW(W0rand, U0);
    double J0 = JWij(W0rand[0], W0rand[1]);

    carrayL a0;
    for (int i = 0; i < L; i++) {
        fiarray fi = f0[i];
        double norm = 0;
        for (int n = 0; n <= nmax; n++) {
            norm += std::norm(fi[n]);
        }
        a0[i] = 0;
        for (int n = 1; n <= nmax; n++) {
            a0[i] += sqrt(n) * conj(fi[n - 1]) * fi[n];
        }
        a0[i] /= norm;
    }

    double E0 = 0;
    arrayL n0;
    for (int i = 0; i < L; i++) {
        double aj = real(J0 * a0[mod(i + 1)]);

        fiarray fi = f0[i];
        arraydim fire;
        arraydim fisq;
        for (int j = 0; j < dim; j++) {
            fire[j] = real(fi[j]);
            fisq[j] = norm(fi[j]);
        }
        E0 += 2 * U0[i] * (fisq[2] + 3 * fisq[3] + 6 * fisq[4] + 10 * fisq[5] + 15 * fisq[6] + 21 * fisq[7]);
        n0[i] = fisq[1] + 2 * fisq[2] + 3 * fisq[3] + 4 * fisq[4] + 5 * fisq[5] + 6 * fisq[6] + 7 * fisq[7];
        E0 += -mui * n0[i];
        E0 += -2 * aj * (fire[0] * fire[1] + sqrt(2) * fire[1] * fire[2] + sqrt(3) * fire[2] * fire[3] + 2 * fire[3] * fire[4] + sqrt(5) * fire[4] * fire[5] + sqrt(6) * fire[5] * fire[6] + sqrt(7) * fire[6] * fire[7]);
        E0 += aj * real(a0[i]);
    }

    printMath(os, "E0Ures", resi, E0);
    E0 /= 2 * Ui;
    printMath(os, "E0res", resi, E0);

    //        double tf = tt.back();
    //        double Wf = (*Wt)(tf);
    farray ff = ft.back();
    arrayL Wfrand;
    for (int i = 0; i < L; i++) {
        Wfrand[i] = Wi * xi[i];
    }
    arrayL Uf;
    UW(Wfrand, Uf);
    double Jf = JWij(Wfrand[0], Wfrand[1]);

    carrayL af;
    for (int i = 0; i < L; i++) {
        fiarray fi = ff[i];
        double norm = 0;
        for (int n = 0; n <= nmax; n++) {
            norm += std::norm(fi[n]);
        }
        af[i] = 0;
        for (int n = 1; n <= nmax; n++) {
            af[i] += sqrt(n) * conj(fi[n - 1]) * fi[n];
        }
        af[i] /= norm;
    }

    complex<double> Ef = 0;
    arrayL nf;
    for (int i = 0; i < L; i++) {
        complex<double> aj = Jf * af[mod(i + 1)];
        complex<double> ajc = conj(aj);

        fiarray fi = ff[i];
        fiarray fic;
        nf[i] = 0;
        for (int n = 0; n <= nmax; n++) {
            fic[n] = conj(fi[n]);
            nf[i] += n * norm(fi[n]);
        }
        Ef += -fi[0] * (aj * fic[1]);
        Ef += -fi[1] * (mui * fic[1] + ajc * fic[0] + sqrt(2) * aj * fic[2]);
        Ef += -fi[2] * (2 * (mui - Uf[i]) * fic[2] + sqrt(2) * ajc * fic[1] + sqrt(3) * aj * fic[3]);
        Ef += -fi[3] * (3 * (mui - 2 * Uf[i]) * fic[3] + sqrt(3) * ajc * fic[2] + 2.0 * aj * fic[4]);
        Ef += -fi[4] * (4 * (mui - 3 * Uf[i]) * fic[4] + 2.0 * ajc * fic[3] + sqrt(5) * aj * fic[5]);
        Ef += -fi[5] * (5 * (mui - 4 * Uf[i]) * fic[5] + sqrt(5) * ajc * fic[4] + sqrt(6) * aj * fic[6]);
        Ef += -fi[6] * (6 * (mui - 5 * Uf[i]) * fic[6] + sqrt(6) * ajc * fic[5] + sqrt(7) * aj * fic[7]);
        Ef += -fi[7] * (7 * (mui - 6 * Uf[i]) * fic[7] + sqrt(7) * ajc * fic[6]);
        Ef += conj(aj) * af[i];
    }

    printMath(os, "EfUres", resi, Ef);
    Ef /= 2 * Ui;
    printMath(os, "Efres", resi, Ef);
    //        printMath(os, "nfres", resi, nf);
    
    complex<double> Q = Ef - E0;
    printMath(os, "Qres", resi, Q);

    double p = 0;
    complex<double> p2 = 1;
    for (int i = 0; i < L; i++) {
        fiarray f0i = f0[i];
        fiarray ffi = ff[i];
        complex<double> overlap = 0;
        for (int n = 0; n <= nmax; n++) {
            overlap += conj(ffi[n]) * f0i[n];
        }
        p += 1 - norm(overlap);
        p2 *= overlap;
    }
    p /= L;
    p2 = 1 - norm(p2);
    
    //        printMath(os, "ffres", resi, ff);
    printMath(os, "pres", resi, p);
    printMath(os, "p2res", resi, p2);

    deque<double> ttres;
    fpoints ftres;
    funcpoints normtres;
    funcpoints ntres;
    funcpoints n2tres;
    funcpoints Ftres;
    funcpoints atres;
    //        deque<pair<double, double> > Wtres;
    deque<double> Wtres;

    carrayL a;
    arrayL n;
    deque<double> fcres;
    deque<complex<double> > Phires;
    deque<double> Ftmeanres;
    //        int k = 0;
    for (int k = 0; k < tt.size(); k++) {
        ttres.push_back(tt[k]);
        for (int i = 0; i < L * dim; i++) {
            //                ftres[i].push_back(ft[k][i]);
        }
        complex<double> Phit = 0;
        double Ftmean = 0;
        for (int i = 0; i < L; i++) {
            fiarray fi = ft[k][i];
            arraydim finorms;
            for (int m = 0; m <= nmax; m++) {
                finorms[m] = norm(fi[m]);
            }
            double normt = 0;
            double nt = 0;
            double n2t = 0;
            for (int m = 0; m <= nmax; m++) {
                normt += finorms[m];
                nt += m * finorms[m];
                n2t += m * (m/* - 1*/) * finorms[m];
            }
            nt /= normt;
            n2t /= normt;

            double Ft = n2t - nt * nt;
            Ftmean += Ft;

            complex<double> at = ai(fi) / normt;
            normtres[i].push_back(normt);
            ntres[i].push_back(nt);
            n2tres[i].push_back(n2t);
            Ftres[i].push_back(Ft);
            atres[i].push_back(at);

            n[i] = nt;
            a[i] = at;
            Phit += at;
        }
        Ftmean /= L;
        Ftmeanres.push_back(Ftmean);
        Phit /= L;
        Phires.push_back(Phit);
        //            Wtres.push_back((*Wt)(tt[k]));
        Wtres.push_back(Wt(tt[k]));
    }
    printMath(os, "ttres", resi, ttres);
    printMath(os, "Ftmeanres", resi, Ftmeanres);
    printMath(os, "Phitres", resi, Phires);
    
//    printMath(os, "f0res", resi, f0);
//    printMath(os, "ffres", resi, ff);

    //        cout << math(ttres.back()) << endl << math(p) << endl << math(E0) << endl << math(Ef) << endl << math(Phires.back()) << endl;

    //    cout << Eigen::SimdInstructionSetsInUse() << endl;

    return 0;
}

