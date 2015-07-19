#include "conditioned_sfs.h"
#include "gauss_legendre.h"

std::mt19937 sfs_gen;

template <typename T>
struct helper_struct
{
    const FunctionEvaluator<T>* f;
    double rate;
    T offset;
};

template <typename T>
T helper(T x, void* obj)
{
    helper_struct<T>* args = (helper_struct<T>*)obj;
    return exp(-((*args->f)(x) - args->offset) * args->rate);
}

std::map<std::array<int, 3>, double> _Wnbj_memo;
double calculate_Wnbj(int n, int b, int j)
{
    switch (j)
    {
        case 2:
            return (double)6 / (n + 1);
        case 3:
            return (double)30 * (n - 2 * b) / (n + 1) / (n + 2);
        default:
            std::array<int, 3> key = {n, b, j};
            if (_Wnbj_memo.count(key) == 0)
            {
                int jj = j - 2;
                double ret = calculate_Wnbj(n, b, jj) * -(1 + jj) * (3 + 2 * jj) * (n - jj) / jj / (2 * jj - 1) / (n + jj + 1);
                ret += calculate_Wnbj(n, b, jj + 1) * (3 + 2 * jj) * (n - 2 * b) / jj / (n + jj + 1);
                _Wnbj_memo[key] = ret;
            }
            return _Wnbj_memo[key];
    }
}

std::map<std::array<int, 2>, long> _binom_memo;
long binom(int n, int k)
{
    assert(k >= 0);
    if (k == 0 or n == k)
        return 1;
    std::array<int, 2> key = {n, k};
    if (_binom_memo.count(key) == 0) 
        _binom_memo[key] = binom(n - 1, k - 1) + binom(n - 1, k);
    return _binom_memo[key];
}

double pnkb_dist(int n, int m, int l1)
{
    // Probability that lineage 1 has size |L_1|=l1 below tau,
    // the time at which 1 and 2 coalesce, when there are k 
    // undistinguished lineages remaining, in a sample of n
    // undistinguished (+2 distinguished) lineages overall.
    double ret = l1 * (n + 2 - l1) / (double) binom(n + 3, m + 3);
    if (m > 0)
        ret *= (n + 1 - l1) * (double)binom(n - l1, m - 1) / m / (m + 1);
    return ret;
}

double pnkb_undist(int n, int m, int l3)
{
    // Probability that undistinguished lineage has size |L_1|=l1 below tau,
    // the time at which 1 and 2 coalesce, when there are k 
    // undistinguished lineages remaining, in a sample of n
    // undistinguished (+2 distinguished) lineages overall.
    assert(m > 0);
    double ret = (n + 3 - l3) * (n + 2 - l3) * (n + 1 - l3) / (double) binom(n + 3, m + 3);
    if (m == 1)
        ret /= 6.0;
    else
        ret *= (n - l3) * (double)binom(n - l3 - 1, m - 2) / (m - 1) / m / (m + 1) / (m + 2);
    return ret;
}


template <typename T>
ConditionedSFS<T>::ConditionedSFS(const PiecewiseExponentialRateFunction<T> eta, int n, 
        MatrixInterpolator moran_interp) : 
    eta(eta), n(n), moran_interp(moran_interp),
    D_subtend_above(n, n), D_not_subtend_above(n, n),
	D_subtend_below(n + 1, n + 1), D_not_subtend_below(n + 1, n + 1),
	Wnbj(n, n), P_dist(n + 1, n + 1), 
    P_undist(n + 1, n), tK(n + 1, n + 1), csfs(3, n + 1), 
    csfs_above(3, n + 1), csfs_below(3, n + 1), ETnk_below(n + 1, n + 1)
{
    long long seed = std::uniform_int_distribution<long long>{}(sfs_gen);
    gen.seed(seed);
    fill_matrices();
}

template <typename T>
double ConditionedSFS<T>::rand_exp(void)
{
    return std::exponential_distribution<double>{1.0}(gen);
}

template <typename T>
double ConditionedSFS<T>::unif(void)
{
    return std::uniform_real_distribution<double>{0.0, 1.0}(gen);
}

template <typename T>
T ConditionedSFS<T>::exp1_conditional(T a, T b)
{
    // If X ~ Exp(1),
    // P(X < x | a <= X <= b) = (e^-a - e^-x) / (e^-a - e^-b)
    // so P^-1(y) = -log(e^-a - (e^-a - e^-b) * y)
    //            = -log(e^-a(1 - (1 - e^-(b-a)) * y)
    //            = a - log(1 - (1 - e^-(b-a)) * y)
    //            = a - log(1 + expm1(-(b-a)) * y)
    if (std::isinf(toDouble(b)))
        return a - log1p(-unif());
    else
        return a - log1p(expm1(-(b - a)) * unif());
}

template <typename T>
void ConditionedSFS<T>::compute_ETnk_below(const Vector<T> &etjj)
{
    ETnk_below.setZero();
    ETnk_below.diagonal() = etjj;
    for (int nn = 2; nn < n + 3; ++nn)
    {
        for (int k = nn - 1; k > 1; --k)
        {
            ETnk_below(nn - 2, k - 2) = ETnk_below(nn - 3, k - 2) -
                (double)(k + 2) * (k - 1) / (nn + 1) / (nn - 2) * ETnk_below(nn - 2, k - 1);
            ETnk_below(nn - 2, k - 2) /= 1.0 - (double)(k + 1) * (k - 2) / (nn + 1) / (nn - 2);
        }
    }
}

template <typename T>
T etnk_recursion(std::map<std::pair<int, int>, T> &memo, const Vector<T> &etjj, int n, int k) 
{
    if (k == n)
        return etjj(k - 2);
    std::pair<int, int> key = {n, k};
    if (memo.count(key) == 0)
    {
        T ret = etnk_recursion(memo, etjj, n - 1, k);
        ret -= (double)(k + 2) * (k - 1) / (n + 1) / (n - 2) * etnk_recursion(memo, etjj, n, k + 1);
        ret /= 1.0 - (double)(k + 1) * (k - 2) / (n + 1) / (n - 2);
        memo[key] = ret;
    }
    return memo[key];
}

template <typename T>
std::thread ConditionedSFS<T>::compute_threaded(int num_samples, T t1, T t2)
{
    return std::thread(&ConditionedSFS::compute, this, num_samples, t1, t2);    
}

constexpr int nc2(int n) { return n * (n - 1) / 2; }

template <typename T>
void ConditionedSFS<T>::compute(int num_samples, T t1, T t2)
{
    // feenableexcept(FE_INVALID | FE_OVERFLOW);
    auto R = eta.getR();
    auto Rinv = eta.getRinv();
    auto feta = eta.geteta();
    int m;
    T y, tau, rate, tp;
	// Mixing constants with adoubles causes problems because the library
	// doesn't know how to allocate the VectorXd of derivatives(). 
	// Here, we do it by hand.
	T zero = eta.zero;
    csfs.fill(zero);
    csfs_above.fill(zero);
    csfs_below.fill(zero);
    std::vector<T> ys(num_samples);
    std::generate(ys.begin(), ys.end(), [&](){return exp1_conditional(t1, t2);});
    std::sort(ys.begin(), ys.end());
    std::vector<int> eis(num_samples);
    std::vector<T> taus = (*Rinv)(ys);
    for (m = 0; m < num_samples; ++m)
    {
        tau = taus[m];
        y = ys[m];
        csfs(1,0) += 2.0 * tau;
        int nsingletons = n;
        T Rt = 0.0; 
        T t = 0.0;
        int coal12 = 1;
        for (int k = n + 2; k > 1; --k)
        {
            if (nsingletons == 0)
                break;
            rate = k * (k - 1) / 2 - coal12;
            Rt += rand_exp() / rate;
            tp = (*Rinv)(Rt);
            if (tp > tau and coal12)
            {
                // Distinguished coalescence happens
                coal12 = 0;
                tp = tau;
                Rt = y;
                csfs(0,1) += nsingletons * (tau - t);
            }
            else
            {
                csfs(0,1) += nsingletons * (tp - t);
                // 3 possibilites:
                // coalescence between singletons
                // coalescence not involving singleton
                // coalescenc singleton <-> non-singleton
                double pcoal1 = nc2(nsingletons);
                double pcoal2 = nc2(k - nsingletons) - coal12;
                double pcoal3 = nc2(k) - coal12 - pcoal1 - pcoal2;
                int coal_type = std::discrete_distribution<int>{pcoal1, pcoal2, pcoal3}(gen);
                switch (coal_type)
                {
                    case 0:
                        nsingletons -= 2; break;
                    case 2:
                        nsingletons -= 1; break;
                }
            }
            t = tp;
        }
    }
    csfs /= num_samples;
    // std::cout << csfs.template cast<double>() << std::endl << std::endl;
}


template <typename T>
void ConditionedSFS<T>::fill_matrices(void)
{
	Matrix<T> I = Matrix<T>::Identity(n, n);
	Matrix<T> I1 = Matrix<T>::Identity(n + 1, n + 1);
    // Construct some matrices that will be used later on
    D_subtend_above.setZero();
    D_subtend_above.diagonal() = Eigen::VectorXd::LinSpaced(n, 1, n).template cast<T>() / (n + 1);
	D_not_subtend_above = I - D_subtend_above;

    D_subtend_below.setZero();
    for (int k = 2; k < n + 3; ++k)
        D_subtend_below.diagonal()(k - 2) = 2. / k;
	D_not_subtend_below = I1 - D_subtend_below;

    tK.setZero();
    tK.diagonal() = Eigen::VectorXd::LinSpaced(n + 1, 2, n + 2).template cast<T>();

    // Calculate the Polanski-Kimmel matrix
    // TODO: this could be sped up by storing the matrix outside of the class
    Wnbj.setZero();
    for (int b = 1; b < n + 1; ++b)
        for (int j = 2; j < n + 2; ++j)
            Wnbj(b - 1, j - 2) = calculate_Wnbj(n + 1, b, j);

    // P_dist(k, b) = probability of state (1, b) when there are k undistinguished lineages remaining
    P_dist.setZero();
    for (int k = 0; k < n + 1; ++k)
        for (int b = 1; b < n - k + 2; ++b)
            P_dist(k, b - 1) = pnkb_dist(n, k, b);

    // P_undist(k, b) = probability of state (0, b + 1) when there are k undistinguished lineages remaining
    P_undist.setZero();
    for (int k = 1; k < n + 1; ++k)
        for (int b = 1; b < n - k + 2; ++b)
            P_undist(k, b - 1) = pnkb_undist(n, k, b);
}

void print_sfs(int n, const std::vector<double> &sfs)
{
    std::vector<double> rsfs(n, 0);
    double x;
    int k = 0;
    for (int i = 0; i < 3; i++)
    {
        printf("%i:\t", i);
        for (int j = 0; j < n; j++)
        {
            x = sfs[k++];
            rsfs[i + j] += x;
            printf("%i:%e ", j, x);
        }
        printf("\n");
    }
    for (int i = 0; i < n; i++)
    {
        printf("%i:%f\n", i, rsfs[i]);
    }
}

template <typename T>
Matrix<T> ConditionedSFS<T>::average_csfs(std::vector<ConditionedSFS<T>> &csfs, double theta)
{
    Matrix<T> ret = Matrix<T>::Zero(csfs[0].matrix().rows(), csfs[0].matrix().cols());
    int m = 0;
    for (const ConditionedSFS<T> &c : csfs)
    {
        ret += c.matrix();
        ++m;
    }
    ret /= (double)m;
    ret *= theta;
    T tauh = ret.sum();
    ret(0, 0) = 1.0 - tauh;
    // ret *= -expm1(-theta * tauh) / tauh;
    // ret(0, 0) = exp(-theta * tauh);
    /*
    T undist = ret(0, 1);
    ret.row(0).fill(ret(0, 0));
    ret(0, 1) = undist;
    ret.row(2).fill(ret(0, 0));
    ret.row(1).fill(ret(1, 0));
    */
    // ret *= theta;
    // ret(0, 0) = 1. - ret.sum();
    return ret;
}

template <typename T>
Matrix<T> ConditionedSFS<T>::calculate_sfs(const PiecewiseExponentialRateFunction<T> &eta, int n, int num_samples, 
        const MatrixInterpolator &moran_interp, double tau1, double tau2, int numthreads, double theta)
{
    // eta.print_debug();
    auto R = eta.getR();
    std::vector<ConditionedSFS<T>> csfs;
    std::vector<std::thread> t;
    T t1 = (*R)(tau1);
    T t2;
    if (std::isinf(tau2))
        t2 = INFINITY;
    else
        t2 = (*R)(tau2);
    std::vector<Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>> _expM;
    if (numthreads == 0)
    {
        csfs.emplace_back(eta, n, moran_interp);
        csfs.back().compute(num_samples, t1, t2);
    }
    else
    {
        for (int i = 0; i < numthreads; ++i)
            csfs.emplace_back(eta, n, moran_interp);
        for (auto &c : csfs)
            t.push_back(c.compute_threaded(num_samples, t1, t2));
        std::for_each(t.begin(), t.end(), [](std::thread &t) {t.join();});
    }
    Eigen::Matrix<T, 3, Eigen::Dynamic> ret = ConditionedSFS<T>::average_csfs(csfs, theta);
    if (ret(0,0) <= 0.0 or ret(0.0) >= 1.0)
    {
        std::cout << ret.template cast<double>() << std::endl << std::endl;
        std::cout << t1 << " " << t2 << std::endl << std::endl;
        std::cerr << "sfs is no longer a probability distribution. branch lengths are too long." << std::endl;
    }
    return ret;
}

void set_seed(long long seed)
{
    sfs_gen.seed(seed);
}

void store_sfs_results(const Matrix<double> &csfs, double* outsfs)
{
    int n = csfs.cols();
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> _outsfs(outsfs, 3, n);
    _outsfs = csfs.cast<double>();
}

void store_sfs_results(const Matrix<adouble> &csfs, double* outsfs, double* outjac)
{
    store_sfs_results(csfs.cast<double>(), outsfs);
    int n = csfs.cols();
    int num_derivatives = csfs(0,1).derivatives().rows();
    Eigen::VectorXd d;
    int m = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < n; ++j)
        {
            d = csfs(i, j).derivatives();
            assert(d.rows() == num_derivatives);
            for (int k = 0; k < num_derivatives; ++k)
                outjac[m++] = d(k);
        }
}

void cython_calculate_sfs(const std::vector<std::vector<double>> params,
        int n, int num_samples, const MatrixInterpolator &moran_interp,
        double tau1, double tau2, int numthreads, double theta, 
        double* outsfs)
{
    RATE_FUNCTION<double> eta(params);
    // eta.print_debug();
    Matrix<double> out = ConditionedSFS<double>::calculate_sfs(eta, n, num_samples, moran_interp, tau1, tau2, numthreads, theta);
    store_sfs_results(out, outsfs);
}

void cython_calculate_sfs_jac(const std::vector<std::vector<double>> params,
        int n, int num_samples, const MatrixInterpolator &moran_interp,
        double tau1, double tau2, int numthreads, double theta, 
        double* outsfs, double* outjac)
{
    RATE_FUNCTION<adouble> eta(params);
    // eta.print_debug();
    Matrix<adouble> out = ConditionedSFS<adouble>::calculate_sfs(eta, n, num_samples, moran_interp, tau1, tau2, numthreads, theta);
    store_sfs_results(out, outsfs, outjac);
}

template class ConditionedSFS<double>;
template class ConditionedSFS<adouble>;
