#include "piecewise_constant_rate_function.h"

constexpr long nC2(int n) { return n * (n - 1) / 2; }

template <typename T>
inline void vec_insert(std::vector<T> &v, const int pos, const T &x)
{
    v.insert(v.begin() + pos, x);
}

template <typename T>
inline T _conv(const adouble x);

template <typename T>
inline std::vector<T> _vconv(const std::vector<adouble> v);

template <>
inline double _conv(const adouble x) { return x.value(); }

template <>
inline adouble _conv(const adouble x) { return x; }

template <>
inline std::vector<adouble> _vconv(const std::vector<adouble> v) { return v; }

template <>
inline std::vector<double> _vconv(const std::vector<adouble> v) 
{ 
    std::vector<double> ret; 
    for (adouble x : v)
        ret.push_back(x.value());
    return ret;
}

template <typename T>
PiecewiseConstantRateFunction<T>::PiecewiseConstantRateFunction(
        const std::vector<std::vector<adouble>> params, 
        const std::vector<double> hidden_states) :
    params(params),
    nder(params[0][0].derivatives().size()),
    K(params[0].size()), ada(_vconv<T>(params[0])),
    s(_vconv<double>(params[1])),
    ts(K + 1), Rrng(K), 
    hidden_states(hidden_states),
    tmax(std::accumulate(s.begin(), s.end(), 0.0))
{
    for (auto &pp : params)
        if (pp.size() != params[0].size())
            throw std::runtime_error("all params must have same size");
    // Final piece is required to be flat.
    ts[0] = 0.;
    Rrng[0] = 0.;
    // These constant values need to have compatible derivative shape
    // with the calculated values.
    // Fix last piece to be constant
    for (int k = 0; k < K; ++k)
    {
        ada[k] = 1. / ada[k];
        ts[k + 1] = ts[k] + s[k];
    }
    ts[K] = INFINITY;

    int ip;
    for (T h : hidden_states)
    {
        if (std::isinf(toDouble(h)))
        {
            hs_indices.push_back(ts.size() - 1);
            continue;
        }
        ip = insertion_point(h, ts, 0, ts.size());
        if (myabs(ts[ip] - h) < 1e-8)
        // if (ts[ip] == h)
            hs_indices.push_back(ip);
        else if (ip + 1 < ts.size() and myabs(ts[ip + 1] - h) < 1e-8)
            hs_indices.push_back(ip + 1);
        else
        {
            vec_insert(ts, ip + 1, h);
            vec_insert<T>(ada, ip + 1, ada[ip]);
            check_nan(ada[ip + 1]);
            check_nan(ts[ip + 1]);
            hs_indices.push_back(ip + 1);
        }
    }
    K = ada.size();
    Rrng.resize(K + 1);
    compute_antiderivative();

    _eta.reset(new PExpEvaluator<T>(ada, ts, Rrng));
    _R.reset(new PExpIntegralEvaluator<T>(ada, ts, Rrng));
    _Rinv.reset(new PExpInverseIntegralEvaluator<T>(ada, ts, Rrng));
}

template <typename T>
inline T _double_integral_below_helper(const int rate, const T &tsm, const T &tsm1, const T &ada, 
        const T &Rrng, const T &log_denom)
{
    const int l1r = 1 + rate;
    T _tsm = tsm, _tsm1 = tsm1;
    T _ada = ada, _Rrng = Rrng; // don't ask
    T z = _tsm - _tsm;
    const T l1rinv = 1 / (z + l1r);
    T diff = _tsm1 - _tsm;
    T _adadiff = _ada * diff;
    if (rate == 0)
    {
        if (tsm1 == INFINITY)
            return exp(-_Rrng - log_denom) / _ada;
        else
            return exp(-_Rrng - log_denom) * (1 - exp(-_adadiff) * (1 + _adadiff)) / _ada;
    }
    if (tsm1 == INFINITY)
        return exp(-l1r * _Rrng - log_denom) * (1 - l1rinv) / (rate * _ada);
    return exp(-l1r * _Rrng - log_denom) * (expm1(-l1r * _adadiff) * l1rinv - expm1(-_adadiff)) / (rate * _ada);
}

template <typename U>
inline U _double_integral_above_helper(const int rate, const int lam, const U &_tsm, 
        const U &_tsm1, const U &_ada, const U &_Rrng, const U &log_coef)
{
    U diff = _tsm1 - _tsm;
    U adadiff = _ada * diff;
    long l1 = lam + 1;
    if (rate == 0)
        return exp(-l1 * _Rrng + log_coef) * (expm1(-l1 * adadiff) + l1 * adadiff) / l1 / l1 / _ada;
    if (l1 == rate)
    {
        if (_tsm1 == INFINITY)
            return exp(-rate * _Rrng + log_coef) / rate / rate / _ada;
        return exp(-rate * _Rrng + log_coef) * (1 - exp(-rate * adadiff) * (1 + rate * adadiff)) / rate / rate / _ada;
    }
    if (_tsm1 == INFINITY)
        return exp(-l1 * _Rrng + log_coef) / l1 / rate / _ada;
    // return -exp(-l1 * _Rrng + log_coef) * (expm1(-l1 * adadiff) / l1 + (exp(-rate * adadiff) - exp(-l1 * adadiff)) / (l1 - rate)) / rate / _ada;
    if (rate < l1)
        return -exp(-l1 * _Rrng + log_coef) * (expm1(-l1 * adadiff) / l1 + (exp(-rate * adadiff) * -expm1(-(l1 - rate) * adadiff) / (l1 - rate))) / rate / _ada;
    else
        return -exp(-l1 * _Rrng + log_coef) * (expm1(-l1 * adadiff) / l1 + (exp(-l1 * adadiff) * expm1(-(rate - l1) * adadiff) / (l1 - rate))) / rate / _ada;
}

template <typename T>
void PiecewiseConstantRateFunction<T>::print_debug() const
{
    std::vector<std::pair<std::string, std::vector<T>>> arys = 
    {{"ada", ada}, {"Rrng", Rrng}};
    std::cout << std::endl;
    for (auto p : arys)
    {
        std::cout << p.first << std::endl;
        for (adouble x : p.second)
            std::cout << x.value() << "::" << x.derivatives().transpose() << std::endl;
        std::cout << std::endl;
    }
}

template <typename T>
void PiecewiseConstantRateFunction<T>::compute_antiderivative()
{
    Rrng[0] = 0.;
    for (int k = 0; k < K; ++k)
        Rrng[k + 1] = Rrng[k] + ada[k] * (ts[k + 1] - ts[k]);
}

template <typename T>
T PiecewiseConstantRateFunction<T>::R_integral(const T a, const T b) const
{
    return R_integral(a, b, 0);
}

template <typename T>
T PiecewiseConstantRateFunction<T>::R_integral(const T a, const T b, const T log_denom) const
{
    // int_a^b exp(-R(t)) dt
    int ip_a = insertion_point(a, ts, 0, ts.size());
    int ip_b = (std::isinf(toDouble(b))) ? ts.size() - 2 : insertion_point(b, ts, 0, ts.size());
    T ret = 0., r, left, right, diff, Rleft;
    for (int i = ip_a; i < ip_b + 1; ++i)
    {
        left = dmax(a, ts[i]);
        right = dmin(b, ts[i + 1]);
        diff = right - left;
        Rleft = R(left);
        r = exp(-(Rleft + log_denom));
        if (!std::isinf(toDouble(diff)))
            r *= -expm1(-diff * ada[i]);
        r /= ada[i];
        check_negative(r);
        check_nan(r);
        ret += r;
    }
    return ret;
}


template <typename T>
inline T _single_integral(const int rate, const T &tsm, const T &tsm1, 
        const T &ada, const T &Rrng, const T &log_coef)
{
    // = int_ts[m]^ts[m+1] exp(-rate * R(t)) dt
    const int c = rate;
    if (rate == 0)
        return exp(log_coef) * (tsm1 - tsm);
    T ret = exp(-c * Rrng + log_coef);
    if (tsm1 < INFINITY)
        ret *= -expm1(-c * ada * (tsm1 - tsm));
    ret /= ada * c;
    check_negative(ret);
    return ret;
}

template <typename T>
void PiecewiseConstantRateFunction<T>::tjj_double_integral_above(const int n, long jj, std::vector<Matrix<T> > &C) const
{
    T tmp;
    long lam = nC2(jj) - 1;
    // Now calculate with hidden state integration limits
    for (unsigned int h = 0; h < hs_indices.size() - 1; ++h)
    {
        C[h].row(jj - 2).setZero();
        T log_denom = -Rrng[hs_indices[h]];
        if (Rrng[hs_indices[h + 1]] != INFINITY)
            log_denom += log(-expm1(-(Rrng[hs_indices[h + 1]] - Rrng[hs_indices[h]])));
        for (int m = hs_indices[h]; m < hs_indices[h + 1]; ++m)
        {
            for (int j = 2; j < n + 2; ++j)
            {
                long rate = nC2(j);
                tmp = _double_integral_above_helper<T>(rate, lam, ts[m], ts[m + 1], ada[m], Rrng[m], -log_denom);
                check_nan(tmp);
                check_negative(tmp);
                try 
                {
                    check_nan(C[h](jj - 2, j - 2));
                    check_negative(C[h](jj - 2, j - 2));
                } catch (std::runtime_error)
                {
                    CRITICAL << m << " " << rate << " " << lam << " " << ts[m] << " " << ts[m + 1] << 
                        " " << ada[m] << " " << Rrng[m] << " " << log_denom;
                    throw;
                }
                C[h](jj - 2, j - 2) += tmp;
                T log_coef = -log_denom, fac;
                long rp = lam + 1 - rate;
                if (rp == 0)
                    fac = Rrng[m + 1] - Rrng[m];
                else
                {
                    if (rp < 0)
                    {
                        if (-rp * (Rrng[m + 1] - Rrng[m]) > 20)
                        {
                            log_coef += -rp * Rrng[m + 1];
                            fac = -1. / rp;
                        }
                        else
                        {
                            log_coef += -rp * Rrng[m];
                            fac = -expm1(-rp * (Rrng[m + 1] - Rrng[m])) / rp;
                        }
                    }
                    else
                    {
                        if (-rp * (Rrng[m] - Rrng[m + 1]) > 20)
                        {
                            log_coef += -rp * Rrng[m];
                            fac = 1. / rp;
                        }
                        else
                        {
                            log_coef += -rp * Rrng[m + 1];
                            fac = expm1(-rp * (Rrng[m] - Rrng[m + 1])) / rp;
                        }
                    }
                }
                for (int k = m + 1; k < K; ++k)
                {
                    T si = _single_integral(rate, ts[k], ts[k + 1], ada[k], Rrng[k], log_coef) * fac;
                    C[h](jj - 2, j - 2) += si;
                    check_nan(C[h](jj - 2, j - 2));
                    check_negative(C[h](jj - 2, j - 2));
                }
                check_nan(C[h](jj - 2, j - 2));
                check_negative(C[h](jj - 2, j - 2));
            }
        }
    }
}

template <typename T>
void PiecewiseConstantRateFunction<T>::tjj_double_integral_below(
        const int n, const int h, Matrix<T> &tgt) const
{
    T log_denom = -Rrng[hs_indices[h]];
    if (Rrng[hs_indices[h + 1]] != INFINITY)
        log_denom += log(-expm1(-(Rrng[hs_indices[h + 1]] - Rrng[hs_indices[h]])));
    for (int m = hs_indices[h]; m < hs_indices[h + 1]; ++m)
    {
        Vector<T> ts_integrals(n + 1);
        T log_coef = -Rrng[m];
        T fac = 1.;
        if (m < K - 1)
            fac = -expm1(-(Rrng[m + 1] - Rrng[m]));
        for (int j = 2; j < n + 3; ++j)
        {
            long rate = nC2(j) - 1;
            ts_integrals(j - 2) = _double_integral_below_helper<T>(rate, ts[m], ts[m + 1], ada[m], Rrng[m], -log_denom);
            for (int k = 0; k < m; ++k)
            {
                T _c = log_coef - log_denom;
                ts_integrals(j - 2) += fac * _single_integral(rate, ts[k], ts[k + 1], ada[k], Rrng[k], _c);
            }
            check_negative(ts_integrals(j - 2));
        }
        tgt.row(h) += ts_integrals.transpose();
    }
}

template <typename T>
T exp1_conditional(T a, T b, std::mt19937 &gen)
{
    // If X ~ Exp(1),
    // P(X < x | a <= X <= b) = (e^-a - e^-x) / (e^-a - e^-b)
    // so P^-1(y) = -log(e^-a - (e^-a - e^-b) * y)
    //            = -log(e^-a(1 - (1 - e^-(b-a)) * y)
    //            = a - log(1 - (1 - e^-(b-a)) * y)
    //            = a - log(1 + expm1(-(b-a)) * y)
    double unif = std::uniform_real_distribution<double>{0.0, 1.0}(gen);
    if (std::isinf(toDouble(b)))
        return a - log1p(-unif);
    else
        return a - log1p(expm1(-(b - a)) * unif);
}

// This helper function exists for cython
template <typename T>
double PiecewiseConstantRateFunction<T>::random_time(const double a, const double b, const long long seed) const
{
    std::mt19937 gen(seed);
    return toDouble(random_time(1., T(a), T(b), gen));
}


template <typename T>
T PiecewiseConstantRateFunction<T>::random_time(const double fac, const T &a, const T &b, std::mt19937 &gen) const
{
    T Rb;
    if (b == INFINITY)
        Rb = INFINITY;
    else
        Rb = R(b);
    return (*getRinv())(exp1_conditional(R(a), Rb, gen) / fac);
}


template <typename T>
std::vector<T> PiecewiseConstantRateFunction<T>::average_coal_times() const
{
    std::vector<T> ret;
    for (int i = 1; i < hidden_states.size(); ++i)
    {
        // discretize by expected coalescent time within each hidden state
        // e_coal = \int_t0^t1 t eta(t) exp(-R(t)) 
        //        = t0 e^(-R(t0)) - t1 e^(-R(t1)) + \int
        T log_denom = -Rrng[hs_indices[i - 1]];
        bool inf = std::isinf(toDouble(ts[hs_indices[i]]));
        if (!inf)
           log_denom += log(-expm1(-(Rrng[hs_indices[i]] - Rrng[hs_indices[i - 1]])));
        T x = hidden_states[i - 1] * exp(-((Rrng[hs_indices[i - 1]]) + log_denom)) +
            R_integral(ts[hs_indices[i - 1]], ts[hs_indices[i]], log_denom);
        if (!inf)
            x -= hidden_states[i] * exp(-((Rrng[hs_indices[i]]) + log_denom));
        ret.push_back(x);
        if (ret.back() > hidden_states[i] or ret.back() < hidden_states[i - 1])
            throw std::runtime_error("erroneous average coalescence time");
    }
    return ret;
}


template class PiecewiseConstantRateFunction<double>;
template class PiecewiseConstantRateFunction<adouble>;