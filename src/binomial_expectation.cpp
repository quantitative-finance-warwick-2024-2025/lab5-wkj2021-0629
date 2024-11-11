#include <cmath>
#include <vector>

double binomial_expectation_loop(double *spread, double p, unsigned int n) {
    std::vector<double> arr(n + 1, 0.0);

    // set up the base nodes
    // loop through each depth
    // loop through each node at depth i

    double answer = arr[0];
    return answer;
}

double binomial_expectation_recursion_h(
    unsigned int i, unsigned int j, double *spread, double p, unsigned int n) {
    // if (i == n)
    // if (i < n)
    return 0.0;
}

double binomial_expectation_recursion(double *spread, double p, unsigned int n) {
    return binomial_expectation_recursion_h(0, 0, spread, p, n);
}

double risk_neutral_probability(double volatility, double delta_t, double rate_of_return) {
    double u = std::exp(volatility * std::sqrt(delta_t));
    double d = 1.0 / u;
    return (std::exp(rate_of_return) - d) / (u - d);
}

double expected_returns(double S0,
                        unsigned int number_of_periods,
                        double volatility,
                        double time_period,
                        double continuous_returns_rate) {
    double delta_t = time_period / (double)number_of_periods;
    double u = std::exp(volatility * std::sqrt(delta_t));
    double d = 1.0 / u;
    double p = (continuous_returns_rate - d) / (u - d);

    std::vector<double> spread(number_of_periods + 1);

    for (unsigned int i = 0; i < number_of_periods + 1; ++i) {
        spread[i] = S0 * std::pow(u, 2.0 * (double)i - (double)number_of_periods);
    }

    return binomial_expectation_loop(spread.data(), p, number_of_periods);
}

void expected_returns(double S0,
                      unsigned int number_of_periods,
                      double volatility,
                      double time_period,
                      size_t spread_points,
                      double *out_spread,
                      double *out_risk_neutral_rate,
                      double *out_probability,
                      double *out_returns) {

    double delta_t = time_period / (double)number_of_periods;
    double u = std::exp(volatility * std::sqrt(delta_t));
    double d = 1.0 / u;

    // populate the output arrays
}
