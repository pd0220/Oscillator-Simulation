// oscillator Monte Carlo simulation with given temperature

// including used headers and libraries
#include <iostream>
#include <math.h>
#include <numeric>
#include <random>

// lambda to calculate squared values
auto sq = [](auto const a) { return a * a; };

// calculate beta * delta energy
// val --> beta * k * a^2
template <typename T>
auto BetaDeltaE(T val, int nPrev, int nNext)
{
    using R = decltype(val * nPrev);
    return (R)0.5 * val * (sq(nNext) - sq(nPrev));
}

// calculate rate
template <typename T>
auto Rate(T val, int nPrev, int nNext)
{
    auto exponent = BetaDeltaE(val, nPrev, nNext);
    // to lesser energy level
    if (exponent < 0)
        return 1;
    // to higher energy level
    else if (exponent > 0)
        return std::exp(-exponent);
    // error
    else
    {
        std::cout << "ERROR\nWrong state numbers are given." << std::endl;
        std::exit(-1);
    }
}

// choose direction (left or right)
bool Direction()
{
    std::random_device rd{};
    std::mt19937 gen(rd());
    // uniform distribution in closed interval [0,1]
    std::uniform_int_distribution<> distr(0, 1);
    return distr(gen);
}

// main function
int main(int, char **)
{

    std::cout << Direction() << std::endl;
}
